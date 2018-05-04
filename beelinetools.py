#!/usr/bin/env python
"""Script that converts the beeline report to one that is readable by Plink."""


from __future__ import print_function

import os
import re
import sys
import logging
import argparse
from collections import namedtuple
from tempfile import NamedTemporaryFile

try:
    from pyplink import PyPlink
    _HAS_PYPLINK = True
except ImportError:
    _HAS_PYPLINK = False

try:
    from six.moves import range, zip
except ImportError:
    pass


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2015, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "MIT"
__version__ = "0.3.0"


# The location tuple
_Location = namedtuple("Location", ["chrom", "pos", "alleles"])
_unknown_location = _Location(chrom=0, pos=0, alleles={})

# The complement (including special alleles)
_complement = {"A": "T", "T": "A", "C": "G", "G": "C",
               "-": "-", "D": "D", "I": "I"}

# The AB encoding
_geno_add_encoding = {"AA": 0, "AB": 1, "BA": 1, "BB": 2}

# The forward/reverse strand formatter
_forward_strand_re = re.compile(r"_[TBPM]_([RF])_")


def main():
    """The main function."""
    # Setting the logging
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    # The parser
    desc = ("Performs different tasks on Illumina's beeline report(s) "
            "(version {}).".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    try:
        # Parsing and checking the options and arguments
        args = parse_args(parser)
        logging.info("Version {}".format(__version__))
        check_args(args)

        # Getting the strand column
        strand_column = None
        if args.beeline_strand == "forward":
            strand_column = args.map_illumina_id
        elif args.beeline_strand == "plus":
            strand_column = args.map_ref_strand
        elif args.beeline_strand == "top":
            strand_column = args.map_illumina_strand

        # Reading the map file
        map_data = read_mapping_info(
            i_filename=args.map_filename,
            delim=args.map_delim,
            map_id=args.map_id,
            map_chr=args.map_chr,
            map_pos=args.map_pos,
            map_allele=args.map_allele,
            allele_strand=args.beeline_strand,
            strand_col=strand_column,
        )

        if args.analysis_type == "convert":
            # Converting the beeline report(s)
            convert_beeline(
                i_filenames=args.i_filenames,
                out_dir=args.output_dir,
                locations=map_data,
                other_opts=args,
            )

        elif args.analysis_type == "extract":
            # Reading the list of samples to keep
            samples_to_keep = read_list(args.samples_to_keep)

            # Extracting from the beeline report(s)
            extract_beeline(
                i_filenames=args.i_filenames,
                out_dir=args.output_dir,
                o_suffix=args.o_suffix,
                locations=map_data,
                samples=samples_to_keep,
                other_opts=args,
            )

        elif args.analysis_type == "sample-split":
            # Split the report per sample
            split_report(
                i_filenames=args.i_filenames,
                out_dir=args.output_dir,
                locations=map_data,
                other_opts=args,
            )

    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        sys.exit(1)

    except ProgramError as e:
        logging.error(e)
        parser.error(e.message)

    except Exception as e:
        logging.error(e)
        raise


def convert_beeline(i_filenames, out_dir, locations, other_opts):
    """Convert beeline report(s) to Plink files.

    Args:
        i_filenames (list): a list of file names (str)
        out_dir (str): the name of the output directory
        locations (dict): a dictionary from marker ID to genomic location
        other_opts(argparse.Namespace): the program options

    """
    # The samples that were already seen
    seen_samples = set()

    for i_filename in i_filenames:
        logging.info("Converting '{}'".format(i_filename))

        # Getting the output filename
        o_filename = os.path.splitext(os.path.basename(i_filename))[0]

        # Opening the file
        i_file = None
        if i_filename != "-":
            i_file = open(i_filename, "r")

        else:
            i_file = sys.stdin
            o_filename = "from_stdin"

        # The output files
        prefix = os.path.join(out_dir, o_filename)
        pedfile = None
        mapfile = None
        bedfile = None
        bimfile = None
        famfile = None
        sample_file = None
        sample_file_end = None
        sample_file_sep = None
        if other_opts.o_format == "ped":
            pedfile = open(prefix + ".ped", "w")
            mapfile = open(prefix + ".map", "w")
            sample_file = pedfile
            sample_file_end = ""
            sample_file_sep = "\t"
        elif other_opts.o_format == "bed":
            bedfile = PyPlink(prefix, "w", "INDIVIDUAL-major")
            bimfile = open(prefix + ".bim", "w")
            famfile = open(prefix + ".fam", "w")
            sample_file = famfile
            sample_file_end = "\n"
            sample_file_sep = " "

        # Reading the file (or STDIN)
        try:
            # The number of markers
            nb_markers = None

            # Reading the assay information
            line = i_file.readline().rstrip("\r\n")
            while line != "[Header]":
                line = i_file.readline().rstrip("\r\n")

            while not line.startswith("[Data]"):
                if line.startswith(other_opts.nb_snps_kw):
                    nb_markers = int(line.rstrip("\r\n").split(",")[-1])
                line = i_file.readline()

            if nb_markers is None:
                raise ProgramError(
                    "{}: invalid header (missing '{}' value)".format(
                        i_filename,
                        other_opts.nb_snps_kw,
                    )
                )

            logging.info("There are {:,d} markers".format(nb_markers))

            # Reading and checking the header
            header = {
                name: i
                for i, name in
                enumerate(i_file.readline().rstrip("\r\n").split(","))
            }
            required_columns = (
                other_opts.beeline_id, other_opts.beeline_sample,
                other_opts.beeline_a1, other_opts.beeline_a2,
            )
            for name in required_columns:
                if name not in header:
                    raise ProgramError(
                        "{}: '{}': missing column".format(i_filename, name),
                    )

            # Creating the list that will contain the marker list
            current_marker_i = 0
            all_markers = [None for i in range(nb_markers)]
            nb_samples = 0

            # The genotypes (if in BED format)
            genotypes = None
            if other_opts.o_format == "bed":
                genotypes = [-1 for i in range(nb_markers)]

            # Reading the first data line
            line = i_file.readline()
            row = line.rstrip("\r\n").split(",")

            # Allele of each markers
            marker_alleles = {}

            while line != "":
                # Getting the marker name and sample id
                sample = row[header[other_opts.beeline_sample]]
                if sample in seen_samples:
                    logging.warning("{}: duplicate sample "
                                    "found".format(sample))
                seen_samples.add(sample)

                # Logging
                logging.info("Processing {}".format(sample))
                print(sample, sample, "0", "0", "0", "-9", sep=sample_file_sep,
                      end=sample_file_end, file=sample_file)

                # Reading the rest of the data for this sample
                current_sample = sample
                while current_sample == sample:
                    # Checking the marker order
                    marker = row[header[other_opts.beeline_id]]

                    # If the index is > than the length, it might be a
                    # duplicated sample...
                    if current_marker_i == len(all_markers):
                        break

                    if all_markers[current_marker_i] is None:
                        all_markers[current_marker_i] = marker
                    if all_markers[current_marker_i] != marker:
                        raise ProgramError(
                            "{}: marker order is not the same for "
                            "sample '{}'".format(i_filename, sample)
                        )

                    # Getting the genotype
                    allele_1 = row[header[other_opts.beeline_a1]]
                    allele_2 = row[header[other_opts.beeline_a2]]
                    genotype = "{} {}".format(allele_1, allele_2)
                    if "-" in genotype:
                        genotype = "0 0"

                    if other_opts.o_format == "ped":
                        pedfile.write("\t" + genotype)
                    else:
                        if marker not in locations:
                            raise ProgramError(
                                "{}: no mapping information".format(marker)
                            )

                        # Is this the first time we have seen this marker?
                        if marker not in marker_alleles:
                            marker_alleles[marker] = locations[marker].alleles

                        # Computing the genotypes
                        genotypes[current_marker_i] = encode_genotype(
                            allele_1, allele_2, marker_alleles[marker],
                        )

                    # Increasing the current marker
                    current_marker_i += 1

                    # Reading the next row
                    line = i_file.readline()
                    if line == "":
                        # End of file
                        break

                    # Splitting and current sample
                    row = line.rstrip("\r\n").split(",")
                    current_sample = row[header[other_opts.beeline_sample]]

                if other_opts.o_format == "ped":
                    pedfile.write("\n")
                else:
                    bedfile.write_genotypes(genotypes)

                # If there is only one marker, there is a problem
                if nb_markers != 1 and current_marker_i == 1:
                    raise ProgramError(
                        "{}: data should be sorted by samples, not by "
                        "markers ('{}' had 1 marker, expecting "
                        "{:,d}".format(i_filename, sample, nb_markers)
                    )

                # Are there any missing marker?
                if current_marker_i != nb_markers:
                    nb_missing = abs(current_marker_i - nb_markers)
                    raise ProgramError(
                        "{}: missing {} marker{} for sample '{}'".format(
                            i_filename,
                            nb_missing,
                            "s" if nb_missing > 1 else "",
                            sample,
                        )
                    )
                current_marker_i = 0
                nb_samples += 1

            # Closing the output files
            logging.info("Done writing {:,d} samples".format(nb_samples))

            # Printing the map file
            for marker in all_markers:
                marker_location = None
                if marker in locations:
                    marker_location = locations[marker]
                else:
                    marker_location = _unknown_location
                    logging.warning("{}: no mapping "
                                    "information".format(marker))
                if other_opts.o_format == "ped":
                    print(marker_location.chrom, marker, "0",
                          marker_location.pos, sep="\t", file=mapfile)
                else:
                    alleles = {
                        v: k for k, v in marker_alleles[marker].items()
                    }
                    print(marker_location.chrom, marker, "0",
                          marker_location.pos, alleles["B"], alleles["A"],
                          sep="\t", file=bimfile)

        finally:
            # Closing the input file
            i_file.close()

            # Closing the output files
            for f in (pedfile, mapfile, bedfile, bimfile, famfile):
                if f is not None:
                    f.close()


def encode_genotype(a1, a2, encoding):
    """Encode a genotype according to its A/B alleles.

    Args:
        geno (str): the genotype to encode
        encoding (dict): the genotype encoding

    Returns:
        int: the encoded genotype

    Note
    ----
        The encoding is as follow, 0 for the AA genotype, 1 for the AB
        genotype, and 2 for the BB genotype. Unknown genotype is encoded as -1.

    """
    return _geno_add_encoding.get(
        encode_allele(a1, encoding) + encode_allele(a2, encoding),
        -1,
    )


def split_report(i_filenames, out_dir, locations, other_opts):
    """Split beeline report(s) into a file for each sample.

    Args:
        i_filenames (list): a list of file names (str)
        out_dir (str): the name of the output directory
        locations (dict): a dictionary from marker ID to genomic location
        other_opts(argparse.Namespace): the program options

    """
    # The samples that were already seen
    seen_samples = set()

    for i_filename in i_filenames:
        logging.info("Splitting '{}'".format(i_filename))

        # Opening the file
        i_file = None
        if i_filename != "-":
            i_file = open(i_filename, "r")
        else:
            i_file = sys.stdin

        # Reading the file (or STDIN)
        try:
            # The number of markers
            nb_markers = None

            # Reading the assay information
            line = i_file.readline()
            while line.rstrip("\r\n") != "[Header]":
                line = i_file.readline()

            # The metadata
            metadata = line

            while not line.startswith("[Data]"):
                if line.startswith(other_opts.nb_snps_kw):
                    nb_markers = int(line.rstrip("\r\n").split(",")[-1])
                line = i_file.readline()
                metadata += line

            if nb_markers is None:
                raise ProgramError(
                    "{}: invalid header (missing '{}' value)".format(
                        i_filename,
                        other_opts.nb_snps_kw,
                    )
                )

            logging.info("There are {:,d} markers".format(nb_markers))

            # Reading and checking the header
            header_row = i_file.readline().rstrip("\r\n").split(",")
            header = {name: i for i, name in enumerate(header_row)}
            required_columns = (
                other_opts.beeline_id, other_opts.beeline_sample,
                other_opts.beeline_a1, other_opts.beeline_a2)
            for name in required_columns:
                if name not in header:
                    raise ProgramError(
                        "{}: '{}': missing column".format(i_filename, name),
                    )

            # Creating the list that will contain the marker list
            current_marker_i = 0
            all_markers = [None for i in range(nb_markers)]
            nb_samples = 0

            # Reading the first data line
            line = i_file.readline()
            row = line.rstrip("\r\n").split(",")

            while line != "":
                # Getting the marker name and sample id
                sample = row[header[other_opts.beeline_sample]]
                if sample in seen_samples:
                    logging.warning("{}: duplicate sample found, output file "
                                    "will be overwritten".format(sample))
                seen_samples.add(sample)

                logging.info("Processing {}".format(sample))
                o_filename = os.path.join(out_dir, sample + ".txt")
                with open(o_filename, "w") as o_file:
                    # Writing the metadata?
                    if other_opts.keep_meta:
                        o_file.write(metadata)

                    # Printing the header row
                    if other_opts.add_mapping:
                        print("Chr", "Pos", sep=other_opts.o_delim,
                              end=other_opts.o_delim, file=o_file)

                    if other_opts.add_ab:
                        print("Allele1 - AB", "Allele2 - AB",
                              sep=other_opts.o_delim, end=other_opts.o_delim,
                              file=o_file)

                    print(*header_row, sep=other_opts.o_delim, file=o_file)

                    # Reading the rest of the data for this sample
                    current_sample = sample
                    while current_sample == sample:
                        # Checking the marker order
                        marker = row[header[other_opts.beeline_id]]

                        # If the index is > than the length, it might be a
                        # duplicated sample...
                        if current_marker_i == len(all_markers):
                            break

                        if all_markers[current_marker_i] is None:
                            all_markers[current_marker_i] = marker
                        if all_markers[current_marker_i] != marker:
                            raise ProgramError(
                                "{}: marker order is not the same for "
                                "sample '{}'".format(i_filename, sample)
                            )

                        # The marker location information
                        chrom, pos, allele_encoding = locations[marker]

                        # Printing
                        if other_opts.add_mapping:
                            print(chrom, pos, sep=other_opts.o_delim,
                                  end=other_opts.o_delim, file=o_file)

                        if other_opts.add_ab:
                            # Adding A/B alleles
                            allele_1 = encode_allele(
                                allele=row[header[other_opts.beeline_a1]],
                                encoding=allele_encoding,
                            )
                            allele_2 = encode_allele(
                                allele=row[header[other_opts.beeline_a2]],
                                encoding=allele_encoding,
                            )
                            print(*sorted([allele_1, allele_2]),
                                  sep=other_opts.o_delim,
                                  end=other_opts.o_delim, file=o_file)

                        # Printing
                        print(*row, sep=other_opts.o_delim, file=o_file)

                        # Increasing the current marker
                        current_marker_i += 1

                        # Reading the next row
                        line = i_file.readline()
                        if line == "":
                            # End of file
                            break

                        # Splitting and current sample
                        row = line.rstrip("\r\n").split(",")
                        current_sample = row[header[other_opts.beeline_sample]]

                    # If there is only one marker, there is a problem
                    if nb_markers != 1 and current_marker_i == 1:
                        raise ProgramError(
                            "{}: data should be sorted by samples, not by "
                            "markers ('{}' had 1 marker, expecting "
                            "{:,d}".format(i_filename, sample, nb_markers)
                        )

                    # Are there any missing marker?
                    if current_marker_i != nb_markers:
                        nb_missing = abs(current_marker_i - nb_markers)
                        raise ProgramError(
                            "{}: missing {} marker{} for sample '{}'".format(
                                i_filename,
                                nb_missing,
                                "s" if nb_missing > 1 else "",
                                sample,
                            )
                        )
                    current_marker_i = 0
                    nb_samples += 1

            # Closing the output files
            logging.info("Done writing {:,d} samples".format(nb_samples))

        finally:
            i_file.close()


def extract_beeline(i_filenames, out_dir, o_suffix, locations, samples,
                    other_opts):
    """Extracts information from beeline report(s).

    Args:
        i_filenames (list): a list of file names (str)
        out_dir (str): the name of the output directory
        o_suffix (str): the suffix to add to the output file(s)
        locations (dict): a dictionary from marker ID to genomic location
        samples (set): a set of samples to keep
        other_opts(argparse.Namespace): the program options

    """
    # The chromosome to extract
    chrom = other_opts.chrom

    # The regex for changing delimiter (if required)
    regex = re.compile(",")

    # Reading all the files
    for i_filename in i_filenames:
        logging.info("Extracting from '{}'".format(i_filename))

        # The number of extracted markers
        nb_extracted_markers = 0

        # Getting the output filename
        o_filename = os.path.basename(os.path.splitext(i_filename)[0])
        o_filename += o_suffix + ".csv"

        # Opening the file
        i_file = None
        if i_filename != "-":
            i_file = open(i_filename, "r")

        else:
            i_file = sys.stdin
            o_filename = "from_stdin" + o_suffix + ".csv"

        # Joining the paths
        o_filename = os.path.join(out_dir, o_filename)

        # Reading the file
        try:
            # Getting to the data
            line = i_file.readline()
            while not line.startswith("[Data]"):
                line = i_file.readline()

            # Reading and checking the header
            header_line = i_file.readline()
            header = {
                name: i for i, name in
                enumerate(header_line.rstrip("\r\n").split(","))
            }
            required_columns = [other_opts.beeline_id]
            if other_opts.samples_to_keep is not None:
                required_columns.append(other_opts.beeline_sample)
            for name in required_columns:
                if name not in header:
                    raise ProgramError(
                        "{}: '{}': missing column".format(i_filename, name),
                    )
            # We will add chromosome (Chr) and position (Position) to the
            # header (to add in the file) if not present
            to_add = ""
            name_to_add = set()
            for name in ("Chr", "Position"):
                if name not in header:
                    to_add += name + ","
                    name_to_add.add(name)
            new_header_line = to_add + header_line

            # Reading the first data line
            line = i_file.readline()

            with open(o_filename, "w") as o_file:
                # Writing the header
                if other_opts.o_delim == ",":
                    o_file.write(new_header_line)
                else:
                    o_file.write(regex.sub(
                        other_opts.o_delim,
                        new_header_line,
                    ))

                while line != "":
                    # Reading the next line
                    row = line.rstrip("\r\n").split(",")

                    # Keep a subset of samples?
                    if other_opts.samples_to_keep is not None:
                        sample = row[header[other_opts.beeline_sample]]
                        if sample not in samples:
                            line = i_file.readline()
                            continue

                    # Getting the marker information
                    marker = row[header[other_opts.beeline_id]]
                    marker_location = locations.get(marker, _unknown_location)
                    if marker_location.chrom in chrom:
                        to_add = ""
                        if "Chr" in name_to_add:
                            to_add += "{},".format(marker_location.chrom)
                        if "Position" in name_to_add:
                            to_add += "{},".format(marker_location.pos)
                        # We need this marker, we write the line
                        if other_opts.o_delim == ",":
                            o_file.write(to_add + line)
                        else:
                            o_file.write(regex.sub(
                                other_opts.o_delim,
                                to_add + line,
                            ))

                        # Updating the number of extracted markers
                        nb_extracted_markers += 1

                    # Reading the next line
                    line = i_file.readline()

        finally:
            i_file.close()

        # Logging
        logging.info("  - {:,d} markers in '{}'".format(
            nb_extracted_markers,
            o_filename,
        ))


def read_mapping_info(i_filename, delim, map_id, map_chr, map_pos, map_allele,
                      allele_strand, strand_col):
    """Reads the mapping information to gather genomic locations.

    Args:
        i_filename (str): the name of the input file
        delim (str): the field delimiter
        map_id (str): the name of the column containing marker ID
        map_chr (str): the name of the column containing the chromosome
        map_pos (str): the name of the column containing the position
        map_allele (str): the name of the column containing the alleles
        allele_strand (str): the strand of the alleles in the beeline file
        strand_col (str): the name of the column containing the stand info

    Returns:
        dict: a dictionary from marker ID to genomic location

    """
    logging.info("Reading mapping information")
    logging.info("  - marker strand is '{}'".format(allele_strand))
    logging.info("  - strand information is '{}'".format(strand_col))

    map_info = {}
    with open(i_filename, "r") as i_file:
        # Reading until the header
        header = get_header(i_file, delim, "[Assay]")

        # Reading the data
        for line in i_file:
            if line.startswith("[Controls]"):
                break

            row = line.rstrip("\r\n").split(delim)

            # Gathering the mapping information
            name = row[header[map_id]]
            chrom = encode_chromosome(row[header[map_chr]])
            pos = int(row[header[map_pos]])

            # Checking if we need to complement the alleles or not
            # Getting the strand information
            complement_required = need_complement(
                row[header[strand_col]], allele_strand,
            )

            # Splitting the alleles
            alleles = row[header[map_allele]].split("/")
            a_allele = alleles[0][1:]
            b_allele = alleles[1][:-1]

            # Complementing if required
            if complement_required:
                a_allele = _complement[a_allele]
                b_allele = _complement[b_allele]

            # Saving the information
            map_info[name] = _Location(
                chrom=chrom,
                pos=pos,
                alleles={a_allele: "A", b_allele: "B"},
            )

    logging.info("  - {:,d} markers".format(len(map_info)))
    return map_info


def need_complement(strand, strand_type):
    """Checks if complementing the alleles is required.

    Args:
        strand (str): the strand value
        strand_type (str): the strand type (forward, top, plus)

    Returns:
        bool: needs complement or not

    """
    if strand_type == "forward":
        value = _forward_strand_re.search(strand)
        if value is None:
            raise ValueError("{}: invalid forward value".format(strand))
        return value.group(1) == "R"

    elif strand_type == "plus":
        if strand != "+" and strand != "-":
            raise ValueError("{}: invalid plus value".format(strand))
        return strand == "-"

    elif strand_type == "top":
        if strand != "TOP" and strand != "BOT":
            raise ValueError("{}: invalid top value".format(strand))
        return strand == "BOT"

    else:
        raise ValueError("{}: invalid strand".format(strand_type))


def read_list(i_filename):
    """Reads a file containing a list of element.

    Args:
        i_filename (str): the name of the file to read (might be ``None``)

    Returns:
        set: a set containing the list of element in the file

    If ``i_filename`` is ``None``, an empty set is returned.

    """
    elements = set()
    if i_filename is not None:
        logging.info("Reading list of samples to keep")
        with open(i_filename, "r") as i_file:
            elements = set(i_file.read().splitlines())
        logging.info("  - {:,d} samples to keep".format(len(elements)))

    return elements


def get_header(f, delim=",", data_delim="[Data]"):
    """Reads a file until the header.

    Args:
        f (file): a file object
        delim (str): the file delimiter
        data_delim (str): the delimiter before data in the file

    Returns:
        dict: the header as a dictionary (column => index)

    Reads until the line equals to ``[Assay]``, meaning that the next line is
    the header line (valid for Illumina manifest files).

    """
    nb_line = 1
    line = f.readline()
    while not line.startswith(data_delim):
        if line == "":
            raise ProgramError("{}: no data in file".format(f.name))

        line = f.readline()
        nb_line += 1

        if nb_line > 1000:
            raise ProgramError("{}: no data after 1000 lines".format(f.name))

    return {
        name: i for i, name in
        enumerate(f.readline().rstrip("\r\n").split(delim))
    }


def encode_chromosome(chrom):
    """Encodes a chromosome.

    Args:
        chrom (str): the chromosome to encode

    Returns:
        int: the encoded chromosome

    """
    chrom = chrom.upper()
    if chrom == "X":
        return 23
    if chrom == "Y":
        return 24
    if chrom == "XY" or chrom == "YX":
        return 25
    if chrom == "M" or chrom == "MT":
        return 26

    try:
        chrom = int(chrom)
    except ValueError:
        return 0

    if chrom < 1 or chrom > 26:
        return 0

    return chrom


def encode_allele(allele, encoding):
    """Encodes an allele according to encoding.

    Args:
        allele (str): the allele to encode
        encoding (dict): the allele encoding

    Returns:
        str: the encoded allele.

    Note
    ----
        The function search for the allele and its complement. If none of this
        is found in the encoding, then "-" is returned.

    """
    if allele == "-":
        return allele

    else:
        return encoding[allele]


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the options to verify

    Note
    ----
        If there is a problem, a :py:class:`ProgramError` is raised.

    """
    # Checking the output format if present
    if args.analysis_type == "convert":
        if (args.o_format == "bed") and (not _HAS_PYPLINK):
            raise ProgramError("BED format requires pyplink module")

    # Checking the input file(s)
    if (len(args.i_filenames) == 1) and (args.i_filenames[0] == "-"):
        # Checking if we have only one input file and its '-' (stdin)
        logging.info("Beeline report will be read from STDIN")

    else:
        for filename in args.i_filenames:
            if not os.path.isfile(filename):
                raise ProgramError("{}: no such file".format(filename))

            # Checking if the file is gzip
            if filename.endswith(".gz"):
                raise ProgramError(
                    "{}: GZIP not yet implemented (use '-')".format(filename)
                )

    # Checking the strand of the alleles
    if args.beeline_strand is None:
        a1 = args.beeline_a1.lower()
        a2 = args.beeline_a2.lower()
        if "forward" in a1 and "forward" in a2:
            args.beeline_strand = "forward"
        elif "top" in a1 and "top" in a2:
            args.beeline_strand = "top"
        elif "plus" in a1 and "plus" in a2:
            args.beeline_strand = "plus"
        else:
            raise ProgramError("Impossible to infer the allele strand from "
                               "the column names")
        logging.info("Inferred that the allele strand is '{}'"
                     "".format(args.beeline_strand))

    # Checking the map file
    if not os.path.isfile(args.map_filename):
        raise ProgramError("{}: no such file".format(args.map_filename))

    # Checking the columns are inside the map file
    with open(args.map_filename, "r") as i_file:
        # Reading the header
        header = get_header(i_file, args.map_delim, "[Assay]")

        # The required columns
        required_columns = [args.map_id, args.map_chr, args.map_pos,
                            args.map_allele]
        if args.beeline_strand == "forward":
            required_columns.append(args.map_illumina_id)
        elif args.beeline_strand == "top":
            required_columns.append(args.map_illumina_strand)
        elif args.beeline_strand == "plus":
            required_columns.append(args.map_ref_strand)

        # Checking the column
        for name in required_columns:
            if name not in header:
                raise ProgramError("{}: missing column '{}'".format(
                    args.map_filename, name,
                ))

    # Checking we can write to the output directory (if required)
    if args.output_dir is not None:
        if not os.path.isdir(args.output_dir):
            raise ProgramError(
                "{}: no such directory".format(args.output_dir)
            )
        try:
            with NamedTemporaryFile(dir=args.output_dir):
                pass
        except OSError:
            raise ProgramError(
                "{}: not writable".format(args.output_dir)
            )
    else:
        args.output_dir = "."

    # These options are specific to the extract analysis
    if args.analysis_type == "extract":
        # Checking the chromosomes
        e_chrom = [encode_chromosome(chrom) for chrom in args.chrom]
        for chrom, o_chrom in zip(e_chrom, args.chrom):
            if chrom == 0:
                raise ProgramError("{}: invalid chromosome".format(o_chrom))
        args.chrom = set(e_chrom)

        # List of samples to keep
        if args.samples_to_keep is not None:
            if not os.path.isfile(args.samples_to_keep):
                raise ProgramError("{}: no such "
                                   "file".format(args.samples_to_keep))


def parse_args(parser):
    """Parses the command line options and arguments.

    Args:
        parser (argparse.ArgumentParser): the argument parser object

    Returns:
        argparse.Namespace: the list of options and arguments

    Note
    ----
        The only check that is done here is by the parser itself. Values are
        verified later by the :py:func:`check_args` function.

    """
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s version {}".format(__version__),
    )

    # Creating a parent parser (for common options between analysis type)
    p_parser = argparse.ArgumentParser(add_help=False)

    p_parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s version {}".format(__version__),
    )

    # The input files
    group = p_parser.add_argument_group("Input Files")
    group.add_argument(
        "-i", "--input", type=str, metavar="FILE", dest="i_filenames",
        required=True, nargs="+",
        help="The name of the input file(s). Use '-' only once to read on the "
             "standard input (STDIN).",
    )
    group.add_argument(
        "-m", "--map", type=str, metavar="FILE", dest="map_filename",
        required=True,
        help="The name of the file containing mapping information.",
    )

    # The Beeline options
    group = p_parser.add_argument_group("Beeline Options")
    group.add_argument(
        "--beeline-id", type=str, metavar="COL", default="SNP Name",
        help="The name of the column containing the marker identification "
             "number for beeline. [%(default)s]",
    )
    group.add_argument(
        "--beeline-sample", type=str, metavar="COL", default="Sample Name",
        help="The name of the column containing the sample identification "
             "number for beeline. [%(default)s]",
    )
    group.add_argument(
        "--beeline-a1", type=str, metavar="COL", default="Allele1 - Forward",
        help="The name of the column containing the first allele for beeline. "
             "[%(default)s]",
    )
    group.add_argument(
        "--beeline-a2", type=str, metavar="COL", default="Allele2 - Forward",
        help="The name of the column containing the second allele for "
             "beeline. [%(default)s]",
    )
    group.add_argument(
        "--beeline-strand", type=str, metavar="STR",
        choices={"forward", "plus", "top"},
        help="The strand for the alleles. This will help in determining the "
             "A/B alleles. If unset, the strand will be infer by the name of "
             "the columns for the two alleles.",
    )

    # The mapping options
    group = p_parser.add_argument_group("Mapping Options")
    group.add_argument(
        "--map-id", type=str, metavar="COL", default="Name",
        help="The name of the column containing the marker identification "
             "numbers. [%(default)s]",
    )
    group.add_argument(
        "--map-chr", type=str, metavar="COL", default="Chr",
        help="The name of the column containing the chromosome. [%(default)s]",
    )
    group.add_argument(
        "--map-pos", type=str, metavar="COL", default="MapInfo",
        help="The name of the column containing the position. [%(default)s]",
    )
    group.add_argument(
        "--map-allele", type=str, metavar="COL", default="SNP",
        help="The name of the column containing the alleles. [%(default)s]",
    )
    group.add_argument(
        "--map-illumina-id", type=str, metavar="COL", default="IlmnID",
        help="The name of the column containing the Illumina ID. This ID "
             "contains information about the marker's Forward strand and help "
             "determining the A/B alleles. [%(default)s]",
    )
    group.add_argument(
        "--map-illumina-strand", type=str, metavar="COL", default="IlmnStrand",
        help="The name of the column containing the Illumina strand. This "
             "information helps in determining the A/B alleles since it "
             "contains the marker's orientation for the TOP alleles. "
             "[%(default)s]",
    )
    group.add_argument(
        "--map-ref-strand", type=str, metavar="COL", default="RefStrand",
        help="The name of the column containing the reference strand. This "
             "information helps in determining the A/B alleles since it "
             "contains the marker's orientation for the Plus alleles. "
             "[%(default)s]",
    )
    group.add_argument(
        "--map-delim", type=str, metavar="SEP", default=",",
        help="The field delimiter. [%(default)s]",
    )
    group.add_argument(
        "--nb-snps-kw", type=str, metavar="KEYWORD", default="Num Used SNPs",
        help="The keyword that describe the number of used markers for the "
             "report(s) (useful if beeline header format changes). "
             "[%(default)s]",
    )

    # The output options
    group = p_parser.add_argument_group("Output Directory")
    group.add_argument(
        "-o", "--output-dir", type=str, metavar="DIR", dest="output_dir",
        help="The output directory (default is working directory).",
    )

    # Adding sub parsers
    subparsers = parser.add_subparsers(
        title="Analysis Type",
        description="The type of analysis to be performed on the Beeline "
                    "(long) report. The following list of analysis is "
                    "available.",
        dest="analysis_type",
    )
    subparsers.required = True

    # The convert parser
    convert_parser = subparsers.add_parser(
        "convert",
        help="Converts the Beeline long report to other format.",
        description="Conversion from the Beeline long report format to other, "
                    "more commonly used, format (such as Plink) "
                    "(version {}).".format(__version__),
        parents=[p_parser],
    )

    # The different format
    group = convert_parser.add_argument_group("Output Format")
    group.add_argument(
        "--format", type=str, metavar="FORMAT", choices={"bed", "ped"},
        default="bed", dest="o_format",
        help="The output format (one of 'bed' or 'ped' for binary or "
             "normal pedfile format from Plink). [%(default)s]",
    )

    # The sample split parser
    sample_split_parser = subparsers.add_parser(
        "sample-split",
        help="Split the report so that there is one file for each sample.",
        description="The long report will be split so that each sample has "
                    "its own separate file (version {}).".format(__version__),
        parents=[p_parser],
    )

    # The split options
    group = sample_split_parser.add_argument_group("Split Options")
    group.add_argument(
        "--keep-metadata", action="store_true", dest="keep_meta",
        help="Keeps the meta data ([Header]) in each of the output reports.",
    )
    group.add_argument(
        "--add-ab", action="store_true", dest="add_ab",
        help="Adds the A/B alleles in the output file (sometimes required by "
             "other software).",
    )
    group.add_argument(
        "--add-mapping", action="store_true", dest="add_mapping",
        help="Adds mapping information (chromosome/Position).",
    )
    group.add_argument(
        "--output-delim", type=str, metavar="SEP", dest="o_delim", default=",",
        help="The output file field delimiter. [%(default)s]",
    )

    # The extract parser
    extract_parser = subparsers.add_parser(
        "extract",
        help="Extract information from the Beeline long report.",
        description="Extracts information from a Beeline long report "
                    "(version {}).".format(__version__),
        parents=[p_parser],
    )

    # The extraction options
    group = extract_parser.add_argument_group("Extraction Options")
    group.add_argument(
        "-c", "--chr", type=str, nargs="+", dest="chrom",
        default=[str(encode_chromosome(str(chrom))) for chrom in range(1, 27)],
        help="The chromosome to extract. [%(default)s]",
    )
    group.add_argument(
        "-k", "--keep", type=str, metavar="FILE", dest="samples_to_keep",
        help="A list of samples to extract.",
    )

    # The output options
    group = extract_parser.add_argument_group("Output Options")
    group.add_argument(
        "-s", "--suffix", type=str, metavar="STR", dest="o_suffix",
        default="_extract",
        help="The suffix to add to the output file(s). [%(default)s]",
    )
    group.add_argument(
        "--output-delim", type=str, metavar="SEP", dest="o_delim", default=",",
        help="The output file field delimiter. [%(default)s]",
    )

    # Parsing and returning the arguments and options
    return parser.parse_args()


class ProgramError(Exception):
    """An Exception raised in case of a problem."""
    def __init__(self, msg):
        """Construction of the ProgramError class."""
        self.message = str(msg)

    def __str__(self):
        return self.message


if __name__ == "__main__":
    main()
