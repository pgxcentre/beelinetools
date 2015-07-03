#!/usr/bin/env python
"""Script that converts the beeline report to one that is readable by Plink."""


from __future__ import print_function

import os
import sys
import logging
import argparse
from collections import namedtuple
from tempfile import NamedTemporaryFile


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2015, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "MIT"
__version__ = "0.1.1"


# The location tuple
_Location = namedtuple("Location", ["chrom", "pos"])
_unknown_location = _Location(chrom=0, pos=0)


def main():
    """The main function."""
    # Setting the logging
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    # The parser
    desc = ("Convert beeline report(s) into Plink readable files "
            "(version {})".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    try:
        # Parsing and checking the options and arguments
        args = parse_args(parser)
        check_args(args)

        # Reading the map file
        map_data = read_mapping_info(
            args.map_filename,
            delim=args.delim,
            id_col=args.id_col,
            chr_col=args.chr_col,
            pos_col=args.pos_col,
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
            # Extracting from the beeline report(s)
            extract_beeline(
                i_filenames=args.i_filenames,
                o_filename=args.o_filename,
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
    for i_filename in i_filenames:
        logging.info("Converting '{}'".format(i_filename))

        # Getting the output filename
        o_filename = os.path.splitext(os.path.basename(i_filename))[0]

        with open(i_filename, "r") as i_file:
            # The number of markers and samples
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
            required_columns = ("SNP Name", "Sample ID", "Allele1 - Forward",
                                "Allele2 - Forward")
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

            # The pedfile
            pedfile = o_filename + ".ped"
            if out_dir is not None:
                pedfile = os.path.join(out_dir, pedfile)
            with open(pedfile, "w") as pedfile:
                while line != "":
                    # Getting the marker name and sample id
                    sample = row[header["Sample ID"]]

                    # Logging
                    logging.info("Processing {}".format(sample))
                    # Printing the starting of the file
                    print(sample, sample, "0", "0", "0", "-9", sep="\t",
                          end="", file=pedfile)

                    # Reading the rest of the data for this sample
                    current_sample = sample
                    while current_sample == sample:
                        # Checking the marker order
                        marker = row[header["SNP Name"]]
                        if all_markers[current_marker_i] is None:
                            all_markers[current_marker_i] = marker
                        if all_markers[current_marker_i] != marker:
                            raise ProgramError(
                                "{}: marker order is not the same for "
                                "sample '{}'".format(i_filename, sample)
                            )

                        # Getting the genotype
                        genotype = "{} {}".format(
                            row[header["Allele1 - Forward"]],
                            row[header["Allele2 - Forward"]],
                        )
                        if "-" in genotype:
                            genotype = "0 0"

                        print("\t" + genotype, sep="\t", end="", file=pedfile)

                        # Increasing the current marker
                        current_marker_i += 1

                        # Reading the next row
                        line = i_file.readline()
                        if line == "":
                            # End of file
                            break

                        # Splitting and current sample
                        row = line.rstrip("\r\n").split(",")
                        current_sample = row[header["Sample ID"]]

                    print("\n", sep="\t", end="", file=pedfile)

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
        map_filename = o_filename + ".map"
        if out_dir is not None:
            map_filename = os.path.join(out_dir, map_filename)
        with open(map_filename, "w") as o_file:
            for marker in all_markers:
                marker_location = None
                if marker in locations:
                    marker_location = locations[marker]
                else:
                    marker_location = _unknown_location
                    logging.warning("{}: no mapping "
                                    "information".format(marker))
                print(marker_location.chrom, marker, "0",
                      marker_location.pos, sep="\t", file=o_file)


def extract_beeline(i_filenames, o_filename, locations, other_opts):
    """Extracts information from beeline report(s).

    Args:
        i_filenames (list): a list of file names (str)
        o_filename (str): the name of the output file
        locations (dict): a dictionary from marker ID to genomic location
        other_opts(argparse.Namespace): the program options

    """
    # The number of markers for all files (and the list of markers)
    all_nb_markers = None
    extracted_markers = set()

    # The chromosome to extract
    chrom = other_opts.chrom

    # Reading all the files
    for i_filename in i_filenames:
        logging.info("Converting '{}'".format(i_filename))

        with open(i_filename, "r") as i_file:
            # The number of markers and samples
            nb_markers = None

            # Reading the assay information (and saving it)
            line = i_file.readline()
            while not line.startswith("[Header]"):
                line = i_file.readline()

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

            # Checking the number of markers
            if all_nb_markers is None:
                all_nb_markers = nb_markers
            if nb_markers != all_nb_markers:
                raise ProgramError(
                    "{}: not same number of markers as other "
                    "report(s)".format(i_filename)
                )
            logging.info("There are {:,d} markers".format(nb_markers))

            # Reading and checking the header
            header_line = i_file.readline()
            header = {
                name: i for i, name in
                enumerate(header_line.rstrip("\r\n").split(","))
            }
            required_columns = ("SNP Name", )
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

            # The number of samples
            nb_samples = 0

            # Reading the first data line
            line = i_file.readline()
            row = line.rstrip("\r\n").split(",")

            with open(o_filename, "w") as o_file:
                # Writing the header
                o_file.write(new_header_line)

                while line != "":
                    # Reading the next line
                    line = i_file.readline()
                    row = line.rstrip("\r\n").split(",")

                    # Getting the marker information
                    marker = row[header["SNP Name"]]
                    marker_location = locations.get(marker, _unknown_location)
                    if marker_location.chrom in chrom:
                        to_add = ""
                        if "Chr" in name_to_add:
                            to_add += "{},".format(marker_location.chrom)
                        if "Position" in name_to_add:
                            to_add += "{},".format(marker_location.pos)
                        # We need this marker, we write the line
                        o_file.write(to_add + line)

                        # Updating the extracted markers
                        extracted_markers.add(marker)


def read_mapping_info(i_filename, delim, id_col, chr_col, pos_col):
    """Reads the mapping information to gather genomic locations.

    Args:
        i_filename (str): the name of the input file
        delim (str): the field delimiter
        id_col (str): the name of the column containing marker ID
        chr_col (str): the name of the column containing the chromosome
        pos_col (str): the name of the column containing the position

    Returns:
        dict: a dictionary from marker ID to genomic location

    """
    logging.info("Reading mapping information")
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
            name = row[header[id_col]]
            chrom = encode_chromosome(row[header[chr_col]])
            pos = int(row[header[pos_col]])

            # Saving the information
            map_info[name] = _Location(chrom=chrom, pos=pos)

    logging.info("  - {:,d} markers".format(len(map_info)))
    return map_info


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


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the options to verify

    Note
    ----
        If there is a problem, a :py:class:`ProgramError` is raised.

    """
    # Checking the input file(s)
    for filename in args.i_filenames:
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))

    # Checking the map file
    if not os.path.isfile(args.map_filename):
        raise ProgramError("{}: no such file".format(args.map_filename))

    # Checking the columns are inside the map file
    with open(args.map_filename, "r") as i_file:
        # Reading the header
        header = get_header(i_file, args.delim, "[Assay]")

        # Checking the column
        for name in (args.id_col, args.chr_col, args.pos_col):
            if name not in header:
                raise ProgramError("{}: missing column '{}'".format(
                    args.map_filename,
                    name,
                ))

    # These options are specific to the convert analysis
    if args.analysis_type == "convert":
        # Checking we can write to the output directory (if required)
        if args.output_dir is not None:
            if not os.path.isdir(args.output_dir):
                raise ProgramError(
                    "{}: no such directory".format(args.output_dir)
                )
            try:
                with NamedTemporaryFile(dir=args.output_dir) as tmp_file:
                    pass
            except OSError:
                raise ProgramError(
                    "{}: not writable".format(args.output_dir)
                )

    # These options are specific to the extract analysis
    if args.analysis_type == "extract":
        # Checking the chromosomes
        e_chrom = [encode_chromosome(chrom) for chrom in args.chrom]
        for chrom, o_chrom in zip(e_chrom, args.chrom):
            if chrom == 0:
                raise ProgramError("{}: invalid chromosome".format(o_chrom))
        args.chrom = set(e_chrom)

        # Checking that we can write the file
        try:
            with open(args.o_filename, "w") as o_file:
                pass
        except OSError:
            raise ProgramError(
                "{}: not writable".format(args.o_filename)
            )


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
        "-v",
        "--version",
        action="version",
        version="%(prog)s version {}".format(__version__),
    )

    # Creating a parent parser (for common options between analysis type)
    p_parser = argparse.ArgumentParser(add_help=False)

    p_parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s version {}".format(__version__),
    )

    # The input files
    group = p_parser.add_argument_group("Input Files")
    group.add_argument(
        "-i",
        "--input",
        type=str,
        metavar="FILE",
        dest="i_filenames",
        required=True,
        nargs="+",
        help="The name of the input file(s).",
    )
    group.add_argument(
        "-m",
        "--map",
        type=str,
        metavar="FILE",
        dest="map_filename",
        required=True,
        help="The name of the file containing mapping information.",
    )

    # The mapping options
    group = p_parser.add_argument_group("Mapping Options")
    group.add_argument(
        "--id-col",
        type=str,
        metavar="COL",
        default="Name",
        help="The name of the column containing the marker identification "
             "numbers [%(default)s]",
    )
    group.add_argument(
        "--chr-col",
        type=str,
        metavar="COL",
        default="Chr",
        help="The name of the column containing the chromosome [%(default)s]",
    )
    group.add_argument(
        "--pos-col",
        type=str,
        metavar="COL",
        default="MapInfo",
        help="The name of the column containing the position [%(default)s]",
    )
    group.add_argument(
        "--delim",
        type=str,
        metavar="SEP",
        default=",",
        help="The field delimiter [%(default)s]",
    )
    group.add_argument(
        "--nb-snps-kw",
        type=str,
        metavar="KEYWORD",
        default="Num Used SNPs",
        help="The keyword that describe the number of used markers for the "
             "report(s) (useful if beeline header format changes) "
             "[%(default)s]",
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
                    "more commonly used, format (such as Plink).",
        parents=[p_parser],
    )

    # The output options
    group = convert_parser.add_argument_group("Output Options")
    group.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="DIR",
        dest="output_dir",
        help="The output directory (default is working directory)",
    )

    # The extract parser
    extract_parser = subparsers.add_parser(
        "extract",
        help="Extract information from the Beeline long report.",
        description="Extracts information from a Beeline long report.",
        parents=[p_parser],
    )

    # The extraction options
    group = extract_parser.add_argument_group("Extraction Options")
    group.add_argument(
        "-c",
        "--chr",
        nargs="+",
        dest="chrom",
        default=[str(encode_chromosome(str(chrom))) for chrom in range(1, 27)],
        help="The chromosome to extract %(default)s",
    )

    # The output options
    group = extract_parser.add_argument_group("Output Options")
    group.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="CSV",
        dest="o_filename",
        default="beeline_extract.csv",
        help="The name of the output file [%(default)s]",
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
