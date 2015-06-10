#!/usr/bin/env python
"""Script that converts the beeline report to one that is readable by Plink."""


from __future__ import print_function

import os
import sys
import logging
import argparse
import tempfile
from collections import namedtuple


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2015, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "MIT"
__version__ = "0.1.0"


# The location tuple
_Location = namedtuple("Location", ["chrom", "pos"])


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

        # Converting the beeline report(s)
        convert_beeline(args.i_filenames, args.output_dir, locations=map_data)

    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        sys.exit(1)

    except ProgramError as e:
        logging.error(e)
        parser.error(e.message)

    except Exception as e:
        logging.error(e)
        raise


def convert_beeline(i_filenames, out_dir, locations):
    """Convert beeline report(s) to Plink files.

    Args:
        i_filenames (list): a list of file names (str)

    """
    for i_filename in i_filenames:
        logging.info("Converting '{}'".format(i_filename))

        # Getting the output filename
        o_filename = os.path.splitext(os.path.basename(i_filename))[0]

        with open(i_filename, "r") as i_file:
            # The number of markers and samples
            nb_markers = None

            # Reading the assay information
            line = i_file.readline().rstrip("\n")
            while line != "[Header]":
                line = i_file.readline().rstrip("\n")

            while not line.startswith("[Data]"):
                if line.startswith("Num Used SNPs"):
                    nb_markers = int(line.rstrip("\n").split(",")[-1])
                line = i_file.readline()

            if nb_markers is None:
                raise ProgramError("{}: invalid header (missing 'Num Used "
                                   "SNPs' value)".format(i_filename))

            logging.info("There are {:,d} markers".format(nb_markers))

            # Reading and checking the header
            header = {
                name: i
                for i, name in
                enumerate(i_file.readline().rstrip("\n").split(","))
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
            row = line.rstrip("\n").split(",")

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
                        row = line.rstrip("\n").split(",")
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
                if marker not in locations:
                    raise ProgramError("{}: no mapping "
                                       "information".format(marker))
                print(locations[marker].chrom, marker, "0",
                      locations[marker].pos, sep="\t", file=o_file)


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

            row = line.rstrip("\n").split(delim)

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
        enumerate(f.readline().rstrip("\n").split(delim))
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

    # Checking we can write to the output directory (if required)
    if args.output_dir is not None:
        if not os.path.isdir(args.output_dir):
            raise ProgramError("{}: no such directory".format(args.output_dir))
        try:
            with tempfile.NamedTemporaryFile(dir=args.output_dir) as tmp_file:
                pass
        except PermissionError:
            raise ProgramError("{}: not writable".format(args.output_dir))


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
        version="$(prog)s version {}".format(__version__),
    )

    # The input files
    group = parser.add_argument_group("Input Files")
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
    group = parser.add_argument_group("Mapping Options")
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

    # The output options
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="DIR",
        dest="output_dir",
        help="The output directory (default is working directory)",
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
