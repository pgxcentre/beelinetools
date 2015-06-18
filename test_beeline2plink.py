"""Tests the beeline2plink script."""


from __future__ import print_function

import os
import re
import shutil
import random
import logging
import platform
import unittest
import collections
from io import StringIO
from collections import defaultdict
from tempfile import mkdtemp, NamedTemporaryFile

import beeline2plink


class _DummyArgs(object):
    pass


class TestBeeline2Plink(unittest.TestCase):
    def setUp(self):
        """Setup the tests."""
        self.tmp_dir = mkdtemp(prefix="beeline2plink_test_")

    def tearDown(self):
        """Finishes the tests."""
        # Cleaning the temporary directory
        shutil.rmtree(self.tmp_dir)

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python < 3.4."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            return AssertLogsContext_Compatibility(self, logger, level)

    def test_encode_chromosome(self):
        """Tests the 'encode_chromosome' function."""
        # Testing all valid chromosome
        for chrom in range(1, 27):
            observed = beeline2plink.encode_chromosome(str(chrom))
            self.assertEqual(chrom, observed)

        # Testing X chromosome
        for chrom in ("x", "X", "23"):
            observed = beeline2plink.encode_chromosome(chrom)
            self.assertEqual(23, observed)

        # Testing the Y chromosome
        for chrom in ("y", "Y", "24"):
            observed = beeline2plink.encode_chromosome(chrom)
            self.assertEqual(24, observed)

        # Testing the pseudo autosomal region
        for chrom in ("XY", "YX", "25"):
            observed = beeline2plink.encode_chromosome(chrom)
            self.assertEqual(25, observed)

        # Testing the mitochondrial chromosome
        for chrom in ("M", "MT", "26"):
            observed = beeline2plink.encode_chromosome(chrom)
            self.assertEqual(26, observed)

    def test_encode_chromosome_invalid(self):
        """Tests the 'encode_chromosome' function for invalid chromosome."""
        # Testing invalid chromosome
        for chrom in ("-9", "0", "27"):
            observed = beeline2plink.encode_chromosome(chrom)
            self.assertEqual(0, observed)

    def test_get_header(self):
        """Tests the 'get_header' function."""
        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False) as f:
            tmp_filename = f.name
            print(*range(6), sep="\n", file=f)
            print("[Data],,,", file=f)
            print("header_1,header_2,header_3", file=f)
            for j in range(6):
                print(*["data_{}_{}".format(i+1, j) for i in range(3)],
                      sep=",", file=f)

        # Reading until the header (with good separator)
        expected = {"header_{}".format(i+1): i for i in range(3)}
        with open(tmp_filename, "r") as i_file:
            observed = beeline2plink.get_header(i_file)
        self.assertEqual(expected, observed)

        # Reading until the header (with bas separator)
        expected = {"header_1,header_2,header_3": 0}
        with open(tmp_filename, "r") as i_file:
            observed = beeline2plink.get_header(i_file, delim="\t")
        self.assertEqual(expected, observed)

    def test_get_header_error(self):
        """Tests the 'get_header' function with errors."""
        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False) as f:
            tmp_filename = f.name
            print(*range(6), sep="\n", file=f)
            print("[Data],,,", file=f)
            print("header_1,header_2,header_3", file=f)
            for j in range(6):
                print(*["data_{}_{}".format(i+1, j) for i in range(3)],
                      sep=",", file=f)

        # This should raise an exception
        with open(tmp_filename, "r") as i_file:
            with self.assertRaises(beeline2plink.ProgramError) as e:
                beeline2plink.get_header(i_file, data_delim="[Assay]")
        self.assertEqual("{}: no data in file".format(tmp_filename),
                         e.exception.message)

        # Writing a long file without data
        with open(tmp_filename, "w") as o_file:
            for i in range(1001):
                print(i, file=o_file)

        # This should raise an exception
        with open(tmp_filename, "r") as i_file:
            with self.assertRaises(beeline2plink.ProgramError) as e:
                beeline2plink.get_header(i_file)
        self.assertEqual(
            "{}: no data after 1000 lines".format(tmp_filename),
            e.exception.message,
        )

    def test_read_mapping_info(self):
        """Tests the 'read_mapping_info' function."""
        tmp_filename = None
        expected = {}
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False) as f:
            tmp_filename = f.name
            print(*range(6), sep="\n", file=f)
            print("[Assay]", file=f)
            print("Name,Chr,MapInfo", file=f)
            for i in range(100):
                marker_name = "marker_{}".format(i + 1)
                chrom = random.randint(1, 26)
                pos = random.randint(1, 3000000)
                print(marker_name, chrom, pos, sep=",", file=f)
                expected[marker_name] = beeline2plink._Location(
                    chrom=chrom,
                    pos=pos,
                )

        # Getting the expected data
        observed = beeline2plink.read_mapping_info(
            i_filename=tmp_filename,
            delim=",",
            id_col="Name",
            chr_col="Chr",
            pos_col="MapInfo",
        )
        self.assertEqual(expected, observed)

        # Adding a '[Controls]' line, everything below should be excluded
        with open(tmp_filename, "a") as o_file:
            print("added_marker,1,1", file=o_file)
            print("[Controls],,,", file=o_file)
            print("skipped_line", file=o_file)
        expected["added_marker"] = beeline2plink._Location(chrom=1, pos=1)

        # Getting the expected data
        observed = beeline2plink.read_mapping_info(
            i_filename=tmp_filename,
            delim=",",
            id_col="Name",
            chr_col="Chr",
            pos_col="MapInfo",
        )
        self.assertEqual(expected, observed)

    def test_convert_beeline(self):
        """Tests the 'convert_beeline' function."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        sample_genotype = defaultdict(dict)
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)
                    sample_genotype[sample_id][marker_id] = genotype

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Creating a temporary file
        tmp_filename_2 = None
        sample_genotype_2 = defaultdict(dict)
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename_2 = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + nb_samples + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)
                    sample_genotype_2[sample_id][marker_id] = genotype

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        beeline2plink.convert_beeline(
            i_filenames=[tmp_filename, tmp_filename_2],
            out_dir=self.tmp_dir,
            locations=mapping_info,
            other_opts=other_options,
        )

        # Checking the two map files
        for map_filename in (tmp_filename, tmp_filename_2):
            map_filename = os.path.splitext(map_filename)[0] + ".map"
            self.assertTrue(os.path.isfile(map_filename))

            # Checking the first map file content
            seen_markers = set()
            with open(map_filename, "r") as i_file:
                for i, line in enumerate(i_file):
                    # Gathering the information
                    marker = "marker_{}".format(i + 1)
                    row = line.rstrip("\n").split("\t")

                    # Comparing the content of the file
                    self.assertTrue(marker in mapping_info)
                    self.assertEqual(4, len(row))
                    self.assertEqual(mapping_info[marker].chrom, int(row[0]))
                    self.assertEqual(marker, row[1])
                    self.assertEqual("0", row[2])
                    self.assertEqual(mapping_info[marker].pos, int(row[3]))

                    # We have compared this marker
                    seen_markers.add(marker)
            self.assertEqual(set(mapping_info.keys()), seen_markers)

        # Checking the two ped files
        zipped = zip(
            (1, nb_samples + 1),
            (tmp_filename, tmp_filename_2),
            (sample_genotype, sample_genotype_2),
        )
        for sample_id_add, ped_filename, sample_geno in zipped:
            ped_filename = os.path.splitext(ped_filename)[0] + ".ped"
            self.assertTrue(os.path.isfile(ped_filename))

            # Checking the ped file content
            seen_samples = set()
            with open(ped_filename, "r") as i_file:
                for i, line in enumerate(i_file):
                    # Gathering the information
                    sample_id = "sample_{}".format(i + sample_id_add)
                    row = line.rstrip("\n").split("\t")

                    # Checking the sample information
                    self.assertEqual(6 + nb_markers, len(row))
                    sample_info = row[:6]
                    self.assertEqual(
                        [sample_id, sample_id, "0", "0", "0", "-9"],
                        sample_info,
                    )

                    # Checking the genotypes
                    seen_markers = set()
                    genotypes = row[6:]
                    self.assertEqual(nb_markers, len(genotypes))
                    for j, genotype in enumerate(genotypes):
                        marker_id = "marker_{}".format(j + 1)
                        self.assertEqual(
                            sample_geno[sample_id][marker_id],
                            genotype,
                        )
                        seen_markers.add(marker_id)

                    self.assertEqual(set(mapping_info.keys()), seen_markers)

                    # We've seen this sample now
                    seen_samples.add(sample_id)
            expected = {
                "sample_{}".format(i + sample_id_add) for i in
                range(nb_samples)
            }
            self.assertEqual(expected, seen_samples)

    def test_convert_beeline_2(self):
        """Tests the 'convert_beeline' function (different nb SNPs kw)."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        sample_genotype = defaultdict(dict)
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)
                    sample_genotype[sample_id][marker_id] = genotype

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Creating a temporary file
        tmp_filename_2 = None
        sample_genotype_2 = defaultdict(dict)
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename_2 = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + nb_samples + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)
                    sample_genotype_2[sample_id][marker_id] = genotype

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num SNPs"
        beeline2plink.convert_beeline(
            i_filenames=[tmp_filename, tmp_filename_2],
            out_dir=self.tmp_dir,
            locations=mapping_info,
            other_opts=other_options,
        )

        # Checking the two map files
        for map_filename in (tmp_filename, tmp_filename_2):
            map_filename = os.path.splitext(map_filename)[0] + ".map"
            self.assertTrue(os.path.isfile(map_filename))

            # Checking the first map file content
            seen_markers = set()
            with open(map_filename, "r") as i_file:
                for i, line in enumerate(i_file):
                    # Gathering the information
                    marker = "marker_{}".format(i + 1)
                    row = line.rstrip("\n").split("\t")

                    # Comparing the content of the file
                    self.assertTrue(marker in mapping_info)
                    self.assertEqual(4, len(row))
                    self.assertEqual(mapping_info[marker].chrom, int(row[0]))
                    self.assertEqual(marker, row[1])
                    self.assertEqual("0", row[2])
                    self.assertEqual(mapping_info[marker].pos, int(row[3]))

                    # We have compared this marker
                    seen_markers.add(marker)
            self.assertEqual(set(mapping_info.keys()), seen_markers)

        # Checking the two ped files
        zipped = zip(
            (1, nb_samples + 1),
            (tmp_filename, tmp_filename_2),
            (sample_genotype, sample_genotype_2),
        )
        for sample_id_add, ped_filename, sample_geno in zipped:
            ped_filename = os.path.splitext(ped_filename)[0] + ".ped"
            self.assertTrue(os.path.isfile(ped_filename))

            # Checking the ped file content
            seen_samples = set()
            with open(ped_filename, "r") as i_file:
                for i, line in enumerate(i_file):
                    # Gathering the information
                    sample_id = "sample_{}".format(i + sample_id_add)
                    row = line.rstrip("\n").split("\t")

                    # Checking the sample information
                    self.assertEqual(6 + nb_markers, len(row))
                    sample_info = row[:6]
                    self.assertEqual(
                        [sample_id, sample_id, "0", "0", "0", "-9"],
                        sample_info,
                    )

                    # Checking the genotypes
                    seen_markers = set()
                    genotypes = row[6:]
                    self.assertEqual(nb_markers, len(genotypes))
                    for j, genotype in enumerate(genotypes):
                        marker_id = "marker_{}".format(j + 1)
                        self.assertEqual(
                            sample_geno[sample_id][marker_id],
                            genotype,
                        )
                        seen_markers.add(marker_id)

                    self.assertEqual(set(mapping_info.keys()), seen_markers)

                    # We've seen this sample now
                    seen_samples.add(sample_id)
            expected = {
                "sample_{}".format(i + sample_id_add) for i in
                range(nb_samples)
            }
            self.assertEqual(expected, seen_samples)

    def test_convert_beeline_error_1(self):
        """Tests the 'convert_beeline' function with missing nb markers."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.convert_beeline(
                i_filenames=[tmp_filename] * 2,
                out_dir=self.tmp_dir,
                locations=mapping_info,
                other_opts=other_options,
            )
        self.assertEqual(
            tmp_filename + ": invalid header (missing 'Num Used SNPs' value)",
            e.exception.message,
        )

    def test_convert_beeline_error_2(self):
        """Tests the 'convert_beeline' function with missing marker map."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        sample_genotype = defaultdict(dict)
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)
                    sample_genotype[sample_id][marker_id] = genotype

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            if i + 1 == 3:
                continue
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function (should raise a warning)
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        with self._my_compatibility_assertLogs(level="WARNING") as cm:
            beeline2plink.convert_beeline(
                i_filenames=[tmp_filename],
                out_dir=self.tmp_dir,
                locations=mapping_info,
                other_opts=other_options,
            )
        self.assertEqual(1, len(cm.output))
        self.assertEqual(
            "WARNING:root:marker_3: no mapping information",
            cm.output[0],
        )

        # Checking the map file
        map_filename = os.path.splitext(tmp_filename)[0] + ".map"
        self.assertTrue(os.path.isfile(map_filename))

        # Checking the map file content
        seen_markers = set()
        with open(map_filename, "r") as i_file:
            for i, line in enumerate(i_file):
                # Gathering the information
                marker = "marker_{}".format(i + 1)
                row = line.rstrip("\n").split("\t")

                # Comparing the content of the file for all markers
                self.assertEqual(4, len(row))
                self.assertEqual(marker, row[1])
                self.assertEqual("0", row[2])

                # Comparing the content of the file specific for each marker
                if marker != "marker_3":
                    self.assertTrue(marker in mapping_info)
                    self.assertEqual(mapping_info[marker].chrom, int(row[0]))
                    self.assertEqual(mapping_info[marker].pos, int(row[3]))

                else:
                    self.assertFalse(marker in mapping_info)
                    self.assertEqual(0, int(row[0]))
                    self.assertEqual(0, int(row[3]))

                # We have compared this marker
                seen_markers.add(marker)
        self.assertEqual(
            set(mapping_info.keys()) | {"marker_3"},
            seen_markers,
        )

        # Checking the ped file
        ped_filename = os.path.splitext(tmp_filename)[0] + ".ped"
        self.assertTrue(os.path.isfile(ped_filename))

        # Checking the ped file content
        seen_samples = set()
        with open(ped_filename, "r") as i_file:
            for i, line in enumerate(i_file):
                # Gathering the information
                sample_id = "sample_{}".format(i + 1)
                row = line.rstrip("\n").split("\t")

                # Checking the sample information
                self.assertEqual(6 + nb_markers, len(row))
                sample_info = row[:6]
                self.assertEqual(
                    [sample_id, sample_id, "0", "0", "0", "-9"],
                    sample_info,
                )

                # Checking the genotypes
                seen_markers = set()
                genotypes = row[6:]
                self.assertEqual(nb_markers, len(genotypes))
                for j, genotype in enumerate(genotypes):
                    marker_id = "marker_{}".format(j + 1)
                    self.assertEqual(
                        sample_genotype[sample_id][marker_id],
                        genotype,
                    )
                    seen_markers.add(marker_id)

                self.assertEqual(
                    set(mapping_info.keys()) | {"marker_3"},
                    seen_markers,
                )

                # We've seen this sample now
                seen_samples.add(sample_id)
        expected = {
            "sample_{}".format(i + 1) for i in
            range(nb_samples)
        }
        self.assertEqual(expected, seen_samples)

    def test_convert_beeline_error_3(self):
        """Tests the 'convert_beeline' function (wrong nb of markers)."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers + 1), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.convert_beeline(
                i_filenames=[tmp_filename] * 2,
                out_dir=self.tmp_dir,
                locations=mapping_info,
                other_opts=other_options,
            )
        self.assertEqual(
            tmp_filename + ": missing 1 marker for sample 'sample_1'",
            e.exception.message,
        )

    def test_convert_beeline_error_4(self):
        """Tests the 'convert_beeline' function (missing column)."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "B Allele Freq", "Log R Ratio", sep=",",
                  file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.convert_beeline(
                i_filenames=[tmp_filename] * 2,
                out_dir=self.tmp_dir,
                locations=mapping_info,
                other_opts=other_options,
            )
        self.assertEqual(
            tmp_filename + ": 'Allele2 - Forward': missing column",
            e.exception.message,
        )

    def test_convert_beeline_error_5(self):
        """Tests the 'convert_beeline' function (invalid order)."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for sample in range(nb_samples):
                sample_id = "sample_{}".format(sample + 1)

                for marker in range(nb_markers):
                    marker_id = "marker_{}".format(marker + 1)

                    # We want to switch markers for sample 2
                    if sample_id == "sample_2":
                        if marker_id == "marker_3":
                            marker_id = "marker_4"
                        elif marker_id == "marker_4":
                            marker_id = "marker_3"

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.convert_beeline(
                i_filenames=[tmp_filename] * 2,
                out_dir=self.tmp_dir,
                locations=mapping_info,
                other_opts=other_options,
            )
        self.assertEqual(
            tmp_filename + ": marker order is not the same for sample "
                           "'sample_2'",
            e.exception.message,
        )

    def test_convert_beeline_error_6(self):
        """Tests the 'convert_beeline' function (invalid sorting)."""
        # The number of samples and of markers for this test
        nb_samples = 3
        nb_markers = 10

        # Creating a temporary file
        tmp_filename = None
        with NamedTemporaryFile("w", dir=self.tmp_dir, delete=False,
                                suffix=".csv") as f:
            tmp_filename = f.name

            # We need a header line
            print("[Header]", file=f)
            print("Some information", file=f)
            print("Num Used Samples,3", file=f)
            print("Some more information", file=f)
            print("Num Used SNPs,{}".format(nb_markers), file=f)
            print("Some more information", file=f)
            print("Some final information", file=f)

            # Needs consistent alleles for 10 markers
            alleles = {}
            for marker in range(nb_markers):
                alleles["marker_{}".format(marker + 1)] = tuple(
                    random.sample(("A", "C", "T", "G"), 2)
                )

            # We write the data
            print("[Data]", file=f)
            print("Sample ID", "SNP Name", "X", "Y",
                  "Allele1 - Forward", "Allele2 - Forward",
                  "B Allele Freq", "Log R Ratio", sep=",", file=f)
            for marker in range(nb_markers):
                marker_id = "marker_{}".format(marker + 1)

                for sample in range(nb_samples):
                    sample_id = "sample_{}".format(sample + 1)

                    # Getting the possible alleles
                    missing = random.random() < 0.1
                    marker_alleles = alleles[marker_id]
                    a1 = "-" if missing else random.choice(marker_alleles)
                    a2 = "-" if missing else random.choice(marker_alleles)
                    genotype = "0 0" if a1 == "-" else "{} {}".format(a1, a2)

                    # Printing the file
                    print(sample_id, marker_id, random.uniform(0, 3),
                          random.uniform(0, 3), a1, a2, random.random(),
                          random.uniform(-10, 10), sep=",", file=f,)

        # Generating mapping information
        mapping_info = {}
        for i in range(nb_markers):
            mapping_info["marker_{}".format(i + 1)] = beeline2plink._Location(
                chrom=random.randint(1, 26),
                pos=random.randint(1, 1000000),
            )

        # Executing the function
        other_options = _DummyArgs()
        other_options.nb_snps_kw = "Num Used SNPs"
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.convert_beeline(
                i_filenames=[tmp_filename] * 2,
                out_dir=self.tmp_dir,
                locations=mapping_info,
                other_opts=other_options,
            )
        self.assertEqual(
            tmp_filename + ": data should be sorted by samples, not by "
                           "markers ('sample_1' had 1 marker, expecting 10",
            e.exception.message,
        )

    def test_check_args(self):
        """Tests the 'check_args' function."""
        # Creating dummy Beeline reports
        beeline_reports = [
            os.path.join(self.tmp_dir, "file_{}.csv".format(i + 1))
            for i in range(10)
        ]
        for filename in beeline_reports:
            with open(filename, "w") as o_file:
                pass

        # Creating a dummy map file
        map_filename = os.path.join(self.tmp_dir, "map_file.csv")
        with open(map_filename, "w") as o_file:
            print(
                "Illumina, Inc.\n"
                "[Heading]\n"
                "Descriptor File Name,HumanOmni25Exome-8v1-1_A.bpm\n"
                "Assay Format,Infinium LCG\n"
                "Date Manufactured,4/22/2014\n"
                "Loci Count ,2583651\n"
                "[Assay]\n"
                "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,"
                "AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,"
                "Species,Source,SourceVersion,SourceStrand,SourceSeq,"
                "TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand\n"
                "Dummy_data",
                file=o_file,
            )

        # Creating dummy options
        args = _DummyArgs()
        args.i_filenames = beeline_reports
        args.map_filename = map_filename
        args.delim = ","
        args.id_col = "Name"
        args.chr_col = "Chr"
        args.pos_col = "MapInfo"
        args.output_dir = self.tmp_dir
        args.nb_snps_kw = "Num Used SNPs"

        # Executing the function
        beeline2plink.check_args(args)

    def test_check_args_error_1(self):
        """Tests the 'check_args' function (missing beeline report(s))."""
        # Creating dummy Beeline reports (missing the third one)
        beeline_reports = [
            os.path.join(self.tmp_dir, "file_{}.csv".format(i + 1))
            for i in range(10)
        ]
        for filename in beeline_reports:
            if not filename.endswith("file_3.csv"):
                with open(filename, "w") as o_file:
                    pass

        # Creating a dummy map file
        map_filename = os.path.join(self.tmp_dir, "map_file.csv")
        with open(map_filename, "w") as o_file:
            print(
                "Illumina, Inc.\n"
                "[Heading]\n"
                "Descriptor File Name,HumanOmni25Exome-8v1-1_A.bpm\n"
                "Assay Format,Infinium LCG\n"
                "Date Manufactured,4/22/2014\n"
                "Loci Count ,2583651\n"
                "[Assay]\n"
                "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,"
                "AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,"
                "Species,Source,SourceVersion,SourceStrand,SourceSeq,"
                "TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand\n"
                "Dummy_data",
                file=o_file,
            )

        # Creating dummy options
        args = _DummyArgs()
        args.i_filenames = beeline_reports
        args.map_filename = map_filename
        args.delim = ","
        args.id_col = "Name"
        args.chr_col = "Chr"
        args.pos_col = "MapInfo"
        args.output_dir = self.tmp_dir
        args.nb_snps_kw = "Num Used SNPs"

        # Executing the function
        self.assertFalse(os.path.isfile(beeline_reports[2]))
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.check_args(args)
        self.assertEqual(
            beeline_reports[2] + ": no such file",
            e.exception.message,
        )

    def test_check_args_error_2(self):
        """Tests the 'check_args' function (missing map file)."""
        # Creating dummy Beeline reports
        beeline_reports = [
            os.path.join(self.tmp_dir, "file_{}.csv".format(i + 1))
            for i in range(10)
        ]
        for filename in beeline_reports:
            with open(filename, "w") as o_file:
                pass

        # Not creating a dummy map file
        map_filename = os.path.join(self.tmp_dir, "map_file.csv")

        # Creating dummy options
        args = _DummyArgs()
        args.i_filenames = beeline_reports
        args.map_filename = map_filename
        args.delim = ","
        args.id_col = "Name"
        args.chr_col = "Chr"
        args.pos_col = "MapInfo"
        args.output_dir = self.tmp_dir
        args.nb_snps_kw = "Num Used SNPs"

        # Executing the function
        self.assertFalse(os.path.isfile(map_filename))
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.check_args(args)
        self.assertEqual(
            map_filename + ": no such file",
            e.exception.message,
        )

    def test_check_args_error_3(self):
        """Tests the 'check_args' function (missing column in map file)."""
        # Creating dummy Beeline reports
        beeline_reports = [
            os.path.join(self.tmp_dir, "file_{}.csv".format(i + 1))
            for i in range(10)
        ]
        for filename in beeline_reports:
            with open(filename, "w") as o_file:
                pass

        # Creating a dummy map file
        map_filename = os.path.join(self.tmp_dir, "map_file.csv")
        with open(map_filename, "w") as o_file:
            print(
                "Illumina, Inc.\n"
                "[Heading]\n"
                "Descriptor File Name,HumanOmni25Exome-8v1-1_A.bpm\n"
                "Assay Format,Infinium LCG\n"
                "Date Manufactured,4/22/2014\n"
                "Loci Count ,2583651\n"
                "[Assay]\n"
                "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,"
                "AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,"
                "Species,Source,SourceVersion,SourceStrand,SourceSeq,"
                "TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand\n"
                "Dummy_data",
                file=o_file,
            )

        # Creating dummy options
        args = _DummyArgs()
        args.i_filenames = beeline_reports
        args.map_filename = map_filename
        args.delim = ","
        args.id_col = "Name"
        args.chr_col = "Chromosome"
        args.pos_col = "MapInfo"
        args.output_dir = self.tmp_dir
        args.nb_snps_kw = "Num Used SNPs"

        # Executing the function
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.check_args(args)
        self.assertEqual(
            map_filename + ": missing column 'Chromosome'",
            e.exception.message,
        )

    def test_check_args_error_4(self):
        """Tests the 'check_args' function (missing output directory)."""
        # Creating dummy Beeline reports
        beeline_reports = [
            os.path.join(self.tmp_dir, "file_{}.csv".format(i + 1))
            for i in range(10)
        ]
        for filename in beeline_reports:
            with open(filename, "w") as o_file:
                pass

        # Creating a dummy map file
        map_filename = os.path.join(self.tmp_dir, "map_file.csv")
        with open(map_filename, "w") as o_file:
            print(
                "Illumina, Inc.\n"
                "[Heading]\n"
                "Descriptor File Name,HumanOmni25Exome-8v1-1_A.bpm\n"
                "Assay Format,Infinium LCG\n"
                "Date Manufactured,4/22/2014\n"
                "Loci Count ,2583651\n"
                "[Assay]\n"
                "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,"
                "AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,"
                "Species,Source,SourceVersion,SourceStrand,SourceSeq,"
                "TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand\n"
                "Dummy_data",
                file=o_file,
            )

        # The missing directory
        missing_directory = os.path.join(self.tmp_dir, "missing_dir")

        # Creating dummy options
        args = _DummyArgs()
        args.i_filenames = beeline_reports
        args.map_filename = map_filename
        args.delim = ","
        args.id_col = "Name"
        args.chr_col = "Chr"
        args.pos_col = "MapInfo"
        args.output_dir = missing_directory
        args.nb_snps_kw = "Num Used SNPs"

        # Executing the function
        self.assertFalse(os.path.isdir(missing_directory))
        with self.assertRaises(beeline2plink.ProgramError) as e:
            beeline2plink.check_args(args)
        self.assertEqual(
            missing_directory + ": no such directory",
            e.exception.message,
        )

    @unittest.skipIf(platform.system() == "Windows",
                     "Not a problem on windows systems")
    def test_check_args_error_5(self):
        """Tests the 'check_args' function (output directory not writable)."""
        # Creating dummy Beeline reports
        beeline_reports = [
            os.path.join(self.tmp_dir, "file_{}.csv".format(i + 1))
            for i in range(10)
        ]
        for filename in beeline_reports:
            with open(filename, "w") as o_file:
                pass

        # Creating a dummy map file
        map_filename = os.path.join(self.tmp_dir, "map_file.csv")
        with open(map_filename, "w") as o_file:
            print(
                "Illumina, Inc.\n"
                "[Heading]\n"
                "Descriptor File Name,HumanOmni25Exome-8v1-1_A.bpm\n"
                "Assay Format,Infinium LCG\n"
                "Date Manufactured,4/22/2014\n"
                "Loci Count ,2583651\n"
                "[Assay]\n"
                "IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,"
                "AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,"
                "Species,Source,SourceVersion,SourceStrand,SourceSeq,"
                "TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand\n"
                "Dummy_data",
                file=o_file,
            )

        # The output directory
        output_directory = os.path.join(self.tmp_dir, "output_dir")
        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

        args = _DummyArgs()
        args.i_filenames = beeline_reports
        args.map_filename = map_filename
        args.delim = ","
        args.id_col = "Name"
        args.chr_col = "Chr"
        args.pos_col = "MapInfo"
        args.output_dir = output_directory
        args.nb_snps_kw = "Num Used SNPs"

        # Executing the function
        try:
            # Changing the permission of the directory
            os.chmod(output_directory, 0o111)

            # Checking the arguments
            with self.assertRaises(beeline2plink.ProgramError) as e:
                beeline2plink.check_args(args)
            self.assertEqual(
                output_directory + ": not writable",
                e.exception.message,
            )

        finally:
            # Changing the permission back
            os.chmod(output_directory, 0o750)


class BaseTestCaseContext_Compatibility:

    def __init__(self, test_case):
        self.test_case = test_case

    def _raiseFailure(self, standardMsg):
        msg = self.test_case._formatMessage(self.msg, standardMsg)
        raise self.test_case.failureException(msg)


_LoggingWatcher = collections.namedtuple("_LoggingWatcher",
                                         ["records", "output"])


class CapturingHandler_Compatibility(logging.Handler):
    """
    A logging handler capturing all (raw and formatted) logging output.
    """

    def __init__(self):
        logging.Handler.__init__(self)
        self.watcher = _LoggingWatcher([], [])

    def flush(self):
        pass

    def emit(self, record):
        self.watcher.records.append(record)
        msg = self.format(record)
        self.watcher.output.append(msg)


class AssertLogsContext_Compatibility(BaseTestCaseContext_Compatibility):
    """A context manager used to implement TestCase.assertLogs()."""

    LOGGING_FORMAT = "%(levelname)s:%(name)s:%(message)s"

    def __init__(self, test_case, logger_name, level):
        BaseTestCaseContext_Compatibility.__init__(self, test_case)
        self.logger_name = logger_name
        if level:
            # Python < 3.4 logging doesn't have a _nameToLevel dictionary
            nameToLevel = {
                'CRITICAL': logging.CRITICAL,
                'ERROR': logging.ERROR,
                'WARN': logging.WARNING,
                'WARNING': logging.WARNING,
                'INFO': logging.INFO,
                'DEBUG': logging.DEBUG,
                'NOTSET': logging.NOTSET,
            }
            self.level = nameToLevel.get(level, level)
        else:
            self.level = logging.INFO
        self.msg = None

    def __enter__(self):
        if isinstance(self.logger_name, logging.Logger):
            logger = self.logger = self.logger_name
        else:
            logger = self.logger = logging.getLogger(self.logger_name)
        formatter = logging.Formatter(self.LOGGING_FORMAT)
        handler = CapturingHandler_Compatibility()
        handler.setFormatter(formatter)
        self.watcher = handler.watcher
        self.old_handlers = logger.handlers[:]
        self.old_level = logger.level
        self.old_propagate = logger.propagate
        logger.handlers = [handler]
        logger.setLevel(self.level)
        logger.propagate = False
        return handler.watcher

    def __exit__(self, exc_type, exc_value, tb):
        self.logger.handlers = self.old_handlers
        self.logger.propagate = self.old_propagate
        self.logger.setLevel(self.old_level)
        if exc_type is not None:
            # let unexpected exceptions pass through
            return False
        if len(self.watcher.records) == 0:
            self._raiseFailure(
                "no logs of level {} or higher triggered on {}"
                .format(logging.getLevelName(self.level), self.logger.name))


if __name__ == "__main__":
    unittest.main()
