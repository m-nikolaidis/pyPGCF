import os
import unittest
import zipfile
from pathlib import Path
from unittest.mock import call, patch

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from context import config, download_genomes, utils

FILETYPE_TO_DIR = {
    "cds": "CDS_fasta_files",
    "protein": "Protein_fasta_files",
    "genome": "Genomic_fasta_files",
}


def calculate_num_records_in_fasta(file: Path) -> int:
    num = 0
    parser = parse(file, "fasta")
    for _ in parser:
        num += 1
    return num


class TestModule(unittest.TestCase):
    def setUp(self):
        self.temp_dir = Path("/tmp/genome_downloader_test")
        taxon = 1234
        assembly_source = "RefSeq"
        assembly_level = "chromosome,complete"
        debug = True
        self.gd = download_genomes.GenomeDownloader(
            taxon=taxon,
            out_dir=self.temp_dir,
            assembly_source=assembly_source,
            assembly_level=assembly_level,
            debug=debug,
        )
        os.makedirs(self.temp_dir, exist_ok=True)

    def tearDown(self):
        for root, dirs, files in os.walk(self.temp_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(self.temp_dir)

    def test_download_hydrated(self):
        res = self.gd.download_hydrated()
        self.assertIs(res, 0)

    @pytest.mark.skip(reason="Expects existing file")
    def test_extract_dataset_zip(self):
        expected_directory = Path(self.gd.out_dir / "ncbi_dataset")
        readme_file = Path(self.gd.out_dir / "README.md")
        self.assertIs(expected_directory.exists(), False)
        self.gd.extract_dataset_zip()
        self.assertIs(expected_directory.exists(), True)
        utils.recursive_unlink(expected_directory)
        utils.recursive_unlink(readme_file)

    @patch("context.utils.check_if_file_exists")
    @patch("zipfile.ZipFile.extractall")
    def test_extract_dataset_zip_mock(self, mock_extractall, mock_check_file_exists):
        with zipfile.ZipFile(self.gd.out_dir / "dataset.zip", "w") as _:
            pass
        self.gd.extract_dataset_zip()
        mock_extractall.assert_called_once_with(self.gd.out_dir)
        mock_check_file_exists.assert_called_once()

    @patch("context.utils.check_if_file_exists")
    def test_extract_dataset_zip_nonexistent(self, mock_check_file_exists):
        with self.assertRaises(FileNotFoundError):
            self.gd.extract_dataset_zip()
        mock_check_file_exists.assert_called_once()

    def test_create_output_directories(self):
        self.gd.create_output_directories()
        expected_directories = [
            self.temp_dir / "Genomic_fasta_files",
            self.temp_dir / "Protein_fasta_files",
            self.temp_dir / "CDS_fasta_files",
        ]
        unexpected_directory = self.temp_dir / "GBFF_files"
        for expected_directory in expected_directories:
            self.assertIs(expected_directory.exists(), True)
            expected_directory.rmdir()
        self.assertIs(unexpected_directory.exists(), False)
        # TODO: When I add option to create either protein or CDS files I should change this test

    @patch("pathlib.Path.mkdir")
    def test_create_output_directories_mock(self, mock_mkdir):
        self.gd.create_output_directories()
        calls = [call(exist_ok=True, parents=True) for _ in FILETYPE_TO_DIR.values()]
        mock_mkdir.assert_has_calls(calls, any_order=True)

    @patch("pathlib.Path.glob")
    def test_save_gbff_dir_to_memory(self, mock_glob):
        mock_glob.return_value = [
            Path("data/GCF_022870945"),
            Path("data/GCF_902375805"),
        ]
        self.gd.save_gbff_dir_to_memory()
        self.assertIn("GCF_022870945", self.gd.gbff_files)
        self.assertIn(
            Path("data/GCF_022870945/genomic.gbff"), self.gd.gbff_files.values()
        )
        self.assertIs(type(self.gd.gbff_files), dict)
        self.assertIs(len(self.gd.gbff_files.keys()), 2)

    def test_process_gbff_files_with_plasmids(self):
        self.gd.gbff_files = {
            "GCF_022870945": config.TEST_DATA_DIR / "GCF_022870945.gbff",
            "GCF_902375805": config.TEST_DATA_DIR / "GCF_902375805.gbff",
        }
        self.gd.keep_plasmids = True
        self.gd.create_output_directories()
        self.gd.process_gbff_files()

        expected_directories = [
            self.gd.out_dir / "Genomic_fasta_files",
            self.gd.out_dir / "Protein_fasta_files",
            self.gd.out_dir / "CDS_fasta_files",
        ]
        suffixes = [".fna", ".faa", ".fna"]
        num_records = {
            "Genomic_fasta_files": {
                "GCF_022870945": 3,
                "GCF_902375805": 10,
            },
            "Protein_fasta_files": {
                "GCF_022870945": 2467,
                "GCF_902375805": 3563,
            },
            "CDS_fasta_files": {
                "GCF_022870945": 2467,
                "GCF_902375805": 3563,
            },
        }
        for expected_directory, suffix in zip(expected_directories, suffixes):
            g1 = "GCF_022870945"
            g2 = "GCF_902375805"
            f1 = expected_directory / (g1 + suffix)
            f2 = expected_directory / (g2 + suffix)
            self.assertIs(f1.exists(), True)
            self.assertIs(f2.exists(), True)
            num_f1_records = calculate_num_records_in_fasta(f1)
            num_f2_records = calculate_num_records_in_fasta(f2)
            exp_num_f1_records = num_records.get(expected_directory.stem, {}).get(g1)
            exp_num_f2_records = num_records.get(expected_directory.stem, {}).get(g2)
            if exp_num_f1_records is None or exp_num_f2_records is None:
                raise ValueError("Not possible")
            self.assertEqual(num_f1_records, exp_num_f1_records)
            self.assertEqual(num_f2_records, exp_num_f2_records)
            utils.recursive_unlink(expected_directory)

    def test_process_gbff_files_without_plasmids(self):
        self.gd.gbff_files = {
            "GCF_022870945": config.TEST_DATA_DIR / "GCF_022870945.gbff",
            "GCF_902375805": config.TEST_DATA_DIR / "GCF_902375805.gbff",
        }
        self.gd.keep_plasmids = False
        self.gd.create_output_directories()
        self.gd.process_gbff_files()

        expected_directories = [
            self.gd.out_dir / "Genomic_fasta_files",
            self.gd.out_dir / "Protein_fasta_files",
            self.gd.out_dir / "CDS_fasta_files",
        ]
        suffixes = [".fna", ".faa", ".fna"]
        num_records = {
            "Genomic_fasta_files": {
                "GCF_022870945": 1,
                "GCF_902375805": 10,
            },
            "Protein_fasta_files": {
                "GCF_022870945": 2440,
                "GCF_902375805": 3563,
            },
            "CDS_fasta_files": {
                "GCF_022870945": 2440,
                "GCF_902375805": 3563,
            },
        }
        for expected_directory, suffix in zip(expected_directories, suffixes):
            g1 = "GCF_022870945"
            g2 = "GCF_902375805"
            f1 = expected_directory / (g1 + suffix)
            f2 = expected_directory / (g2 + suffix)
            self.assertIs(f1.exists(), True)
            self.assertIs(f2.exists(), True)
            num_f1_records = calculate_num_records_in_fasta(f1)
            num_f2_records = calculate_num_records_in_fasta(f2)
            exp_num_f1_records = num_records.get(expected_directory.stem, {}).get(g1)
            exp_num_f2_records = num_records.get(expected_directory.stem, {}).get(g2)
            if exp_num_f1_records is None or exp_num_f2_records is None:
                raise ValueError("Not possible")
            self.assertEqual(num_f1_records, exp_num_f1_records)
            self.assertEqual(num_f2_records, exp_num_f2_records)
            utils.recursive_unlink(expected_directory)

    def test_calculate_GC_and_N_base_multiple_records(self):
        records = [
            SeqRecord(name="1", seq=Seq("AGCT")),
            SeqRecord(name="2", seq=Seq("CGAT")),
            SeqRecord(name="3", seq=Seq("NNNN")),
        ]
        result = self.gd.calculate_GC_and_N(records)
        self.assertEqual(result, (12, 33.33, 33.333))

    def test_calculate_GC_and_N_single_record(self):
        records = [SeqRecord(Seq("AAGG" * 25))]
        result = self.gd.calculate_GC_and_N(records)
        self.assertEqual(result, (100, 0.0, 50.0))

    def test_calculate_GC_N_no_sequence(self):
        records = []
        result = self.gd.calculate_GC_and_N(records)
        self.assertEqual(result, (0, 0.0, 0.0))

    def test_calculate_GC_and_N_all_known_nucleotides(self):
        records = [
            SeqRecord(name="1", seq=Seq("AAA")),
            SeqRecord(name="2", seq=Seq("CCC")),
            SeqRecord(name="3", seq=Seq("GGG")),
            SeqRecord(name="4", seq=Seq("TTT")),
        ]
        result = self.gd.calculate_GC_and_N(records)
        self.assertEqual(result, (12, 0.0, 50.0))

    def test_calculate_GC_and_N_all_unknown_nucleotides(self):
        records = [
            SeqRecord(name="1", seq=Seq("NNN")),
            SeqRecord(name="1", seq=Seq("NNNNNN")),
        ]
        result = self.gd.calculate_GC_and_N(records)
        self.assertEqual(result, (9, 100.0, 0.0))

    @patch("Bio.SeqIO.write")
    def test_gbff_file_to_genomic_fasta(self, mock_write):
        records = [SeqRecord(Seq("ATGC"))]
        self.gd.gbff_file_to_genomic_fasta(records, Path("/tmp/output.fna"))
        mock_write.assert_called_once_with(records, Path("/tmp/output.fna"), "fasta")

    @patch("Bio.SeqIO.write")
    def test_gbff_file_to_protein_fasta(self, mock_write):
        feature = SeqFeature(
            FeatureLocation(1, 12),
            type="CDS",
            qualifiers={
                "locus_tag": ["L1"],
                "translation": ["MTGK"],
                "product": ["protein"],
            },
        )
        records = [SeqRecord(Seq("ATGC"), features=[feature])]
        result = self.gd.gbff_file_to_protein_fasta(records, Path("/tmp/protein.faa"))
        self.assertEqual(result, 1)
        self.assertTrue(mock_write.called)

    def test_gather_gbff_source_info(self):
        feature = SeqFeature(
            FeatureLocation(1, 100), qualifiers={"organism": ["test_organism"]}
        )
        records = [SeqRecord(Seq("ATGC"), features=[feature]), SeqRecord(Seq("ATGC"))]
        result = self.gd.gather_gbff_source_info(records)
        self.assertIn("organism", result)

    def test_get_16S_from_gbff_real_data(self):
        self.gd.gbff_files = {
            "GCF_022870945": config.TEST_DATA_DIR / "GCF_022870945.gbff",
            "GCF_902375805": config.TEST_DATA_DIR / "GCF_902375805.gbff",
        }
        expected_sequences = {
            "GCF_022870945": ["LVJ81_RS08520", 1542],
            "GCF_902375805": ["FXY68_RS00045", 1538],
        }
        observed_sequences = {}
        for genome, gbff_file in self.gd.gbff_files.items():
            parser = parse(gbff_file, "genbank")
            records = [r for r in parser]
            rrna_16S_record = self.gd.get_16S_from_gbff(records)
            observed_sequences[genome] = [
                rrna_16S_record.name,
                len(rrna_16S_record.seq),
            ]
        for key in expected_sequences:
            exp_name = expected_sequences.get(key)[0]
            obs_name = observed_sequences.get(key)[0]
            exp_length = expected_sequences.get(key)[1]
            obs_length = observed_sequences.get(key)[1]
            self.assertEqual(exp_name, obs_name)
            self.assertEqual(exp_length, obs_length)

    def test_get_16S_from_gbff_mock_feature(self):
        feature = SeqFeature(
            FeatureLocation(1, 100),
            type="rRNA",
            qualifiers={"product": ["16S ribosomal RNA"], "locus_tag": ["L1"]},
        )
        records = [SeqRecord(Seq("ATGC" * 25), features=[feature])]
        result = self.gd.get_16S_from_gbff(records)
        self.assertEqual(result.id, "L1")
        self.assertEqual(len(result.seq), 99)

    @patch("pandas.DataFrame.to_excel")
    def test_write_annotations(self, mock_to_excel):
        self.gd.gbff_source_info = {"genome": {"dummy_key": "dummy_value"}}
        self.gd.write_annotations()
        self.assertTrue(mock_to_excel.called)

    def test_write_annotations_empty_input(self):
        self.gd.gbff_files = {
            "GCF_022870945": Path("../test_data/GCF_022870945.gbff"),
            "GCF_902375805": Path("../test_data/GCF_902375805.gbff"),
        }
        self.gd.write_annotations()
        expected_file = self.gd.out_dir / "Genome_information.xlsx"
        self.assertIs(expected_file.exists(), True)
        df = pd.read_excel(expected_file, index_col=0)
        self.assertIs(df.empty, True)

    @pytest.mark.skip(reason="Gives too many requests error")
    def test_pipeline(self):
        self.gd.download_hydrated()
        self.gd.extract_dataset_zip()
        self.gd.create_output_directories()
        self.gd.save_gbff_dir_to_memory()
        self.gd.process_gbff_files()
        self.gd.write_annotation()
        self.gd.write_16S_fasta()

    def test_calculate_perc_pseudogenes_has_pseudo(self):
        feature = SeqFeature(
            FeatureLocation(1, 100),
            type="CDS",
            qualifiers={"product": ["pseudo"], "pseudo": ["True"]},
        )
        records = [SeqRecord(Seq("ATGC" * 25), features=[feature])]
        result = self.gd.calculate_perc_pseudogenes(records)
        self.assertEqual(result, 100.0)

    def test_calculate_perc_pseudogenes_no_pseudo(self):
        feature = SeqFeature(
            FeatureLocation(1, 100),
            type="CDS",
            qualifiers={
                "product": ["Product1"],
            },
        )
        records = [SeqRecord(Seq("ATGC" * 25), features=[feature])]
        result = self.gd.calculate_perc_pseudogenes(records)
        self.assertEqual(result, 00.0)

    def test_calculate_perc_pseudogenes_empty_list(self):
        records = []
        result = self.gd.calculate_perc_pseudogenes(records)
        self.assertEqual(result, 100.0)

    def test_identify_plasmids_from_gbff(self):
        feature = SeqFeature(
            FeatureLocation(1, 100), qualifiers={"plasmid": ["test_plasmid"]}
        )
        records = [SeqRecord(Seq("ATGC"), features=[feature], description="plasmid")]
        result = self.gd.identify_plasmids_from_gbff(records)
        self.assertIn(records[0].name, result)

    def test_identify_plasmids_from_gbff_no_plasmid(self):
        feature = SeqFeature(FeatureLocation(1, 100), qualifiers={})
        records = [SeqRecord(Seq("ATGC"), features=[feature])]
        result = self.gd.identify_plasmids_from_gbff(records)
        self.assertEqual(0, len(result))

    def test_identify_plasmids_from_gbff_empty(self):
        records = []
        result = self.gd.identify_plasmids_from_gbff(records)
        self.assertEqual(0, len(result))

    @patch("Bio.SeqIO.write")
    def test_write_16S_fasta(self, mock_write):
        self.gd.rrna_16S_records = {"genome1": SeqRecord(Seq("ATGC"))}
        self.gd.write_16S_fasta()
        self.assertTrue(mock_write.called)
        self.assertTrue(self.gd.out_dir / "16S_sequences.fna")

    @patch("context.utils.recursive_unlink")
    def test_remove_dataset_archive(self, mock_recursive_unlink):
        self.gd.remove_dataset_archive()
        mock_recursive_unlink.assert_called()


#


# def test_write_annotations(self):
#     self.gd.gbff_files = {
#         "GCF_022870945": config.TEST_DATA_DIR / "GCF_022870945.gbff",
#         "GCF_902375805": config.TEST_DATA_DIR / "GCF_902375805.gbff",
#     }
#
#     self.gbff_source_info = {
#         "GCF_022870945": {
#             "organism": "Vitreoscilla stercoraria",
#             "mol_type": "genomic DNA",
#             "strain": "SAG 1488-6",
#             "isolation_source": "Dung",
#             "host": "Bos taurus",
#             "culture_collection": "DSM:513",
#             "type_material": "type strain of Vitreoscilla stercoraria",
#             "db_xref": "taxon:61",
#             "collection_date": "1990",
#             "Perc_pseudogenes": 0,
#             "Genome_len": 1000000,
#             "N%": 0,
#             "GC%": 55,
#         },
#         "GCF_902375805": {
#             "organism": "Vitreoscilla massiliensis",
#             "mol_type": "genomic DNA",
#             "isolate": "MGYG-HGUT-01520",
#             "isolation_source": "human gut",
#             "db_xref": "taxon:1689272",
#             "note": "contig: NZ_LN898212.1",
#             "Perc_pseudogenes": 0,
#             "Genome_len": 1000000,
#             "N%": 0,
#             "GC%": 55,
#         },
#     }

# @patch("Bio.SeqIO.write")
# def test_gbff_file_to_cds_fasta(self, mock_write):
#     feature = SeqFeature(
#         FeatureLocation(1, 100),
#         type="CDS",
#         qualifiers={
#             "locus_tag": ["L1"],
#             "translation": ["MTGK"],
#             "product": ["protein"],
#         },
#     )
#     records = [SeqRecord(Seq("ATGC"), features=[feature])]
#     result = self.gd.gbff_file_to_cds_fasta(records, Path("/tmp/cds.fna"))
#     self.assertEqual(result, 1)
#     self.assertTrue(mock_write.called)

# No longer needed
# @patch("context.utils.is_valid_assembly_source")
# def test_init_invalid_source(self, mock_is_valid_assembly_source):
#     mock_is_valid_assembly_source.return_value = False
#     with self.assertRaises(ValueError):
#         download_genomes.GenomeDownloader(taxon="test_taxon", out_dir=self.temp_dir)
