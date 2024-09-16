import unittest
from unittest.mock import patch
from pathlib import Path

# Assuming the module you provided is named cli
from context import cli


class TestValidateDirectory(unittest.TestCase):
    @patch("context.utils.check_if_dir_exists")
    @patch("context.utils.check_if_dir_is_empty")
    def test_validate_directory_exists_and_non_empty(self, mock_is_empty, mock_exists):
        # Mock behavior for existing and non-empty directory
        mock_exists.return_value = True
        mock_is_empty.return_value = False

        path = Path("/some/existing/directory")
        result = cli.validate_directory(path)

        # Check that validate_directory returns True
        self.assertTrue(result)
        mock_exists.assert_called_once_with(path)
        mock_is_empty.assert_called_once_with(path)

    @patch("context.utils.check_if_dir_exists")
    @patch("context.utils.check_if_dir_is_empty")
    def test_validate_directory_does_not_exist(self, mock_is_empty, mock_exists):
        # Mock behavior for non-existing directory
        mock_exists.return_value = False
        mock_is_empty.return_value = False

        path = Path("/non/existent/directory")
        result = cli.validate_directory(path)

        # Check that validate_directory returns False
        self.assertFalse(result)
        mock_exists.assert_called_once_with(path)
        mock_is_empty.assert_not_called()

    @patch("context.utils.check_if_dir_exists")
    @patch("context.utils.check_if_dir_is_empty")
    def test_validate_directory_empty(self, mock_is_empty, mock_exists):
        # Mock behavior for existing but empty directory
        mock_exists.return_value = True
        mock_is_empty.return_value = True

        path = Path("/empty/directory")
        result = cli.validate_directory(path)

        # Check that validate_directory returns False
        self.assertFalse(result)
        mock_exists.assert_called_once_with(path)
        mock_is_empty.assert_called_once_with(path)


class TestArgumentParsing(unittest.TestCase):
    def setUp(self):
        self.parser = cli.setup_parser()

    def test_add_species_demarcation_subparser(self):
        # Simulate CLI input for the species_demarcation module
        args = vars(
            self.parser.parse_args(
                [
                    "species_demarcation",
                    "-in",
                    "input_genomes_dir",
                    "-o",
                    "output_dir",
                    "--fastani_cores",
                    "4",
                    "--kmer",
                    "12",
                    "--fraglen",
                    "15000",
                    "--minfraction",
                    "0.3",
                    "--inflation",
                    "1.5",
                    "--mcl_cores",
                    "2",
                ]
            )
        )

        self.assertEqual(args["module"], "species_demarcation")
        self.assertEqual(args["in"], "input_genomes_dir")
        self.assertEqual(args["o"], "output_dir")
        self.assertEqual(args["fastani_cores"], 4)
        self.assertEqual(args["kmer"], 12)
        self.assertEqual(args["fraglen"], 15000)
        self.assertEqual(args["minfraction"], 0.3)
        self.assertEqual(args["inflation"], "1.5")
        self.assertEqual(args["mcl_cores"], 2)

    def test_add_orthologues_subparser(self):
        # Simulate CLI input for the orthologues module
        args = vars(
            self.parser.parse_args(
                [
                    "orthologues",
                    "-in",
                    "fasta_dir",
                    "-out",
                    "output_dir",
                    "-ref",
                    "reference_strain",
                    "--cores",
                    "8",
                    "--evalue",
                    "1e-5",
                    "--dmnd_sensitivity",
                    "sensitive",
                ]
            )
        )

        self.assertEqual(args["module"], "orthologues")
        self.assertEqual(args["in"], "fasta_dir")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["ref"], "reference_strain")
        self.assertEqual(args["cores"], 8)
        self.assertEqual(args["evalue"], "1e-5")
        self.assertEqual(args["dmnd_sensitivity"], "sensitive")

    def test_add_core_subparser(self):
        # Simulate CLI input for the core module
        args = vars(
            self.parser.parse_args(
                [
                    "core",
                    "-in",
                    "orthology_matrix",
                    "-out",
                    "output_dir",
                    "--core_perc",
                    "95",
                ]
            )
        )

        self.assertEqual(args["module"], "core")
        self.assertEqual(args["in"], "orthology_matrix")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["core_perc"], "95")

    def test_add_phylogenomic_subparser(self):
        # Simulate CLI input for the phylogenomic module
        args = vars(
            self.parser.parse_args(
                [
                    "phylogenomic",
                    "-in",
                    "orthology_matrix",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                ]
            )
        )

        self.assertEqual(args["module"], "phylogenomic")
        self.assertEqual(args["in"], "orthology_matrix")
        self.assertEqual(args["fasta_dir"], "fasta_dir")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["input_type"], "prot")
        self.assertNotEqual(args["cores"], 0)
        self.assertEqual(args["method"], "IQTree")
        self.assertEqual(args["no_keep_fasta"], False)
        self.assertEqual(args["tree_model"], "TEST")

    def test_add_phylogenomic_subparser_fasttree(self):
        # Simulate CLI input for the phylogenomic module
        args = vars(
            self.parser.parse_args(
                [
                    "phylogenomic",
                    "-in",
                    "orthology_matrix",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                    "--method",
                    "Fasttree",
                ]
            )
        )

        self.assertEqual(args["module"], "phylogenomic")
        self.assertEqual(args["in"], "orthology_matrix")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["input_type"], "prot")
        self.assertNotEqual(args["cores"], 0)
        self.assertEqual(args["method"], "Fasttree")
        self.assertEqual(args["no_keep_fasta"], False)
        self.assertEqual(args["tree_model"], "TEST")

    def test_add_eggnog_subparser(self):
        # Simulate CLI input for the phylogenomic module
        args = vars(
            self.parser.parse_args(
                [
                    "eggnog",
                    "-in",
                    "core.xlsx",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                    "--cores",
                    "10",
                ]
            )
        )

        self.assertEqual(args["module"], "eggnog")
        self.assertEqual(args["in"], "core.xlsx")
        self.assertEqual(args["fasta_dir"], "fasta_dir")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["input_type"], "prot")
        self.assertEqual(args["cores"], 10)
        self.assertEqual(args.get("method"), None)

    def test_add_eggnog_subparser_ommit_required(self):
        # Simulate CLI input for the phylogenomic module
        with self.assertRaises(SystemExit):
            _ = vars(
                self.parser.parse_args(
                    [
                        "eggnog",
                        "-in",
                        "core.xlsx",
                        "-out",
                        "output_dir",
                        "--cores",
                        "10",
                    ]
                )
            )

    def test_add_smbgc_subparser(self):
        # Simulate CLI input for the phylogenomic module
        args = vars(
            self.parser.parse_args(
                [
                    "smbgc",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                    "--cores",
                    "10",
                ]
            )
        )

        self.assertEqual(args["module"], "smbgc")
        self.assertEqual(args["fasta_dir"], "fasta_dir")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["cores"], 10)
        self.assertEqual(args.get("method"), None)

    def test_add_smbgc_subparser_noexistent_arg(self):
        # Simulate CLI input for the smbgc module
        with self.assertRaises(SystemExit):
            _ = vars(
                self.parser.parse_args(
                    [
                        "smbgc",
                        "-in",
                        "input_file",
                        "-out",
                        "output_dir",
                        "-fasta_dir",
                        "fasta_dir",
                        "--cores",
                        "10",
                    ]
                )
            )

    def test_add_download_subparser(self):
        # Simulate CLI input for the phylogenomic module
        args = vars(
            self.parser.parse_args(
                [
                    "download",
                    "-out",
                    "output_dir",
                    "-taxon",
                    "taxon1",
                    "--keep_plasmids",
                ]
            )
        )

        self.assertEqual(args["module"], "download")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["taxon"], "taxon1")
        self.assertEqual(args["keep_plasmids"], True)
        self.assertEqual(args.get("cores"), None)

    def test_add_virulence_subparser(self):
        # Simulate CLI input for the phylogenomic module
        args = vars(
            self.parser.parse_args(
                [
                    "virulence",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                ]
            )
        )

        self.assertEqual(args["module"], "virulence")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["fasta_dir"], "fasta_dir")

    def test_add_cazy_subparser(self):
        # Simulate CLI input for the cazy module
        args = vars(
            self.parser.parse_args(
                [
                    "cazy",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                ]
            )
        )

        self.assertEqual(args["module"], "cazy")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["fasta_dir"], "fasta_dir")
        self.assertEqual(args["db"], None)

    def test_add_amr_subparser(self):
        # Simulate CLI input for the amr module
        args = vars(
            self.parser.parse_args(
                [
                    "amr",
                    "-out",
                    "output_dir",
                    "-fasta_dir",
                    "fasta_dir",
                ]
            )
        )

        self.assertEqual(args["module"], "amr")
        self.assertEqual(args["out"], "output_dir")
        self.assertEqual(args["fasta_dir"], "fasta_dir")
