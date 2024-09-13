import unittest
from io import StringIO
from pathlib import Path
from tempfile import mkstemp
from unittest.mock import mock_open, patch, MagicMock

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from context import eggnog

# Mock data for testing
FASTA_DATA = """>protein1
MKTIIALSYIFCLVFAQ
>protein2
MDAIKKKMQLEDKVEELLSK
"""

EGGNOG_CSV_DATA = """#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n
protein1\t224308.BSU02610\t4.66e-175\t488\tCOG2304@1|root,COG2304@2|Bacteria,1UZ85@1239|Firmicutes,4HCUJ@91061|Bacilli,1ZE5Y@1386|Bacillus\t91061|Bacilli\tT\tvWA found in TerF C terminus\tycbR\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tTerD,vWA-TerF-like\n
protein2\t224308.BSU40690\t1.27e-99\t290\t2BWFV@1|root,32QWV@2|Bacteria,1V8UN@1239|Firmicutes,4HK47@91061|Bacilli,1ZHNT@1386|Bacillus\t91061|Bacilli\t-\t-\tyybC\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tDUF2798\n
"""

EGGNOG_HEADERS = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]

# Mock for COG_CATEGORY_ANNOTATION
COG_CATEGORY_ANNOTATION = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: Prophages, transposons",
    "R": "General function, prediction only",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}


######
class TestEggNOGRunner(unittest.TestCase):
    def setUp(self):
        self.runner = eggnog.EggNOGRunner(
            fasta_files_list=[Path("/fake/fasta1.fasta")],
            out_dir=Path("/fake/out_dir"),
            cores=4,
            pident=0.8,
            qcov=0.9,
            scov=0.85,
            input_type="nucl",
            database_dir=Path("/fake/db"),
            debug=False,
        )

    @patch("pathlib.Path.mkdir")
    def test_setup_directories(self, mock_mkdir):
        # Initialize EggNOGRunner instance

        # Test directory setup
        self.runner.setup_directories()

        # Check that the directory was created with the correct arguments
        self.assertEqual(mock_mkdir.call_count, 2)

    @patch("pathlib.Path.rename")
    @patch("pathlib.Path.with_suffix")
    @patch("pathlib.Path.exists", return_value=True)
    def test_clean_unnecessary_output(self, mock_exists, mock_with_suffix, mock_rename):
        # Mock files to remove
        mock_with_suffix.return_value = Path(
            "/fake/out_dir/eggNOG/eggnog_mapper/fasta1.csv.emapper.annotations"
        )

        # Test cleaning unnecessary files
        self.runner.clean_unessecary_output()

        # Check if with_suffix was called to remove unwanted files
        self.assertEqual(mock_with_suffix.call_count, 3)

        # Check if file renaming occurred
        mock_rename.assert_called_once_with(
            Path("/fake/out_dir/eggNOG/eggnog_mapper/fasta1_eggNOG_results.csv")
        )

    @pytest.mark.skip("Long running task. Not implemented properly")
    @patch("context.utils.multiprocess_dispatch")
    @patch("context.eggnog.EggNOGRunner.create_eggnog_cmd")
    @patch("context.eggnog.EggNOGRunner.setup_directories")
    @patch("context.eggnog.EggNOGRunner.clean_unessecary_output")
    def test_execute_eggnog_mapper(
        self, mock_clean_output, mock_setup_dirs, mock_create_cmd, mock_dispatch
    ):
        # Initialize EggNOGRunner instance

        # Mock the multiprocess_dispatch return value
        mock_dispatch.return_value = True

        # Run the full process
        self.runner.execute_eggnog_mapper()

        # Assert setup_directories, create_eggnog_cmd and clean_unessecary_output were called
        mock_setup_dirs.assert_called_once()
        mock_create_cmd.assert_called_once()
        mock_dispatch.assert_called_once()
        mock_clean_output.assert_called_once()

    def test_create_eggnog_cmd_nucl(self):
        # Use nucl as input type
        self.runner.input_type = "cds"
        self.runner.nucl = True
        # Run the command creation method
        self.runner.create_eggnog_cmd()

        # Expected command for nucl type
        expected_cmd = " ".join(
            [
                "emapper.py",
                "--cpu 4",
                "--pident 0.8",
                "--query_cov 0.9",
                "--subject_cov 0.85",
                "-o /fake/out_dir/eggNOG/eggnog_mapper/fasta1_eggNOG_results.csv",
                "-i /fake/fasta1.fasta",
                "--data_dir /fake/db",
                "--itype CDS --translate",
                "> /dev/null",
            ]
        )

        self.assertEqual(self.runner.eggnog_cmds[0], expected_cmd)

    def test_create_eggnog_cmd_prot(self):
        # Initialize EggNOGRunner instance with prot input type
        # Run the command creation method
        self.runner.input_type = "prot"
        self.runner.nucl = False
        self.runner.create_eggnog_cmd()

        # Expected command for prot type
        expected_cmd = " ".join(
            [
                "emapper.py",
                "--cpu 4",
                "--pident 0.8",
                "--query_cov 0.9",
                "--subject_cov 0.85",
                "-o /fake/out_dir/eggNOG/eggnog_mapper/fasta1_eggNOG_results.csv",
                "-i /fake/fasta1.fasta",
                "--data_dir /fake/db",
                "> /dev/null",
            ]
        )

        self.assertEqual(self.runner.eggnog_cmds[0], expected_cmd)


class TestEggNOGParser(unittest.TestCase):
    def setUp(self):
        self.parser = eggnog.EggNOGParser(
            fasta_files_list=[Path("/fake/fasta_file.fasta")],
            core_protein_table_file_list=None,
            out_dir=Path("/fake/out_dir"),
            debug=True,
        )
        tmpf = Path(mkstemp()[1])
        self.tmpdf = pd.read_csv(StringIO(EGGNOG_CSV_DATA), sep="\t")
        self.tmpdf.to_csv(tmpf, sep="\t", index=False)
        self.tmpf = tmpf
        self.expected_df = self.tmpdf.copy()
        self.expected_df.columns = EGGNOG_HEADERS
        self.expected_df = self.expected_df.set_index("query")

    def tearDown(self):
        self.tmpf.unlink()

    @patch("builtins.open", new_callable=mock_open, read_data=FASTA_DATA)
    @patch("Bio.SeqIO.parse")
    def test_load_all_reference_proteins(self, mock_seqio_parse, mock_file):
        # Mock the parsed records from SeqIO.parse
        mock_seqio_parse.return_value = [
            SeqRecord(Seq("MKTIIALSYIFCLVFAQ"), id="protein1"),
            SeqRecord(Seq("MDAIKKKMQLEDKVEELLSK"), id="protein2"),
        ]

        ref_proteins = self.parser.load_all_reference_proteins(
            Path("/fake/fasta_file.fasta")
        )
        expected = {"protein1": ["S"], "protein2": ["S"]}
        self.assertEqual(ref_proteins, expected)

    @patch("pandas.read_csv")
    def test_load_eggnogmapper_csv(self, mock_read_csv):
        # Simulate reading a CSV

        # Test if _load_eggnogmapper_csv properly loads CSV
        with self.assertRaises(FileNotFoundError):
            result_df = self.parser._load_eggnogmapper_csv(
                Path("/fake/eggnog_results.csv")
            )
            self.assertFalse(result_df.empty)

        result_df = self.parser._load_eggnogmapper_csv(self.tmpf)
        self.assertEqual(self.expected_df.equals(result_df), True)

    @pytest.mark.skip(reason="Not implemented properly yet")
    @patch("pandas.read_csv")
    @patch("context.eggnog.EggNOGParser.load_all_reference_proteins")
    @patch("builtins.open", new_callable=mock_open, read_data=FASTA_DATA)
    @patch("Bio.SeqIO.parse")
    def test_parse_eggnog_raw_results(
        self, mock_seqio_parse, mock_file, mock_load_proteins, mock_read_csv
    ):
        # Mock the protein loader and CSV reading
        mock_load_proteins.return_value = {"protein1": ["S"], "protein2": ["S"]}

        mock_seqio_parse.return_value = [
            SeqRecord(Seq("MKTIIALSYIFCLVFAQ"), id="protein1"),
            SeqRecord(Seq("MDAIKKKMQLEDKVEELLSK"), id="protein2"),
        ]
        mock_read_csv.return_value = self.parser._load_eggnogmapper_csv(self.tmpf)

        self.parser.parse_eggnog_raw_results()  # TODO: This gives error

        expected_proteome_cog_data = {
            self.tmpf.stem: {"protein1": ["K"], "protein2": ["K"]}
        }

        self.assertEqual(self.parser.proteome_cog_data, expected_proteome_cog_data)

    def test_turn_category_to_list(self):
        self.assertEqual(self.parser.turn_category_to_list("-"), ["S"])
        self.assertEqual(self.parser.turn_category_to_list("SK"), ["S", "K"])
        self.assertEqual(self.parser.turn_category_to_list(""), [])
        with self.assertRaises(TypeError):
            self.assertEqual(self.parser.turn_category_to_list(None), [])

    def test_calculate_cog_category_counts_per_proteome(self):
        cog_data_dict = {
            "Proteome1": {"p1": ["S"]},
            "Proteome2": {"p1": ["K"], "p2": ["K", "L"]},
        }

        expected_counts = {
            "Proteome1": {c: 0 for c in COG_CATEGORY_ANNOTATION},
            "Proteome2": {c: 0 for c in COG_CATEGORY_ANNOTATION},
        }
        expected_counts["Proteome1"]["S"] = 1
        expected_counts["Proteome1"]["Total"] = 1  # Total proteins
        expected_counts["Proteome2"]["K"] = 2
        expected_counts["Proteome2"]["L"] = 1
        expected_counts["Proteome2"]["Total"] = 2  # Total proteins

        counts = self.parser.calculate_cog_category_counts_per_proteome(cog_data_dict)
        self.assertEqual(expected_counts, counts)

    @patch("pandas.ExcelWriter")
    def test_write_counts_dictionary_to_excel(self, mock_excel_writer):
        counts_dict = {
            "Proteome1": {"S": 1, "K": 0, "Total": 1},
            "Proteome2": {"S": 1, "K": 2, "Total": 2},
        }

        tmpexcel = Path(mkstemp(suffix=".xlsx")[1])
        # Test writing to Excel
        self.parser.write_counts_dictionary_to_excel(counts_dict, tmpexcel)

        # Check if excel file was written properly
        xlreader = pd.ExcelFile(tmpexcel)
        observed_counts_df = xlreader.parse(sheet_name="counts", index_col=0)
        expected_counts_df = pd.DataFrame.from_dict(counts_dict, orient="index")

        observed_percent_df = xlreader.parse(sheet_name="percent", index_col=0)
        expected_percent_df = (
            expected_counts_df.div(expected_counts_df["Total"], axis=0) * 100
        )
        expected_percent_df = expected_percent_df.astype(int)

        self.assertEqual(xlreader.sheet_names, ["counts", "percent"])
        self.assertEqual(expected_counts_df.equals(observed_counts_df), True)
        self.assertEqual(expected_percent_df.equals(observed_percent_df), True)
        tmpexcel.unlink()

    @pytest.mark.skip(
        reason="Can't figure out how to make it work in unit testing. Should add integration test"
    )
    @patch("pandas.read_excel")
    @patch("context.eggnog.EggNOGParser._get_protein_subsets_from_core_prot_df")
    @patch("pandas.ExcelFile")
    @patch("pathlib.Path")
    def test_read_core_protein_tables(
        self, mock_path, mock_excel, mock_get_protein_subsets, mock_read_excel
    ):
        # Mock excel file
        wb = mock_excel.book()
        wb.create_sheet("Sheet1")
        ws = wb["Sheet1"]
        ws.loc["p1", "Core"] = 1
        ws.loc["p2", "Core"] = 1
        ws.loc["p3", "Core"] = 0
        # Mock Excel data
        mock_df = MagicMock(return_value=wb)
        mock_read_excel.return_value = mock_df
        #
        # # Mock the protein subset data
        # mock_get_protein_subsets.return_value = (["p1", "p2"], ["p3"])
        #
        self.parser.core_protein_table_file_list = [mock_path, mock_path]
        #
        # # Mock proteome_cog_data to simulate parsed eggNOG data
        self.parser.proteome_cog_data = {
            "Proteome1": {"p1": ["L"], "p2": ["K"], "p3": ["T"]},
        }
        #
        self.parser.read_core_protein_tables()
        #
        # # Check if read_excel is called correctly
        # mock_read_excel.assert_called_with(
        #     Path("/fake/core_protein_table1.xlsx"), index_col=0
        # )
        #
        # # Ensure that the internal methods for processing protein subsets are called
        # mock_get_protein_subsets.assert_called_once_with(mock_df)
        #
        # # Check if the cog counts were calculated (we'll assume calculate_cog_category_counts_per_proteome is tested separately)
        # self.assertIn("Proteome1", parser.core_protein_cog_counts)
        # self.assertIn("Proteome1", parser.fingerprint_protein_cog_counts)

    def test_get_protein_subsets_from_core_prot_df(self):
        # Mock DataFrame
        df = pd.DataFrame(
            {"core": [1, 1, 0], "fingerprint": [1, 0, 0]}, index=["p1", "p2", "p3"]
        )

        core_proteins, fingerprints = (
            self.parser._get_protein_subsets_from_core_prot_df(df)
        )

        # Check the correct identification of core proteins and fingerprints
        self.assertEqual(core_proteins, ["p1", "p2"])
        self.assertEqual(fingerprints, ["p1"])

    @patch("context.eggnog.EggNOGParser.perform_hypergeom_test")
    @patch("context.utils.dict_to_dataframe")
    def test_compare_core_and_fingerprints_against_background(
        self, mock_dict_to_df, mock_perform_hypergeom_test
    ):
        # Need to mock the output of perform_hypergeom_test
        # These are not correct hypergeometric output values
        mock_perform_hypergeom_test.return_value = (0.5, 1)

        # Mock dict_to_dataframe return
        mock_df = MagicMock()
        mock_dict_to_df.return_value = mock_df

        # Mock data for proteome_cog_counts and core_protein_cog_counts
        self.parser.proteome_cog_counts = {
            "Proteome1": {"Total": 100, "L": 20, "K": 30}
        }
        self.parser.core_protein_cog_counts = {
            "Proteome1": {"Total": 50, "L": 10, "K": 15}
        }
        self.parser.fingerprint_protein_cog_counts = {
            "Proteome1": {"Total": 10, "L": 3, "K": 9}
        }

        result = self.parser.compare_core_and_fingerprints_against_background()

        # Check that the result is returned as expected
        self.assertIn("core", result)
        self.assertIn("fingerprint", result)

        # Ensure perform_hypergeom_test was called with the correct parameters
        mock_perform_hypergeom_test.assert_any_call(
            population_size=100,
            population_successes=20,
            sample_size=50,
            sample_successes=10,
        )

        mock_perform_hypergeom_test.assert_any_call(
            population_size=100,
            population_successes=30,
            sample_size=50,
            sample_successes=15,
        )

    def test_compare_core_and_fingerprints_against_background_no_fingerprints(self):
        # Mock data for proteome_cog_counts and core_protein_cog_counts
        self.parser.proteome_cog_counts = {
            "Proteome1": {"Total": 100, "L": 20, "K": 30},
            "Proteome2": {"Total": 100, "L": 20, "K": 30},
        }
        self.parser.core_protein_cog_counts = {
            "Proteome1": {"Total": 50, "L": 10, "K": 15}
        }
        self.parser.fingerprint_protein_cog_counts = {}

        result = self.parser.compare_core_and_fingerprints_against_background()

        # Check that the result is returned as expected
        expected_core_result_df = pd.DataFrame.from_dict(
            {
                "Proteome1": {
                    "L_pvalue": 0.598436,
                    "L_fold_change": 1.0,
                    "K_pvalue": 0.586242,
                    "K_fold_change": 1.0,
                }
            },
            orient="index",
        ).round(3)
        observed_core_result_df = result["core"].round(3)
        print(observed_core_result_df)

        self.assertIn("core", result)
        self.assertNotIn("fingerprint", result)
        self.assertEqual(expected_core_result_df.equals(observed_core_result_df), True)

    def test_compare_core_and_fingerprints_against_background_no_fingerprints_many_proteomes(
        self,
    ):
        # Mock data for proteome_cog_counts and core_protein_cog_counts
        self.parser.proteome_cog_counts = {
            "Proteome1": {"Total": 100, "L": 20, "K": 30},
            "Proteome2": {"Total": 100, "L": 20, "K": 30},
        }
        self.parser.core_protein_cog_counts = {
            "Proteome1": {"Total": 50, "L": 10, "K": 15},
            "Proteome2": {"Total": 50, "L": 10, "K": 15},
        }
        self.parser.fingerprint_protein_cog_counts = {}

        result = self.parser.compare_core_and_fingerprints_against_background()

        # Check that the result is returned as expected
        expected_core_result_df = pd.DataFrame.from_dict(
            {
                "Proteome1": {
                    "L_pvalue": 0.598436,
                    "L_fold_change": 1.0,
                    "K_pvalue": 0.586242,
                    "K_fold_change": 1.0,
                },
                "Proteome2": {
                    "L_pvalue": 0.598436,
                    "L_fold_change": 1.0,
                    "K_pvalue": 0.586242,
                    "K_fold_change": 1.0,
                },
            },
            orient="index",
        ).round(3)
        observed_core_result_df = result["core"].round(3)

        self.assertIn("core", result)
        self.assertNotIn("fingerprint", result)
        self.assertEqual(expected_core_result_df.equals(observed_core_result_df), True)

    def test_highlight_pvalue_on_output(self):
        # Create a mock row with pvalues and fold changes
        row = pd.Series(
            {
                "L_pvalue": 0.01,
                "L_fold_change": 1.5,
                "K_pvalue": 0.05,
                "K_fold_change": 0.5,
            }
        )

        result = self.parser.highlight_pvalue_on_output(row)

        # Expected colorings
        expected = [
            "background-color:green",
            "background-color:green",
            "background-color:red",
            "background-color:red",
        ]

        self.assertEqual(result, expected)

    def test_highlight_pvalue_on_output_no_highlight(self):
        # Create a mock row with pvalues and fold changes
        row = pd.Series(
            {
                "L_pvalue": 0.1,
                "L_fold_change": 1.5,
                "K_pvalue": 0.5,
                "K_fold_change": 0.5,
            }
        )

        result = self.parser.highlight_pvalue_on_output(row)

        # Expected colorings
        expected = [""] * 4

        self.assertEqual(result, expected)

    def test_p_value_and_fold_change_calculated_correctly(self):
        # Calculate expected values
        population_size = 1000
        population_successes = 50
        sample_size = 300
        sample_successes = 2

        # Perform the hypergeometric test
        pval, fold_change = self.parser.perform_hypergeom_test(
            population_size=population_size,
            population_successes=population_successes,
            sample_size=sample_size,
            sample_successes=sample_successes,
        )

        expected_p_value = 2.9e-6
        self.assertAlmostEqual(pval, expected_p_value)

        # Check if fold change is calculated correctly
        expected_fold_change = (sample_successes / sample_size) / (
            population_successes / population_size
        )
        self.assertAlmostEqual(fold_change, expected_fold_change)

    def test_perform_hypergeom_test_invalid_sample_successes(self):
        # Calculate expected values
        population_size = 1000
        population_successes = 50
        sample_size = 0
        sample_successes = 1

        with self.assertRaises(ValueError):
            # Perform the hypergeometric test
            pval, fold_change = self.parser.perform_hypergeom_test(
                population_size=population_size,
                population_successes=population_successes,
                sample_size=sample_size,
                sample_successes=sample_successes,
            )

    def test_perform_hypergeom_test_succ_more_than_exp(self):
        population_size = 1000
        population_successes = 50

        # Successes are more than expected
        sample_size = 100
        sample_successes = 10
        pval, fold_change = self.parser.perform_hypergeom_test(
            population_size=population_size,
            population_successes=population_successes,
            sample_size=sample_size,
            sample_successes=sample_successes,
        )

        # Check if p-value is calculated correctly
        expected_p_value = 0.02144032938140247
        self.assertAlmostEqual(pval, expected_p_value)

    def test_perform_hypergeom_test_succ_less_than_exp(self):
        population_size = 1000
        population_successes = 50

        # Successes are more than expected
        sample_size = 100
        sample_successes = 4
        pval, fold_change = self.parser.perform_hypergeom_test(
            population_size=population_size,
            population_successes=population_successes,
            sample_size=sample_size,
            sample_successes=sample_successes,
        )

        # Check if p-value is calculated correctly
        expected_p_value = 0.42691543272377
        self.assertAlmostEqual(pval, expected_p_value)

    def test_fold_change_one(self):
        # Calculate expected values
        population_size = 1000
        population_successes = 50
        sample_size = 100
        sample_successes = 5

        # Perform the hypergeometric test
        pval, fold_change = self.parser.perform_hypergeom_test(
            population_size=population_size,
            population_successes=population_successes,
            sample_size=sample_size,
            sample_successes=sample_successes,
        )

        self.assertEqual(fold_change, 1)

    def test_perform_hypergeom_test_sample_size_larger(self):
        # Calculate expected values
        population_size = 1
        population_successes = 1
        sample_size = 3
        sample_successes = 2

        # Perform the hypergeometric test
        with self.assertRaises(ValueError):
            pval, fold_change = self.parser.perform_hypergeom_test(
                population_size=population_size,
                population_successes=population_successes,
                sample_size=sample_size,
                sample_successes=sample_successes,
            )

    def test_fold_change_zero(self):
        # Calculate expected values
        population_size = 100
        population_successes = 5
        sample_size = 20
        sample_successes = 0

        # Perform the hypergeometric test
        pval, fold_change = self.parser.perform_hypergeom_test(
            population_size=population_size,
            population_successes=population_successes,
            sample_size=sample_size,
            sample_successes=sample_successes,
        )

        self.assertEqual(fold_change, 0)
        self.assertAlmostEqual(pval, 0.3193094419898543)
