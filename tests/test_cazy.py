import unittest
from pathlib import Path
from unittest.mock import patch
from pandas import read_excel, DataFrame
from pandas.testing import assert_frame_equal

from context import cazy

# import sys
# sys.path.insert(0, f"{__file__}/../pypgcf")
# print(sys.path)


class TestCAZYBuilder(unittest.TestCase):
    def setUp(self):
        self.database_dir = Path("/tmp/fake/path/to/database")
        self.verbose = False
        self.cores = 4
        self.cazy_builder = cazy.CAZY_builder(
            self.database_dir, self.cores, self.verbose
        )

    @patch("context.cazy.execute_command")
    def test_setup_valid(self, mock_execute_command):
        mock_execute_command.return_value = 0
        self.cazy_builder.setup()
        build_cmd = f"dbcan_build --cores 4 --db-dir {self.database_dir}/CAZY --clean"
        mock_execute_command.assert_called_once_with(build_cmd)

    @patch("context.cazy.execute_command")
    def test_setup_invalid(self, mock_execute_command):
        mock_execute_command.side_effect = Exception("Command failed")
        with self.assertRaises(Exception):
            self.cazy_builder.setup()

    def test_validate_if_database_exists(self):
        # This is a placeholder test since the method is not implemented
        self.assertIsNone(self.cazy_builder.validate_if_database_exists())


class TestCAZYAnalyzer(unittest.TestCase):
    def setUp(self):
        self.cores = 4
        self.evalue = 1e-5
        self.database_dir = Path("/tmp/database")
        self.results_dir = Path("/tmp/results")
        self.fasta_dir = Path("/tmp/fasta")
        self.dmnd_sensitivity = "sensitive"
        self.protein = True
        self.verbose = False
        self.cazy_analyzer = cazy.CAZY_analyzer(
            self.cores,
            self.evalue,
            self.database_dir,
            self.results_dir,
            self.fasta_dir,
            self.dmnd_sensitivity,
            self.protein,
            self.verbose,
        )

    @patch("context.cazy.calc_avail_dispatchers", return_value=2)
    @patch("context.cazy.multiprocess_dispatch")
    @patch("context.cazy.create_temporary_file", return_value=Path("/fake/temp/file"))
    @patch("context.cazy.translate_fasta_records", return_value=["record1", "record2"])
    @patch("context.cazy.seqrecords_to_fasta")
    def test_execute_cazy_search_valid(
        self,
        mock_seqrecords_to_fasta,
        mock_translate_fasta_records,
        mock_create_temporary_file,
        mock_multiprocess_dispatch,
        mock_calc_avail_dispatchers,
    ):
        self.fasta_dir.mkdir(exist_ok=True)
        with open(self.fasta_dir / "test.fa", "w") as wf:
            wf.write(
                ">AAC07180\nMVALRNNWADISRKILEKRRFETKDLVKRLKVITGHDIHIQNYPVETPRVAFNPSIHVFENRLRIYARVV MGYYTYTSAIAEFDIDLEELYNPERKTYEANLTVLPNIKYDLWGVEDPRVYEIDGKLFMTYTGRTVNYFR TDIRTERTLPVTARYENGQWKKIAVFRMPEDIRSFVVSDKNAFLVKTDKLMLYHRLHMLNEKFYLAVCNV PEEVLYTNEFKEIEIGENITIMEEAPFETKIGWATPPVKVGEENLVLIHGVDKELTAYRVFAVLMNKEGY FTAVTPFYILEPKKIYEVYGDRPFVVFPCGIQRLENKLLISYGGADTVVVIGEIDLEELMNILYENRID"
            )
        self.cazy_analyzer.execute_cazy_search()
        mock_multiprocess_dispatch.assert_called_once()

    @patch("context.cazy.calc_avail_dispatchers", return_value=2)
    @patch("context.cazy.multiprocess_dispatch")
    @patch("context.cazy.create_temporary_file", return_value=Path("/fake/temp/file"))
    @patch("context.cazy.translate_fasta_records", return_value=["record1", "record2"])
    @patch("context.cazy.seqrecords_to_fasta")
    def test_execute_cazy_search_invalid(
        self,
        mock_seqrecords_to_fasta,
        mock_translate_fasta_records,
        mock_create_temporary_file,
        mock_multiprocess_dispatch,
        mock_calc_avail_dispatchers,
    ):
        self.cazy_analyzer.fasta_dir = Path("/invalid/path/to/fasta")
        self.assertRaises(FileNotFoundError, self.cazy_analyzer.execute_cazy_search)

    # self.cazy_analyzer.execute_cazy_search()

    def test_clean_hmm_output_valid(self):
        input_str = "GH5_1(1-300)+GH10_2(301-600)"
        expected_output = "GH"
        result = self.cazy_analyzer._clean_hmm_output(input_str)
        self.assertEqual(result, expected_output)

    def test_clean_hmm_output_invalid(self):
        input_str = ""
        expected_output = ""
        result = self.cazy_analyzer._clean_hmm_output(input_str)
        self.assertEqual(result, expected_output)

    def test_clean_dmnd_dbcansub_output_valid(self):
        input_str = "GH5_1+GH10_2"
        expected_output = "GH"
        result = self.cazy_analyzer._clean_dmnd_dbcansub_output(input_str)
        self.assertEqual(result, expected_output)

    def test_clean_dmnd_dbcansub_output_valid_diff_families(self):
        input_str = "GH5_1+GT10"
        expected_output = "GT;GH"
        result = self.cazy_analyzer._clean_dmnd_dbcansub_output(input_str)
        self.assertEqual(result, expected_output)

    def test_clean_dmnd_dbcansub_output_invalid(self):
        input_str = ""
        expected_output = ""
        result = self.cazy_analyzer._clean_dmnd_dbcansub_output(input_str)
        self.assertEqual(result, expected_output)

    # @patch("context.cazy.dict_to_dataframe")
    # @patch("pandas.DataFrame.to_excel")
    # @patch(
    #     "builtins.open",
    #     new_callable=unittest.mock.mock_open,
    #     read_data="header\nprotein\tec\tdmnd_r\tnum_tools\n",
    # )
    def test_parse_results_valid(
        self,
        # self, mock_open, mock_to_excel, mock_dict_to_dataframe
    ):
        # self.results_dir.mkdir(exist_ok=True)
        subdir = self.results_dir / "CAZY" / "CAZY_search" / "GCF_0000000"
        subdir.mkdir(exist_ok=True, parents=True)
        with open(subdir / "overview.txt", "w") as wf:
            wf.write("Gene ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n")
            wf.write("WP_016686262.1\t-\tGH23(468-593)\tGH23_e952\tGH23\t3\n")
            wf.write("WP_016686347.1\t2.4.1.1:158\tGT35(107-814)\tGT35_e0\tGT35\t3\n")

        self.cazy_analyzer.parse_results()
        excel_file = self.results_dir / "CAZY" / "CAZY_families.xlsx"
        self.assertTrue(excel_file.exists())

        res_df = read_excel(excel_file, index_col=0)
        exp_df = DataFrame.from_dict(
            {"GCF_0000000": {"AA": 0, "CBM": 0, "CE": 0, "GH": 1, "GT": 1, "PL": 0}},
            orient="index",
        )
        self.assertTrue(res_df.equals(exp_df))


#
#     @patch("context.cazy.dict_to_dataframe")
#     @patch("pandas.DataFrame.to_excel")
#     @patch(
#         "builtins.open",
#         new_callable=unittest.mock.mock_open,
#         read_data="header\nprotein\tec\tdmnd_r\tnum_tools\n",
#     )
#     def test_parse_results_invalid(
#         self, mock_open, mock_to_excel, mock_dict_to_dataframe
#     ):
#         with self.assertRaises(Exception):
#             self.cazy_analyzer.results_dir = Path("/invalid/path/to/results")
#             self.cazy_analyzer.parse_results()
