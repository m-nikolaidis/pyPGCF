import unittest
from pandas import DataFrame, read_excel
from context import cazy, utils
from unittest.mock import MagicMock, patch, mock_open
from pathlib import Path


class TestCAZYAnalyzer(unittest.TestCase):
    def setUp(self):
        self.cores = 3
        self.evalue = 1e-5
        self.database_dir = Path("/tmp/database")
        self.results_dir = Path("/tmp/results")
        self.fasta_files = [Path("/tmp/f1.fa"), Path("/tmp/f2.fa")]
        self.protein = True
        self.verbose = False
        self.available_cores = 6
        self.concurrent = True
        self.cazy_analyzer = cazy.CAZY_analyzer(
            cores=self.cores,
            evalue=self.evalue,
            database_dir=self.database_dir,
            out_dir=self.results_dir,
            fasta_files_list=self.fasta_files,
            input_type=self.protein,
            verbose=self.verbose,
            concurrent=self.concurrent,
            available_cores=self.available_cores,
        )

    # @patch("context.cazy.calc_avail_dispatchers", return_value=2)
    # @patch("context.cazy.multiprocess_dispatch")
    # @patch("context.cazy.create_temporary_file", return_value=Path("/fake/temp/file"))
    # @patch("context.cazy.translate_fasta_records", return_value=["record1", "record2"])
    # @patch("context.cazy.seqrecords_to_fasta")
    # def test_execute_cazy_search_valid(
    #     self,
    #     mock_seqrecords_to_fasta,
    #     mock_translate_fasta_records,
    #     mock_create_temporary_file,
    #     mock_multiprocess_dispatch,
    #     mock_calc_avail_dispatchers,
    # ):
    #     # self.fasta_files.mkdir(exist_ok=True)
    #     with open(self.fasta_files / "test.fa", "w") as wf:
    #         wf.write(
    #             ">AAC07180\nMVALRNNWADISRKILEKRRFETKDLVKRLKVITGHDIHIQNYPVETPRVAFNPSIHVFENRLRIYARVV MGYYTYTSAIAEFDIDLEELYNPERKTYEANLTVLPNIKYDLWGVEDPRVYEIDGKLFMTYTGRTVNYFR TDIRTERTLPVTARYENGQWKKIAVFRMPEDIRSFVVSDKNAFLVKTDKLMLYHRLHMLNEKFYLAVCNV PEEVLYTNEFKEIEIGENITIMEEAPFETKIGWATPPVKVGEENLVLIHGVDKELTAYRVFAVLMNKEGY FTAVTPFYILEPKKIYEVYGDRPFVVFPCGIQRLENKLLISYGGADTVVVIGEIDLEELMNILYENRID"
    #         )
    #     self.cazy_analyzer.execute_cazy_search()
    #     mock_multiprocess_dispatch.assert_called_once()

    # @patch("context.cazy.calc_avail_dispatchers", return_value=2)
    # @patch("context.cazy.multiprocess_dispatch")
    # @patch("context.cazy.create_temporary_file", return_value=Path("/fake/temp/file"))
    # @patch("context.cazy.translate_fasta_records", return_value=["record1", "record2"])
    # @patch("context.cazy.seqrecords_to_fasta")
    # def test_execute_cazy_search_invalid(
    #     self,
    #     mock_seqrecords_to_fasta,
    #     mock_translate_fasta_records,
    #     mock_create_temporary_file,
    #     mock_multiprocess_dispatch,
    #     mock_calc_avail_dispatchers,
    # ):
    #     self.cazy_analyzer.fasta_files = []
    #     self.assertRaises(FileNotFoundError, self.cazy_analyzer.execute_cazy_search)

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
        # expected_output = "GH;GT"
        result = self.cazy_analyzer._clean_dmnd_dbcansub_output(input_str)
        self.assertEqual(result, expected_output)

    def test_clean_dmnd_dbcansub_output_invalid(self):
        input_str = ""
        expected_output = ""
        result = self.cazy_analyzer._clean_dmnd_dbcansub_output(input_str)
        self.assertEqual(result, expected_output)

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

    # TODO: Use the following function to mock the behaviour of the one above
    # @patch("pandas.DataFrame.to_excel", MagicMock())
    # @patch(
    #     "builtins.open",
    #     new_callable=mock_open,
    #     # read_data="header\nprot\tec\t_\t_\tdmnd_r\tnum_tools\nprotein\t-\t-\t-\tGH1_1+GH3_2\t3\n",
    # )
    # @patch(
    #     "pathlib.Path.glob",
    #     MagicMock(return_value=[Path("/dummy/out_dir/CAZY_search/genome1")]),
    # )
    # @patch("pathlib.Path.exists", MagicMock(return_value=True))
    # def test_parse_results(self, mock_open_file):
    #     self.cazy_analyzer.parse_results()
    #
    #     mock_open_file.assert_called_once_with(
    #         Path("/dummy/out_dir/CAZY/CAZY_search/genome1/overview.txt"), "r"
    #     )
    #     utils.dict_to_dataframe.assert_called()
    #     utils.dict_to_dataframe.return_value.to_excel.assert_called_with(
    #         Path("/dummy/out_dir/CAZY/CAZY_families.xlsx")
    #     )


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
