from pathlib import Path
import unittest
from context import species_demarcation

from unittest.mock import patch, mock_open, MagicMock
import pandas as pd
import pandas.testing as pd_testing
import pytest


class TestSpeciesDemarcator(unittest.TestCase):
    def setUp(self):
        self.in_dir = Path("/path/to/input")
        self.out_dir = Path("/path/to/output")
        self.fastani_cores = 4
        self.kmer = 16
        self.fraglen = 3000
        self.minfrac = 0.2
        self.inflation = 1.5
        self.mcl_cores = 4
        self.debug = True

        self.demarcator = species_demarcation.SpeciesDemarcator(
            in_dir=self.in_dir,
            out_dir=self.out_dir,
            fastani_cores=self.fastani_cores,
            kmer=self.kmer,
            fraglen=self.fraglen,
            minfrac=self.minfrac,
            inflation=self.inflation,
            mcl_cores=self.mcl_cores,
            debug=self.debug,
        )

    @patch("pathlib.Path.mkdir")
    def test_create_directories(self, mock_mkdir):
        self.demarcator.create_directories()
        mock_mkdir.assert_called_once_with(exist_ok=True, parents=True)

    @patch("builtins.open", new_callable=mock_open)
    def test_create_input_for_fastani(self, mock_open):
        files = [Path(f"/path/to/genome_{i}.fa") for i in range(5)]
        tmp_file = Path("/path/to/tmp_file.txt")
        self.demarcator.create_input_for_fastani(files, tmp_file)

        mock_open.assert_called_once_with(str(tmp_file), "w")
        mock_open().write.assert_has_calls(
            [unittest.mock.call(str(f) + "\n") for f in files]
        )

    @pytest.mark.skip(reason="Needs existing file")
    @patch("context.utils.execute_command")
    def test_perform_fastani(self, mock_execute_command):
        org_list = Path("/path/to/org_list.txt")
        fout = Path("/path/to/output.txt")
        mock_execute_command.return_value = 0

        self.demarcator.perform_fastani(org_list, fout)
        cmd = (
            f"fastANI --ql {org_list} --rl {org_list} -t {self.fastani_cores} -k {self.kmer} "
            f"--fragLen {self.fraglen} --minFraction {self.minfrac} -o {fout}"
        )
        mock_execute_command.assert_called_once_with(cmd)

    @patch("context.utils.execute_command")
    def test_perform_fastani_failure(self, mock_execute_command):
        org_list = Path("/path/to/org_list.txt")
        fout = Path("/path/to/output.txt")
        mock_execute_command.return_value = 1

        with self.assertRaises(RuntimeError):
            self.demarcator.perform_fastani(org_list, fout)

    @patch("pandas.read_csv")
    @patch("pandas.DataFrame.to_csv")
    def test_prepare_input_for_mcl(self, mock_to_csv, mock_read_csv):
        mock_df = pd.DataFrame(
            {
                "query": ["q1", "q2"],
                "target": ["t1", "t1"],
                "ANI": [99, 80],
                "query_length": [3000, 3000],
                "target_length": [3000, 1000],
            }
        )
        mock_read_csv.return_value = mock_df

        input_file = Path("/path/to/fastani_output.tsv")
        expected_fout = input_file.parent / "fastANI_for_mcl.txt"

        actual_fout = self.demarcator.prepare_input_for_mcl(input_file)

        # filtered df
        expected_df = pd.DataFrame(
            {
                "query": ["q1"],
                "target": ["t1"],
                "ANI": [99],
                "query_length": [3000],
                "target_length": [3000],
            }
        )
        actual_df = mock_df[mock_df["ANI"] >= 95]

        self.assertEqual(expected_fout, actual_fout)
        self.assertTrue(mock_to_csv.called)
        self.assertEqual(expected_df.equals(actual_df), True)

    @pytest.mark.skip(reason="Needs existing files to work")
    @patch("pathlib.Path.unlink")
    @patch("pathlib.Path.rename")
    def test_clean_mcl(self, mock_rename, mock_unlink):
        outdir = Path("/path/to/output")
        expected_to_remove = [
            outdir / "fastANI_for_mcl.txt",
            outdir / "fastANI_mcx_mtrx.txt",
            outdir / "fastANI_annot.tab",
            outdir / "fastANI_mcl_out.txt",
            outdir / "FastANI_input.txt",
        ]
        expected_rename = (
            outdir / "fastANI_mcx_dump.txt",
            outdir / "fastANI_clusters.tsv",
        )

        self.demarcator.clean_mcl(outdir)

        mock_unlink.assert_has_calls(
            [unittest.mock.call(f, missing_ok=True) for f in expected_to_remove]
        )
        mock_rename.assert_called_once_with(*expected_rename)

    @pytest.mark.skip(reason="Needs existing files to work")
    @patch("utils.execute_command")
    def test_run_mcl(self, mock_execute_command):
        mock_execute_command.return_value = 0

        fastani_for_mcl = Path("/path/to/fastANI_for_mcl.txt")
        actual_fout = self.demarcator.run_mcl(fastani_for_mcl)
        expected_fout = fastani_for_mcl.parent / "fastANI_clusters.tsv"

        self.assertEqual(actual_fout, expected_fout)
        self.assertTrue(mock_execute_command.called)

    @pytest.mark.skip(reason="Needs existing files to work")
    @patch("builtins.open", new_callable=mock_open, read_data="genome1\t0\n")
    @patch("pandas.DataFrame.to_excel")
    def test_parse_mcx_output(self, mock_to_excel, mock_open):
        fastani_mcl_output = Path("/path/to/fastANI_clusters.tsv")

        self.demarcator.parse_mcx_output(fastani_mcl_output)

        self.assertTrue(mock_open.called)
        self.assertTrue(mock_to_excel.called)

    @patch("context.species_demarcation.SpeciesDemarcator.create_directories")
    @patch("context.species_demarcation.SpeciesDemarcator.create_input_for_fastani")
    @patch("context.species_demarcation.SpeciesDemarcator.perform_fastani")
    @patch("context.species_demarcation.SpeciesDemarcator.prepare_input_for_mcl")
    @patch("context.species_demarcation.SpeciesDemarcator.run_mcl")
    @patch("context.species_demarcation.SpeciesDemarcator.parse_mcx_output")
    @patch("context.species_demarcation.check_if_file_exists", return_value=True)
    def test_assign_species(
        self,
        mock_check_if_file_exists,
        mock_parse_mcx_output,
        mock_run_mcl,
        mock_prepare_input_for_mcl,
        mock_perform_fastani,
        mock_create_input_for_fastani,
        mock_create_directories,
    ):
        self.demarcator.assign_species()
        mock_create_directories.assert_called_once()
        mock_create_input_for_fastani.assert_called_once()
        mock_perform_fastani.assert_called_once()
        mock_prepare_input_for_mcl.assert_called_once()
        mock_run_mcl.assert_called_once()
        mock_parse_mcx_output.assert_called_once()


# @dataclass
# class Data:
#      cores = 4
#      kmer = 16
#      fraglen = 3000
#      minfraction = 0.2
#      mcl_inflation = 2
#      data_dir = Path(__file__).parent / "../data/genomes/"
#      outdir = Path(__file__).parent / "../data/species_demarcation/"
#      outdir.mkdir(exist_ok=True)
#
#
# class TestModule(unittest.TestCase):
#
#     def test_assign_species(self):
#         data = Data()
#         demarcator = species_demarcation.SpeciesDemarcator(
#             in_dir=data.data_dir,
#             out_dir=data.outdir,
#             fastani_cores=data.cores,
#             kmer=data.kmer,
#             fraglen=data.fraglen,
#             minfrac=data.minfraction,
#             inflation=data.mcl_inflation,
#             mcl_cores=data.cores,
#             debug=True,
#         )
#         demarcator.assign_species()
#
#     def test_create_input_for_fastani(self):
#         data = Data()
#         demarcator = species_demarcation.SpeciesDemarcator(
#             in_dir=data.data_dir,
#             out_dir=data.outdir,
#             fastani_cores=data.cores,
#             kmer=data.kmer,
#             fraglen=data.fraglen,
#             minfrac=data.minfraction,
#             inflation=data.mcl_inflation,
#             mcl_cores=data.cores,
#         )
#         files_for_fastani = list(data.data_dir.glob("*.fna"))
#         org_list = data.outdir / "org_list.txt"
#         demarcator.create_input_for_fastani(files_for_fastani, org_list)
#         self.assertIs(org_list.exists(), True)
#
#     def test_perform_fastani(self):
#         data = Data()
#         demarcator = species_demarcation.SpeciesDemarcator(
#             in_dir=data.data_dir,
#             out_dir=data.outdir,
#             fastani_cores=data.cores,
#             kmer=data.kmer,
#             fraglen=data.fraglen,
#             minfrac=data.minfraction,
#             inflation=data.mcl_inflation,
#             mcl_cores=data.cores,
#         )
#         org_list = data.outdir / "org_list.txt"
#         files_for_fastani = list(data.data_dir.glob("*.fna"))
#         fastani_fout = data.outdir / "FastANI.tsv"
#         demarcator.create_input_for_fastani(files_for_fastani, org_list)
#         demarcator.perform_fastani(org_list, fastani_fout)
#         self.assertIs(fastani_fout.exists(), True)
#
#     def test_prepare_input_for_mcl(self):
#         data = Data()
#         demarcator = species_demarcation.SpeciesDemarcator(
#             in_dir=data.data_dir,
#             out_dir=data.outdir,
#             fastani_cores=data.cores,
#             kmer=data.kmer,
#             fraglen=data.fraglen,
#             minfrac=data.minfraction,
#             inflation=data.mcl_inflation,
#             mcl_cores=data.cores,
#         )
#         fastani_fout = data.outdir / "FastANI.tsv"
#         mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
#         self.assertIs(mcl_input.exists(), True)
#
#     def test_run_mcl(self):
#         data = Data()
#         demarcator = species_demarcation.SpeciesDemarcator(
#             in_dir=data.data_dir,
#             out_dir=data.outdir,
#             fastani_cores=data.cores,
#             kmer=data.kmer,
#             fraglen=data.fraglen,
#             minfrac=data.minfraction,
#             inflation=data.mcl_inflation,
#             mcl_cores=data.cores,
#             debug=True
#         )
#         fastani_fout = data.outdir / "FastANI.tsv"
#         mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
#         mcl_output = demarcator.run_mcl(mcl_input)
#         self.assertIs(mcl_output.exists(), True)
#
#     def test_parse_mcx_output(self):
#         data = Data()
#         demarcator = species_demarcation.SpeciesDemarcator(
#             in_dir=data.data_dir,
#             out_dir=data.outdir,
#             fastani_cores=data.cores,
#             kmer=data.kmer,
#             fraglen=data.fraglen,
#             minfrac=data.minfraction,
#             inflation=data.mcl_inflation,
#             mcl_cores=data.cores,
#         )
#         fastani_fout = data.outdir / "FastANI.tsv"
#         mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
#         mcl_output = demarcator.run_mcl(mcl_input)
#         demarcator.parse_mcx_output(mcl_output)
#         fastani_from_mcl = data.outdir / "FastANI_species_clusters.xlsx"
#         self.assertIs(fastani_from_mcl.exists(), True)
#
#     def test_read_fastani_output(self):
#         ...
#         # data = Data()
#         # demarcator = species_demarcation.SpeciesDemarcator(
#         #     in_dir=data.data_dir,
#         #     out_dir=data.outdir,
#         #     fastani_cores=data.cores,
#         #     kmer=data.kmer,
#         #     fraglen=data.fraglen,
#         #     minfrac=data.minfraction,
#         #     inflation=data.mcl_inflation,
#         #     mcl_cores=data.cores,
#         # )
#         # fastani_fout = data.outdir / "FastANI.tsv"
#         # mcl_input = demarcator.prepare_input_for_mcl(fastani_fout)
#         # mcl_output = demarcator.run_mcl(mcl_input)
#         # demarcator.parse_mcx_output(mcl_output)
#         # fastani_from_mcl = data.outdir / "FastANI_species_clusters.xlsx"
#         # df = pd.read_excel(fastani_from_mcl, index_col=0)
#         # test_df = pd.DataFrame.from_dict(
#         #     {
#         #         "GCA_000769555": {"FastANI_species": "C0"},
#         #         "GCA_002220285": {"FastANI_species": "C1"},
#         #         "GCA_000009045": {"FastANI_species": "C2"},
#         #     },
#         #     orient="index",
#         # )
#         # print(test_df)
#         # print(df)
#         # self.assertTrue(df.equals(test_df))
#         # return None
