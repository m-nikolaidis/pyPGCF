import unittest
from unittest.mock import patch, mock_open
from pathlib import Path
import pandas as pd
from context import virulence
import tempfile
import pytest


class TestVFAnalyzer(unittest.TestCase):
    def setUp(self):
        self.fasta_files_list = [Path("file1.fasta"), Path("file2.fasta")]
        self.database_dir = Path("/db")
        self.out_dir = Path("/output")
        self.input_type = "prot"
        self.available_cores = 4
        self.blast_cores = 2
        self.evalue = 0.001
        self.dmnd_sensitivity = "sensitive"
        self.concurrent = True
        self.debug = False

    @patch("context.virulence.calc_avail_dispatchers", return_value=2)
    def test_init(self, mock_calc_avail_dispatchers):
        analyzer = virulence.VF_analyzer(
            fasta_files_list=self.fasta_files_list,
            database_dir=self.database_dir,
            out_dir=self.out_dir,
            input_type=self.input_type,
            available_cores=self.available_cores,
            blast_cores=self.blast_cores,
            evalue=self.evalue,
            dmnd_sensitivity=self.dmnd_sensitivity,
            concurrent=self.concurrent,
            debug=self.debug,
        )

        self.assertEqual(analyzer.cores, self.blast_cores)
        self.assertTrue(analyzer.protein)
        self.assertEqual(analyzer.database_file, self.database_dir / "VFDB_setA_pro")
        self.assertEqual(analyzer.concurrent_jobs, 2)

    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="VF1\tCategory1\tOrigin1\tDescription1\n",
    )
    def test_read_vf_desc_filt(self, mock_open):
        analyzer = virulence.VF_analyzer(
            fasta_files_list=self.fasta_files_list,
            database_dir=self.database_dir,
            out_dir=self.out_dir,
            input_type=self.input_type,
            available_cores=self.available_cores,
            blast_cores=self.blast_cores,
            evalue=self.evalue,
            dmnd_sensitivity=self.dmnd_sensitivity,
            concurrent=self.concurrent,
            debug=self.debug,
        )
        vf_desc = analyzer._read_vf_desc_filt()
        mock_open.assert_called_once_with(self.database_dir / "vfdb_desc.tsv", "r")
        self.assertEqual(
            vf_desc["VF1"],
            {
                "VF_Category": "Category1",
                "VF_Origin": "Origin1",
                "VF_Desc": "Description1",
            },
        )

    # @pytest.mark.skip(reason="Not implemented")
    @patch("context.virulence.create_diamond_blastp_cmd", return_value="blastp_cmd")
    @patch("context.virulence.multiprocess_dispatch", return_value=None)
    @patch("pathlib.Path.mkdir")
    def test_execute_homology_search(
        self, mock_mkdir, mock_multiprocess_dispatch, mock_create_diamond_blastp_cmd
    ):
        analyzer = virulence.VF_analyzer(
            fasta_files_list=self.fasta_files_list,
            database_dir=self.database_dir,
            out_dir=self.out_dir,
            input_type=self.input_type,
            available_cores=self.available_cores,
            blast_cores=self.blast_cores,
            evalue=self.evalue,
            dmnd_sensitivity=self.dmnd_sensitivity,
            concurrent=self.concurrent,
            debug=self.debug,
        )
        analyzer.execute_homology_search()
        mock_mkdir.assert_called_once_with(exist_ok=True)
        mock_create_diamond_blastp_cmd.assert_called()
        mock_multiprocess_dispatch.assert_called_once()

    @patch("pandas.read_csv")
    @patch("pandas.concat")
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="VF1\tCategory1\tOrigin1\tDescription1\n",
    )
    @patch("context.virulence.keep_best_homology_hit", return_value=pd.DataFrame())
    @patch(
        "context.virulence.virulence.VF_analyzer._read_vf_desc_filt",
        return_value={
            "VF1": {"VF_Category": "Cat1", "VF_Origin": "Org1", "VF_Desc": "Desc1"}
        },
    )
    @pytest.mark.skip(reason="Not implemented")
    @patch("pandas.DataFrame.to_excel")
    @patch("context.virulence.Path.glob")
    def test_parse_results(
        self,
        mock_glob,
        mock_to_excel,
        # mock_read_vf_desc_filt,
        # mock_keep_best_homology_hit,
        # mock_open,
        # mock_concat,
        # mock_read_csv,
    ):
        temp_dir = tempfile.TemporaryDirectory()
        analyzer = virulence.VF_analyzer(
            fasta_files_list=self.fasta_files_list,
            database_dir=self.database_dir,
            out_dir=Path(temp_dir.name),
            input_type=self.input_type,
            available_cores=self.available_cores,
            blast_cores=self.blast_cores,
            evalue=self.evalue,
            dmnd_sensitivity=self.dmnd_sensitivity,
            concurrent=self.concurrent,
            debug=self.debug,
        )

        mock_glob.return_value = [Path("file1.txt"), Path("file2.txt")]
        mock_df = pd.DataFrame(
            {
                "qseqid": ["VF1"],
                "sseqid": ["SEQ1"],
                "qcovhsp": [60],
                "pident": [70],
                "evalue": [0.0001],
                "bitscore": [50],
            }
        )
        # mock_read_csv.return_value = mock_df
        # mock_concat.return_value = mock_df
        # analyzer.parse_results()
        #
        # mock_read_csv.assert_called()
        # mock_concat.assert_called_once()
        # mock_to_excel.assert_called_once_with(
        #     Path(temp_dir.name) / "VF_results_50pident50qcov.xlsx"
        # )
