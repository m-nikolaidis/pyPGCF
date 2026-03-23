import unittest
import pytest
from pathlib import Path
from unittest.mock import patch, call
import pandas as pd
from context import amr


class TestAMRAnalyzer(unittest.TestCase):
    @patch("context.amr.calc_avail_dispatchers")
    def test_init(self, mock_calc_avail_dispatchers):
        mock_calc_avail_dispatchers.return_value = 4

        analyzer = amr.AMR_analyzer(
            fasta_files_list=[Path("/test/fasta1"), Path("/test/fasta2")],
            database_dir=Path("/test/database"),
            out_dir=Path("/test/out"),
            input_type="prot",
            available_cores=8,
            blast_cores=2,
            concurrent=True,
            debug=True,
        )

        self.assertTrue(analyzer.protein)
        self.assertEqual(analyzer.out_dir, Path("/test/out/AMR"))
        self.assertEqual(len(analyzer.fasta_files), 2)
        self.assertTrue(analyzer.debug)
        self.assertEqual(analyzer.concurrent_jobs, 4)
        self.assertEqual(analyzer.database_dir, Path("/test/database/AMR"))

        # Test with input_type "nucl" and database_dir as None
        analyzer = amr.AMR_analyzer(
            fasta_files_list=[Path("/test/fasta1")],
            database_dir=None,
            out_dir=Path("/test/out"),
            input_type="nucl",
            available_cores=8,
            blast_cores=2,
            concurrent=False,
            debug=False,
        )
        self.assertFalse(analyzer.protein)
        self.assertEqual(analyzer.database_dir, None)
        self.assertEqual(analyzer.concurrent_jobs, 1)  # Non-concurrent case

    @patch("context.amr.Path.mkdir")
    @patch("context.amr.multiprocess_dispatch")
    def test_search_amr(self, mock_multiprocess_dispatch, mock_mkdir):
        analyzer = amr.AMR_analyzer(
            fasta_files_list=[Path("/test/fasta1"), Path("/test/fasta2")],
            database_dir=Path("/test/database"),
            out_dir=Path("/test/out"),
            input_type="prot",
            available_cores=8,
            blast_cores=2,
            concurrent=True,
            debug=False,
        )

        analyzer.search_amr()

        mock_mkdir.assert_called_once_with(exist_ok=True)
        expected_cmds = [
            "amrfinder --plus -i -1 -o /test/out/AMR/amrfinder/fasta1.txt --threads 2 -d /test/database/AMR -p /test/fasta1 --quiet",
            "amrfinder --plus -i -1 -o /test/out/AMR/amrfinder/fasta2.txt --threads 2 -d /test/database/AMR -p /test/fasta2 --quiet",
        ]
        mock_multiprocess_dispatch.assert_called_once_with(
            "system",
            expected_cmds,
            3,
            show_progress=True,
            description="Scanning for AMR genes",
        )

    @pytest.mark.skip(reason="Not implemented")
    @patch("context.amr.read_csv")
    @patch("context.amr.concat")
    @patch("context.amr.Path.glob")
    @patch("context.amr.DataFrame.to_excel")
    def test_parse_results(self, mock_to_excel, mock_glob, mock_concat, mock_read_csv):
        mock_glob.return_value = [
            Path("/test/out/AMR/res1.txt"),
            Path("/test/out/AMR/res2.txt"),
        ]

        df1 = pd.DataFrame(
            {"Element type": ["AMR", "non-AMR"], "Data": [1, 2]}
        ).set_index("Element type")
        df2 = pd.DataFrame({"Element type": ["AMR"], "Data": [3]}).set_index(
            "Element type"
        )
        mock_read_csv.side_effect = [df1, df2]

        analyzer = amr.AMR_analyzer(
            fasta_files_list=[Path("/test/fasta1"), Path("/test/fasta2")],
            database_dir=Path("/test/database"),
            out_dir=Path("/test/out"),
            input_type="prot",
            available_cores=8,
            blast_cores=2,
            concurrent=True,
            debug=False,
        )

        analyzer.parse_results()

        expected_calls = [
            call(Path("/test/out/AMR/res1.txt"), sep="\t", index_col=0),
            call(Path("/test/out/AMR/res2.txt"), sep="\t", index_col=0),
        ]
        mock_read_csv.assert_has_calls(expected_calls)

        mock_concat.assert_called_once()

        mock_to_excel.assert_called_once_with(
            Path("/test/out/AMR/AMRfinder_results.xlsx")
        )
