import unittest
import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

from context import orthologues


class TestOrthologuesIdentifier(unittest.TestCase):
    def setUp(self):
        self.fasta_files_list = [Path("file1.fasta"), Path("file2.fasta")]
        self.out_dir = Path("output")
        self.ref = "reference"
        self.ref_list = None
        self.input_type = "protein"
        self.available_cores = 4
        self.blast_cores = 2
        self.concurrent = True
        self.evalue = 0.001
        self.dmnd_sensitivity = "sensitive"
        self.no_filter = False

        self.ortho_id = orthologues.Orthologues_identifier(
            fasta_files_list=self.fasta_files_list,
            out_dir=self.out_dir,
            ref=self.ref,
            ref_list=self.ref_list,
            input_type=self.input_type,
            available_cores=self.available_cores,
            blast_cores=self.blast_cores,
            concurrent=self.concurrent,
            evalue=self.evalue,
            dmnd_sensitivity=self.dmnd_sensitivity,
            no_filter=self.no_filter,
        )

    @patch("pathlib.Path.mkdir")
    def test_create_directories(self, mock_mkdir):
        self.ortho_id._create_directories(self.ref)
        self.assertEqual(mock_mkdir.call_count, 3)
        mock_mkdir.assert_any_call(exist_ok=True, parents=True)

    def test_get_blast_binaries(self):
        self.ortho_id._get_blast_binaries()
        self.assertEqual(self.ortho_id.blast_bin, "diamond")

        self.ortho_id.input_type = "CDS"
        self.ortho_id._get_blast_binaries()
        self.assertEqual(self.ortho_id.blast_bin, "blastn")
        self.assertEqual(self.ortho_id.blast_db_bin, "makeblastdb")

    @patch("builtins.open", new_callable=mock_open, read_data="ref1\nref2\n")
    @patch("context.orthologues.Orthologues_identifier._create_directories")
    @patch("context.orthologues.Orthologues_identifier._get_blast_binaries")
    def test_setup(self, mock_get_blast_binaries, mock_create_directories, mock_open):
        self.ortho_id.ref_list = "dummy_ref_list.txt"
        self.ortho_id.setup()
        mock_create_directories.assert_called()
        mock_get_blast_binaries.assert_called()
        self.assertEqual(mock_create_directories.call_count, 3)

    @patch("context.orthologues.execute_command")
    @patch("tqdm.tqdm", lambda x, **kwargs: x)  # To bypass the progress bar
    def test_create_blast_db(self, mock_execute_command):
        self.ortho_id.fasta_files = [Path("file1.fasta")]
        self.ortho_id.blast_bin = "diamond"
        self.ortho_id.create_blast_db()
        self.assertEqual(mock_execute_command.call_count, 1)
        dbfile = self.ortho_id.blast_db_dir / "file1"
        mock_execute_command.assert_called_once_with(
            f"diamond makedb --in file1.fasta --quiet --db {dbfile} --threads 2"
        )

    @pytest.mark.skip(reason="Need to also add integration test for this function")
    @patch("context.orthologues.multiprocess_dispatch")
    def test_perform_reciprocal_blast(self, mock_multiprocess_dispatch):
        self.ortho_id.blast_db_dir.mkdir(parents=True, exist_ok=True)
        self.ortho_id.ref = "file1"
        self.ortho_id.blast_bin = "diamond"
        ref_fasta = self.ortho_id.perform_reciprocal_blast(self.ortho_id.ref)
        self.assertIsNotNone(ref_fasta)
        self.assertIn("file1", ref_fasta.name)
        self.assertEqual(mock_multiprocess_dispatch.call_count, 1)

    @patch("context.orthologues.multiprocess_dispatch")
    def test_perform_reciprocal_blast_wrong_ref(self, mock_multiprocess_dispatch):
        self.ortho_id.blast_db_dir.mkdir(parents=True, exist_ok=True)
        self.ortho_id.ref = "FAKE"
        self.ortho_id.blast_bin = "diamond"
        with self.assertRaises(FileNotFoundError):
            _ = self.ortho_id.perform_reciprocal_blast(self.ortho_id.ref)

    def test_get_best_subject(self):
        # Create a sample DataFrame
        df = pd.DataFrame(
            {
                "RefSeq": ["r1", "r1", "r2", "r2", "r3"],
                "Query": ["q1", "q2", "q3", "q4", "q5"],
                "Bitscore": [20, 10, 40, 38, 50],
            }
        )
        expected_df = pd.DataFrame(
            {
                "RefSeq": ["r1", "r2", "r3"],
                "Query": ["q1", "q3", "q5"],
                "Bitscore": [20, 40, 50],
            }
        )

        # Group by 'Query' and get the index of the row with maximum Bitscore
        result = self.ortho_id._get_best_subject(df)
        self.assertEqual(result.equals(expected_df), True)

    def test_create_reciprocal_matrix(self):
        # Create sample DataFrames for reference vs query and query vs reference
        ref_vs_query_df = pd.DataFrame(
            {
                "RefSeq": ["r1", "r2", "r3"],
                "QuerySeq": ["q1", "q2", "q3"],
                "Pident": [100, 60, 70],
                "Evalue": [0.01, 0.02, 0.03],
                "Bitscore": [1000, 200, 300],
            }
        )

        query_vs_ref_df = pd.DataFrame(
            {
                "QuerySeq": ["q1", "q2", "q4"],
                "RefSeq": ["r1", "r2", "r3"],
                "Pident": [60, 50, 70],
                "Evalue": [0.02, 0.01, 0.03],
                "Bitscore": [200, 100, 300],
            }
        )

        result = self.ortho_id._create_reciprocal_matrix(
            ref_vs_query_df, query_vs_ref_df
        )

        expected_df = pd.DataFrame(
            {
                "RefSeq": ["r1", "r2"],
                "QuerySeq": ["q1", "q2"],
                "Pident": [100, 60],
            }
        )
        self.assertEqual(result.shape, (2, 3))
        self.assertEqual(result.equals(expected_df), True)

    def test_orthologue_filter(self):
        # Create a sample DataFrame
        df = pd.DataFrame(
            {
                "RefSeq": ["k"] * 100,
                "QuerySeq": ["q"] * 100,
                "Pident": [100] * 99 + [5],
                "Evalue": [0.05] * 100,
                "Bitscore": [100] * 100,
            }
        )

        result = self.ortho_id._orthologue_filter(df)
        self.assertEqual(result.shape, (99, 5))
        self.assertEqual(result.Pident.min(), 100)

    @pytest.mark.skip(reason="Not implemented")
    @patch("pandas.read_csv")
    @patch("pandas.DataFrame.to_csv")
    @patch("tqdm.tqdm", lambda x, **kwargs: x)  # To bypass the progress bar
    def test_parse_blast_results(self, mock_to_csv, mock_read_csv):
        mock_read_csv.return_value = pd.DataFrame(
            {
                "RefSeq": ["seq1", "seq2"],
                "QuerySeq": ["query1", "query2"],
                "Pident": [99.9, 98.5],
                "Evalue": [0.0, 0.001],
                "Bitscore": [500, 450],
            }
        )

        self.ortho_id.parse_blast_results(self.ref)
        self.assertEqual(mock_read_csv.call_count, 2)
        self.assertEqual(mock_to_csv.call_count, 1)

    @pytest.mark.skip(reason="Not implemented")
    @patch("Bio.SeqIO.parse")
    @patch("pandas.DataFrame.from_dict")
    @patch("pandas.DataFrame.to_csv")
    def test_create_orthology_matrix(
        self, mock_to_csv, mock_from_dict, mock_seqio_parse
    ):
        mock_seqio_parse.return_value = [MagicMock(id="gene1"), MagicMock(id="gene2")]
        mock_from_dict.return_value = pd.DataFrame(
            {
                self.ref: ["gene1", "gene2"],
                "genome1": ["seq1", "seq2"],
                "genome2": ["X", "seq4"],
            }
        )

        self.ortho_id.create_orthology_matrix(self.ref, Path("reference.fasta"))
        self.assertEqual(mock_seqio_parse.call_count, 1)
        self.assertEqual(mock_from_dict.call_count, 1)
        self.assertEqual(mock_to_csv.call_count, 1)

    @patch("context.orthologues.recursive_unlink")
    @patch("context.orthologues.Orthologues_identifier.create_orthology_matrix")
    @patch("context.orthologues.Orthologues_identifier.parse_blast_results")
    @patch("context.orthologues.Orthologues_identifier.perform_reciprocal_blast")
    @patch("context.orthologues.Orthologues_identifier.create_blast_db")
    @patch("context.orthologues.Orthologues_identifier.setup")
    def test_calculate_orthologues(
        self,
        mock_setup,
        mock_create_blast_db,
        mock_perform_reciprocal_blast,
        mock_parse_blast_results,
        mock_create_orthology_matrix,
        mock_recursive_unlink,
    ):
        self.ortho_id.calculate_orthologues()

        mock_setup.assert_called_once()
        mock_create_blast_db.assert_called_once()
        self.assertEqual(mock_perform_reciprocal_blast.call_count, 1)
        self.assertEqual(mock_parse_blast_results.call_count, 1)
        self.assertEqual(mock_create_orthology_matrix.call_count, 1)
        mock_recursive_unlink.assert_called_once()
