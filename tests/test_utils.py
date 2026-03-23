from pathlib import Path
import unittest
from unittest.mock import patch, call, MagicMock
import tempfile
import pandas as pd

from pandas.testing import assert_frame_equal
from pandas import DataFrame

from context import utils


def _helper(x: int):
    return x + 5


class TestUtilityModule(unittest.TestCase):
    @patch("pathlib.Path.is_file", return_value=True)
    @patch("pathlib.Path.unlink")
    def test_recursive_unlink_file(self, mock_unlink, mock_is_file):
        input_path = Path("/fake/file")
        utils.recursive_unlink(input_path)
        mock_unlink.assert_called_once_with()

    @patch("pathlib.Path.is_file", return_value=False)
    @patch("pathlib.Path.glob")
    @patch("pathlib.Path.rmdir")
    @patch("pathlib.Path.unlink")
    def test_recursive_unlink_directory(
        self, mock_unlink, mock_rmdir, mock_glob, mock_is_file
    ):
        # Mock the directory contents
        mock_dir_content = [Path("/fake/dir/file1"), Path("/fake/dir/file2")]
        mock_glob.return_value = mock_dir_content

        input_path = Path("/fake/dir")
        utils.recursive_unlink(input_path)

        mock_glob.assert_called_once_with("*")
        mock_unlink.assert_has_calls([call(), call()])
        mock_rmdir.assert_called_once_with()

    def test_recursive_unlink_complex_dir(self):
        tmp_dir = Path(tempfile.mkdtemp())
        self.assertIs(tmp_dir.exists(), True)

        # Create child directories and populate them with 1 file
        child_dirs = [tmp_dir / str(i) for i in range(3)]
        for c in child_dirs:
            c.mkdir()
            tmp_file = c / "file.txt"
            tmp_file.touch()
            self.assertIs(c.exists(), True)
            self.assertIs(tmp_file.exists(), True)

        utils.recursive_unlink(tmp_dir)
        self.assertIs(tmp_dir.exists(), False)

    def test_calc_avail_dispatchers(self):
        # basic
        self.assertEqual(utils.calc_avail_dispatchers(8, 2), 3)

        # avoid throttle
        self.assertEqual(utils.calc_avail_dispatchers(8, 2, avoid_throttle=True), 3)
        self.assertEqual(utils.calc_avail_dispatchers(8, 2, avoid_throttle=False), 4)

        # usemax
        self.assertEqual(
            utils.calc_avail_dispatchers(8, 2, avoid_throttle=True, usemax=True), 4
        )

        # more jobs than available
        with self.assertRaises(ValueError):
            utils.calc_avail_dispatchers(4, 8)

        # exact division
        self.assertEqual(utils.calc_avail_dispatchers(8, 4), 1)

    def test_join_with_common_index(self):
        df1 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=[0, 1, 2])
        df2 = DataFrame({"C": [7, 8, 9], "D": [10, 11, 12]}, index=[0, 1, 2])
        result = utils.join_pandas_dataframes(df1, df2)
        expected = DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9], "D": [10, 11, 12]},
            index=[0, 1, 2],
        )
        assert_frame_equal(result, expected)

    def test_join_with_different_index(self):
        df1 = DataFrame({"A": [1, 2], "B": [3, 4]}, index=[0, 1])
        df2 = DataFrame({"C": [5, 6], "D": [7, 8]}, index=[1, 2])
        result = utils.join_pandas_dataframes(df1, df2)
        expected = DataFrame(
            {
                "A": [1.0, 2.0, None],
                "B": [3.0, 4.0, None],
                "C": [None, 5.0, 6.0],
                "D": [None, 7.0, 8.0],
            },
            index=[0, 1, 2],
        )
        assert_frame_equal(result, expected)

    def test_join_with_overlapping_columns(self):
        df1 = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=[0, 1, 2])
        df2 = DataFrame({"A": [7, 8, 9], "D": [10, 11, 12]}, index=[0, 1, 2])
        result = utils.join_pandas_dataframes(df1, df2)
        expected = DataFrame(
            {"A_l": [1, 2, 3], "B": [4, 5, 6], "A_r": [7, 8, 9], "D": [10, 11, 12]},
            index=[0, 1, 2],
        )
        assert_frame_equal(result, expected)

    def test_execute_command(self):
        cmd_pass = "pwd"
        cmd_fail = "LLL"
        ret_pass = utils.execute_command(cmd_pass)
        ret_fail = utils.execute_command(cmd_fail)

        self.assertEqual(ret_pass, 0)
        self.assertNotEqual(ret_fail, 0)

    @patch("os.system")
    def test_execute_command_mock(self, mock_system):
        mock_system.return_value = 0  # Mock return for system command
        ret = utils.execute_command('echo "Hello World"')
        self.assertEqual(ret, 0)

    def test_multiprocess_dispatch_callable_with_progress(self):
        inputs = list(range(5000))
        res = utils.multiprocess_dispatch(_helper, inputs, 2, True, description="MP")
        # TODO: pytest doesn't show stdout
        total = sum(res)
        self.assertEqual(total, 12522500)

    def test_multiprocess_dispatch_callable_no_progress(self):
        inputs = list(range(5000))
        res = utils.multiprocess_dispatch(_helper, inputs, 2, False, description="MP")
        # TODO: pytest doesn't show stdout
        total = sum(res)
        self.assertEqual(total, 12522500)

    def test_multiprocess_dispatch_system_with_progress(self):
        cmds = []
        for _ in range(100):
            cmds.append("echo TEST")
        res = utils.multiprocess_dispatch("system", cmds, 2, True, description="Echo")
        total = len(res)
        self.assertEqual(total, 100)

    def test_multiprocess_dispatch_system_no_progress(self):
        cmds = []
        for _ in range(100):
            cmds.append("echo TEST")
        res = utils.multiprocess_dispatch("system", cmds, 1, False, description="Echo")
        total = len(res)
        self.assertEqual(total, 100)

    def test_create_diamond_blastp_cmd(self):
        cmd = utils.create_diamond_blastp_cmd(
            Path("query.fasta"),
            Path("db.dmnd"),
            Path("out.txt"),
            "sensitive",
            0.001,
            4,
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        )
        expected_cmd = (
            "diamond blastp --query query.fasta --quiet --db db.dmnd --outfmt "
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
            "--out out.txt --evalue 0.001 --threads 4 --sensitive"
        )
        self.assertEqual(cmd, expected_cmd)

    def test_create_blastn_cmd(self):
        cmd = utils.create_blastn_cmd(
            fasta_file=Path("query.fasta"),
            database_f=Path("db"),
            out_file=Path("out.txt"),
            blast_evalue=0.001,
            blast_cores=4,
            outfmt="6 qseqid sseqid pident evalue bitscore",
        )
        expected_cmd = "blastn -query query.fasta -db db -outfmt 6 qseqid sseqid pident evalue bitscore -out out.txt -evalue 0.001 -num_threads 4"
        self.assertEqual(cmd, expected_cmd)

    def test_create_blastdb_cmd(self):
        self.assertEqual(
            utils.create_blastdb_cmd(Path("input.fasta"), "nucl"),
            "makeblastdb -dbtype nucl -in input.fasta -out input.fasta > /dev/null",
        )
        self.assertEqual(
            utils.create_blastdb_cmd(Path("input.fasta"), "prot"),
            "makeblastdb -dbtype prot -in input.fasta -out input.fasta > /dev/null",
        )
        self.assertEqual(utils.create_blastdb_cmd(Path("input.fasta"), "invalid"), "")

    def test_create_diamonddb_cmd(self):
        self.assertEqual(
            utils.create_diamonddb_cmd(Path("input.fasta")),
            "diamond makedb --in input.fasta --db input.fasta > /dev/null",
        )

    @patch("pathlib.Path.exists")
    def test_check_if_dir_exists(self, mock_exists):
        mock_exists.return_value = True
        self.assertTrue(utils.check_if_dir_exists(Path("some_dir")))

        mock_exists.return_value = False
        self.assertFalse(utils.check_if_dir_exists(Path("some_dir")))

    @patch("pathlib.Path.exists")
    def test_check_if_file_exists(self, mock_exists):
        mock_exists.return_value = True
        self.assertTrue(utils.check_if_file_exists(Path("some_file")))

        mock_exists.return_value = False
        self.assertFalse(utils.check_if_file_exists(Path("some_file")))

    @patch("zipfile.ZipFile")
    @patch("gzip.GzipFile")
    @patch("context.utils.check_if_file_exists")
    def test_unzip_file(self, mock_check_exists, mock_gzipfile, mock_zipfile):
        mock_check_exists.return_value = True

        file = Path("somefile.zip")
        with self.assertRaises(ValueError):
            utils.unzip_file(file, "rar")

        mock_check_exists.return_value = False
        with self.assertRaises(FileNotFoundError):
            utils.unzip_file(file, "zip")

    @patch("tempfile.mkstemp")
    def test_create_temporary_file(self, mock_mkstemp):
        temp_path = utils.create_temporary_file()
        self.assertEqual(temp_path.exists(), True)

    def test_dict_to_dataframe(self):
        d = {"a": 1, "b": 2}
        df = utils.dict_to_dataframe(d)
        expected_df = DataFrame.from_dict(d, orient="index")
        pd.testing.assert_frame_equal(df, expected_df)

    @patch("requests.get")
    def test_download_file(self, mock_get):
        mock_response = MagicMock()
        mock_get.return_value = mock_response
        mock_response.iter_content.return_value = [b"1234", b"5678"]

        size = utils.download_file("http://example.com", "output.txt")
        self.assertEqual(size, 8)
        Path("output.txt").unlink()

    @patch("requests.head")
    def test_get_remote_file_size(self, mock_head):
        mock_response = MagicMock()
        mock_head.return_value = mock_response
        mock_response.headers = {"content-length": "1234"}
        size = utils.get_remote_file_size("http://example.com")
        self.assertEqual(size, 1234)

        mock_response.headers = {}
        size = utils.get_remote_file_size("http://example.com")
        self.assertEqual(size, 0)

    # TODO:
    # @patch("Bio.SeqIO.parse")
    # @patch("Bio.Seq")
    # def test_translate_fasta_records(self, mock_seq, mock_seqio_parse):
    #     mock_item = MagicMock()
    #     mock_item.translate.return_value = "M"
    #     mock_item.seq = "ATG"
    #     mock_seqio_parse.return_value = [mock_item]
    #
    #     records = utils.translate_fasta_records(Path("fakepath"))
    #     self.assertEqual(len(records), 1)
    #     self.assertEqual(records[0].seq, "ATGC")

    #
    #     @patch('Bio.SeqIO.write')
    #     def test_seqrecords_to_fasta(self, mock_seqio_write):
    #         records = [MagicMock(id="test")]
    #         utility_module.seqrecords_to_fasta(records, Path("output.fasta"))
    #         mock_seqio_write.assert_called_once_with(records, Path("output.fasta"), "fasta")
    #
