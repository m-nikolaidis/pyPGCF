from pathlib import Path
import unittest
from unittest.mock import patch, call
import tempfile

from pandas.testing import assert_frame_equal
from pandas import DataFrame

from context import utils


class TestModule(unittest.TestCase):
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

    def test_calc_avail_dispatchers_basic(self):
        self.assertEqual(utils.calc_avail_dispatchers(8, 2), 3)

    def test_calc_avail_dispatchers_avoid_throttle(self):
        self.assertEqual(utils.calc_avail_dispatchers(8, 2, avoid_throttle=True), 3)
        self.assertEqual(utils.calc_avail_dispatchers(8, 2, avoid_throttle=False), 4)

    def test_calc_avail_dispatchers_usemax(self):
        self.assertEqual(
            utils.calc_avail_dispatchers(8, 2, avoid_throttle=True, usemax=True), 4
        )

    def test_calc_avail_dispatchers_more_cores_per_job_than_available(self):
        with self.assertRaises(ValueError):
            utils.calc_avail_dispatchers(4, 8)

    def test_calc_avail_dispatchers_exact_division(self):
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
