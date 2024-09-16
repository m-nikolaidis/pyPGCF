import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
from context import databases, config
import pytest

# class TestCAZYBuilder(unittest.TestCase):
#     def setUp(self):
#         self.database_dir = Path("/tmp/fake/path/to/database")
#         self.verbose = False
#         self.cores = 4
#         self.cazy_builder = cazy.CAZY_builder(
#             database_dir=self.database_dir, cores=self.cores, verbose=self.verbose
#         )
#
#     @patch("context.cazy.execute_command")
#     def test_setup_valid(self, mock_execute_command):
#         mock_execute_command.return_value = 0
#         self.cazy_builder.setup()
#         build_cmd = f"dbcan_build --cores 4 --db-dir {self.database_dir}/CAZY --clean"
#         mock_execute_command.assert_called_once_with(build_cmd)
#
#     @patch("context.cazy.execute_command")
#     def test_setup_invalid(self, mock_execute_command):
#         mock_execute_command.side_effect = Exception("Command failed")
#         with self.assertRaises(Exception):
#             self.cazy_builder.setup()
#
#     def test_validate_if_database_exists(self):
#         # This is a placeholder test since the method is not implemented
#         self.assertIsNone(self.cazy_builder.validate_if_database_exists())


@pytest.mark.skip(reason="Command doesn't exist in test environment")
class TestCAZYInstaller(unittest.TestCase):
    @patch("context.utils.execute_command")
    @patch("pathlib.Path.mkdir")
    @patch("os.system")
    def test_install_database(self, mock_system, mock_mkdir, mock_execute_command):
        mock_execute_command.return_value = 0
        installer = databases.CAZY_installer(database_dir=None, cores=4, debug=True)
        _ = installer.install_database()
        mock_system.assert_called_with(
            f"dbcan_build --cores 4 --db-dir {config.DB_BASE_DIR}/CAZY --clean"
        )
        mock_mkdir.assert_called_with(exist_ok=True, parents=True)


@pytest.mark.skip(
    reason="Will add in next patch. Command doesn't exist in test environment"
)
class TestVFInstaller(unittest.TestCase):
    @patch("context.utils.execute_command")
    @patch("context.utils.download_file")
    @patch("context.utils.get_remote_file_size")
    @patch("pathlib.Path.mkdir")
    @patch("Bio.SeqIO.parse", return_value=[])
    @patch("builtins.open", new_callable=unittest.mock.mock_open)
    def test_install_database(
        self,
        mock_open,
        mock_SeqIO_parse,
        mock_mkdir,
        mock_get_remote_file_size,
        mock_download_file,
        mock_execute_command,
    ):
        mock_download_file.return_value = 1024
        mock_get_remote_file_size.return_value = 1024
        mock_execute_command.return_value = 0
        # mock_file_handle = mock_open.return_value.__enter__.return_value

        installer = databases.VF_installer(database_dir=Path("/tmp/DB/"), debug=True)
        installer.install_database()

        # Checking directory creation
        mock_mkdir.assert_called_with(exist_ok=True, parents=True)

        # Check file download and size validation
        # self.assertEqual(mock_download_file.call_count, 2)
        # self.assertEqual(mock_get_remote_file_size.call_count, 2)
        #
        # # Check file write operations
        # mock_file_handle.write.assert_any_call(
        #     "# Database downloaded and built: {}\n".format(
        #         datetime.now().strftime("%D")
        #     )
        # )
        #
        # # Ensure command execution for both Prot and Nucl
        # self.assertEqual(mock_execute_command.call_count, 2)


#
#
@pytest.mark.skip(
    reason="Will add in next patch. Command doesn't exist in test environment"
)
class TestAMRInstaller(unittest.TestCase):
    @patch("pypgcf.utils.execute_command")
    @patch("pypgcf.utils.download_file")
    @patch("pypgcf.utils.get_remote_file_size")
    @patch("pathlib.Path.mkdir")
    def test_install_database(
        self,
        mock_mkdir,
        mock_get_remote_file_size,
        mock_download_file,
        mock_execute_command,
    ):
        mock_download_file.return_value = 1024
        mock_get_remote_file_size.return_value = 1024
        mock_execute_command.return_value = 0

        installer = databases.AMR_installer(database_dir=None, debug=True)
        installer.install_database()

        # Checking directory creation
        mock_mkdir.assert_called_with(exist_ok=True, parents=True)

        # Check file downloads for each file
        self.assertEqual(mock_download_file.call_count, len(installer.files))
        self.assertEqual(mock_get_remote_file_size.call_count, len(installer.files))

        # Ensure database indexing command is executed correctly
        self.assertEqual(
            mock_execute_command.call_count,
            len(
                [
                    f
                    for f in installer.files
                    if not (
                        f.endswith(".tab")
                        or f.endswith(".txt")
                        or f.endswith(".LIB")
                        or f.endswith("suppress")
                    )
                ]
            ),
        )


@pytest.mark.skip(
    reason="Will add in next patch. Command doesn't exist in test environment"
)
class TestSMBGCInstaller(unittest.TestCase):
    @patch("pypgcf.utils.execute_command")
    @patch("pathlib.Path.mkdir")
    def test_install_database(self, mock_mkdir, mock_execute_command):
        mock_execute_command.return_value = 0

        installer = databases.SMBGC_installer(database_dir=None, debug=True)
        installer.install_database()

        # Verify that the command was built correctly
        mock_execute_command.assert_called_with("download-antismash-databases")

        self.assertEqual(mock_execute_command.call_count, 1)
        if installer.debug:
            self.assertIsNone(installer.install_database())


@pytest.mark.skip(
    reason="Will add in next patch. Command doesn't exist in test environment"
)
class TestEGGNOGInstaller(unittest.TestCase):
    @patch("pypgcf.utils.execute_command")
    @patch("pathlib.Path.mkdir")
    def test_install_database(self, mock_mkdir, mock_execute_command):
        mock_execute_command.return_value = 0

        installer = databases.EGGNOG_installer(database_dir=None, debug=True)
        installer.install_database()

        # Verify the creation of the correct command string
        cmd = f"download_eggnog_data.py --data_dir {config.DB_BASE_DIR / 'site-packages' / 'data'} -y "
        mock_execute_command.assert_called_with(cmd)

        self.assertEqual(mock_execute_command.call_count, 1)
