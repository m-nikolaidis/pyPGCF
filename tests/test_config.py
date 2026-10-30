import unittest
from pathlib import Path
from unittest.mock import patch

from pypgcf import config


class TestDatabaseDir(unittest.TestCase):
    def test_cli_value_wins_over_environment(self):
        with patch.dict("os.environ", {config.DB_DIR_ENV: "/env"}):
            self.assertEqual(
                config.database_dir("cazy", "/cli"), Path("/cli/CAZY")
            )

    def test_environment_used_when_no_cli_value(self):
        with patch.dict("os.environ", {config.DB_DIR_ENV: "/env"}):
            self.assertEqual(config.database_dir("amr"), Path("/env/AMR"))
            self.assertEqual(config.database_dir("smbgc"), Path("/env/smBGC"))

    def test_defaults_keep_earlier_install_locations(self):
        with patch.dict("os.environ", {}, clear=True):
            self.assertEqual(
                config.database_dir("virulence"), config.DB_BASE_DIR / "Virulence"
            )
            self.assertIsNone(config.database_dir("smbgc"))
            self.assertEqual(
                config.database_dir("eggnog"),
                config.DB_BASE_DIR / "site-packages" / "data",
            )


if __name__ == "__main__":
    unittest.main()
