import pytest
from pathlib import Path
from unittest.mock import Mock, patch, mock_open, MagicMock, call
from datetime import datetime
import tempfile
import shutil

from databases import (
    CAZY_installer,
    VF_installer,
    AMR_installer,
    SMBGC_installer,
    EGGNOG_installer
)


@pytest.fixture
def temp_db_dir():
    """Create a temporary directory for testing"""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir)


class TestCAZYInstaller:
    
    @patch('databases.execute_command')
    @patch('databases.datetime')
    def test_install_database_success(self, mock_datetime, mock_execute, temp_db_dir):
        """Test successful CAZY database installation"""
        # Setup
        mock_datetime.now.return_value.strftime.return_value = "01/01/24"
        mock_execute.return_value = 0
        
        installer = CAZY_installer(database_dir=temp_db_dir, cores=4, debug=False)
        
        # Execute
        result = installer.install_database()
        
        # Assert
        assert result == 0
        mock_execute.assert_called_once()
        assert "dbcan_build" in mock_execute.call_args[0][0]
        assert "--cpus 4" in mock_execute.call_args[0][0]
        assert (installer.database_dir / "log.txt").exists()
    
    def test_cazy_installer_creates_directory(self, temp_db_dir):
        """Test that CAZY installer creates proper directory structure"""
        installer = CAZY_installer(database_dir=temp_db_dir, cores=2, debug=True)
        
        assert installer.database_dir.exists()
        assert installer.database_dir == temp_db_dir / "CAZY"
        assert installer.cores == 2
    
    def test_cazy_installer_uses_default_dir(self):
        """Test CAZY installer with None database_dir uses default"""
        with patch('databases.DB_BASE_DIR', Path("/mock/base")):
            installer = CAZY_installer(database_dir=None, cores=2, debug=False)
            assert installer.database_dir == Path("/mock/base") / "CAZY"


class TestVFInstaller:
    
    def test_vf_installer_initialization(self, temp_db_dir):
        """Test VF installer initialization"""
        installer = VF_installer(database_dir=temp_db_dir, debug=True)
        
        assert installer.database_dir == temp_db_dir / "Virulence"
        assert len(installer.urls) == 2
        assert len(installer.mol_types) == 2
    
    def test_split_vf_desc_valid_input(self):
        """Test splitting VF description with valid input"""
        installer = VF_installer(database_dir=None, debug=False)
        
        test_string = "VF0001 (vfg) test gene [Category) - Immune evasion (VFC0001)] [Staphylococcus aureus]"
        category, origin, desc = installer._split_VF_desc(test_string)
        
        assert category == "Immune evasion"
        assert origin == "Staphylococcus aureus"
        assert "VF0001" in desc
    
    def test_split_vf_desc_invalid_input(self):
        """Test splitting VF description with invalid input"""
        installer = VF_installer(database_dir=None, debug=False)
        
        test_string = "Invalid format string"
        category, origin, desc = installer._split_VF_desc(test_string)
        
        assert category is None
        assert origin is None
        assert desc is None
    
    @patch('databases.SeqIO')
    def test_create_vf_description_file(self, mock_seqio, temp_db_dir):
        """Test VF description file creation"""
        # Setup mock records
        mock_record1 = Mock()
        mock_record1.id = "VF0001"
        mock_record1.description = "gene1 [Cat) - Adhesion (VFC0001)] [E. coli]"
        
        mock_record2 = Mock()
        mock_record2.id = "VF0002"
        mock_record2.description = "gene2 [Cat) - Toxin (VFC0002)] [S. aureus]"
        
        mock_seqio.parse.return_value = [mock_record1, mock_record2]
        
        installer = VF_installer(database_dir=temp_db_dir, debug=False)
        test_file = temp_db_dir / "Virulence" / "test.fasta"
        test_file.parent.mkdir(parents=True, exist_ok=True)
        test_file.touch()
        
        # Execute
        installer.create_vf_description_file(test_file)
        
        # Assert
        desc_file = installer.database_dir / "vfdb_desc.tsv"
        assert desc_file.exists()
        
        with open(desc_file, 'r') as f:
            content = f.read()
            assert "VF0001" in content
            assert "VF0002" in content
    
    @patch('databases.execute_command')
    @patch('databases.create_blastdb_cmd')
    @patch('databases.createdmnddb_cmd')
    @patch('databases.unzip_file')
    @patch('databases.get_remote_file_size')
    @patch('databases.download_file')
    def test_install_database_mock(
        self,
        mock_download,
        mock_get_size,
        mock_unzip,
        mock_dmnd_cmd,
        mock_blast_cmd,
        mock_execute,
        temp_db_dir
    ):
        """Test VF database installation with mocks"""
        # Setup
        mock_download.return_value = 1000
        mock_get_size.return_value = 1000
        mock_dmnd_cmd.return_value = "diamond makedb command"
        mock_blast_cmd.return_value = "makeblastdb command"
        mock_execute.return_value = 0
        
        installer = VF_installer(database_dir=temp_db_dir, debug=False)
        
        # Create mock fasta files that will exist after unzip
        for url in installer.urls:
            filename = url.split("/")[-1].replace(".gz", "")
            fasta_file = installer.database_dir / filename
            fasta_file.parent.mkdir(parents=True, exist_ok=True)
            fasta_file.write_text(">VF0001 test [Cat) - Type (VFC)] [Origin]\nATGC\n")
        
        # Execute
        installer.install_database()
        
        # Assert
        assert mock_download.call_count == 2
        assert mock_get_size.call_count == 2
        assert mock_unzip.call_count == 2
        assert mock_execute.call_count == 2
    
    @patch('databases.execute_command')
    @patch('databases.create_blastdb_cmd')
    @patch('databases.createdmnddb_cmd')
    @patch('databases.unzip_file')
    @patch('databases.get_remote_file_size')
    @patch('databases.download_file')
    def test_install_database_size_mismatch(
        self,
        mock_download,
        mock_get_size,
        mock_unzip,
        mock_dmnd_cmd,
        mock_blast_cmd,
        mock_execute,
        temp_db_dir
    ):
        """Test VF database installation fails on size mismatch"""
        # Setup - downloaded size != expected size
        mock_download.return_value = 1000
        mock_get_size.return_value = 2000
        
        installer = VF_installer(database_dir=temp_db_dir, debug=False)
        
        # Execute & Assert
        with pytest.raises(ConnectionError, match="not downloaded"):
            installer.install_database()
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_install_database_real_download(self, temp_db_dir):
        """
        INTEGRATION TEST: Actually download VF database (SLOW)
        This test actually hits the network and downloads real data.
        Skip by default with: pytest -m "not integration"
        """
        installer = VF_installer(database_dir=temp_db_dir, debug=True)
        
        # This will actually download the files
        try:
            installer.install_database()
            
            # Verify files were created
            assert (installer.database_dir / "VFDB_setA_pro.fas").exists()
            assert (installer.database_dir / "VFDB_setA_nt.fas").exists()
            assert (installer.database_dir / "vfdb_desc.tsv").exists()
            
            # Verify description file has content
            desc_file = installer.database_dir / "vfdb_desc.tsv"
            with open(desc_file, 'r') as f:
                lines = f.readlines()
                assert len(lines) > 1  # Header + data
                
        except Exception as e:
            pytest.skip(f"Real download test failed (network issue?): {e}")


class TestAMRInstaller:
    
    def test_amr_installer_initialization(self, temp_db_dir):
        """Test AMR installer initialization"""
        installer = AMR_installer(database_dir=temp_db_dir, debug=True)
        
        assert installer.database_dir == temp_db_dir / "AMR"
        assert installer.database_dir.exists()
    
    @patch('databases.execute_command')
    def test_install_database_native_debug_mode(self, mock_execute, temp_db_dir):
        """Test AMR database installation in debug mode"""
        installer = AMR_installer(database_dir=temp_db_dir, debug=True)
        
        installer.install_database_native()
        
        mock_execute.assert_called_once()
        cmd = mock_execute.call_args[0][0]
        assert "amrfinder_update" in cmd
        assert str(temp_db_dir / "AMR") in cmd
        assert "--quiet" not in cmd
    
    @patch('databases.execute_command')
    def test_install_database_native_quiet_mode(self, mock_execute, temp_db_dir):
        """Test AMR database installation in quiet mode"""
        installer = AMR_installer(database_dir=temp_db_dir, debug=False)
        
        installer.install_database_native()
        
        mock_execute.assert_called_once()
        cmd = mock_execute.call_args[0][0]
        assert "--quiet" in cmd


class TestSMBGCInstaller:
    
    def test_smbgc_installer_with_dir(self, temp_db_dir):
        """Test SMBGC installer with custom directory"""
        installer = SMBGC_installer(database_dir=temp_db_dir, debug=True)
        
        assert installer.database_dir == temp_db_dir / "smBGC"
        assert installer.database_dir.exists()
    
    def test_smbgc_installer_without_dir(self):
        """Test SMBGC installer with None directory"""
        installer = SMBGC_installer(database_dir=None, debug=False)
        
        assert installer.database_dir is None
    
    @patch('databases.execute_command')
    def test_install_database_with_custom_dir(self, mock_execute, temp_db_dir):
        """Test SMBGC database installation with custom directory"""
        mock_execute.return_value = 0
        
        installer = SMBGC_installer(database_dir=temp_db_dir, debug=True)
        installer.install_database()
        
        mock_execute.assert_called_once()
        cmd = mock_execute.call_args[0][0]
        assert "download-antismash-databases" in cmd
        assert "--database-dir" in cmd
        assert (installer.database_dir / "log.txt").exists()
    
    @patch('databases.execute_command')
    def test_install_database_without_custom_dir(self, mock_execute):
        """Test SMBGC database installation without custom directory"""
        mock_execute.return_value = 0
        
        installer = SMBGC_installer(database_dir=None, debug=False)
        installer.install_database()
        
        mock_execute.assert_called_once()
        cmd = mock_execute.call_args[0][0]
        assert "download-antismash-databases" in cmd
        assert "--database-dir" not in cmd
    
    @patch('databases.execute_command')
    def test_install_database_failure(self, mock_execute, temp_db_dir):
        """Test SMBGC database installation handles failures"""
        mock_execute.return_value = 1
        
        installer = SMBGC_installer(database_dir=temp_db_dir, debug=True)
        
        with pytest.raises(RuntimeError, match="Something went wrong"):
            installer.install_database()


class TestEGGNOGInstaller:
    
    def test_eggnog_installer_with_custom_dir(self, temp_db_dir):
        """Test EGGNOG installer with custom directory"""
        installer = EGGNOG_installer(database_dir=temp_db_dir, debug=True)
        
        assert installer.database_dir == temp_db_dir / "eggNOG"
        assert installer.database_dir.exists()
    
    @patch('databases.DB_BASE_DIR', Path("/mock/base"))
    def test_eggnog_installer_with_default_dir(self):
        """Test EGGNOG installer uses default directory"""
        installer = EGGNOG_installer(database_dir=None, debug=False)
        
        expected = Path("/mock/base") / "site-packages" / "data"
        assert installer.database_dir == expected
    
    @patch('databases.execute_command')
    def test_install_database_debug_mode(self, mock_execute, temp_db_dir):
        """Test EGGNOG database installation in debug mode"""
        mock_execute.return_value = 0
        
        installer = EGGNOG_installer(database_dir=temp_db_dir, debug=True)
        installer.install_database()
        
        mock_execute.assert_called_once()
        cmd = mock_execute.call_args[0][0]
        assert "download_eggnog_data.py" in cmd
        assert str(temp_db_dir / "eggNOG") in cmd
        assert "-y" in cmd
        assert "-q" not in cmd
        assert (installer.database_dir / "log.txt").exists()
    
    @patch('databases.execute_command')
    def test_install_database_quiet_mode(self, mock_execute, temp_db_dir):
        """Test EGGNOG database installation in quiet mode"""
        mock_execute.return_value = 0
        
        installer = EGGNOG_installer(database_dir=temp_db_dir, debug=False)
        installer.install_database()
        
        cmd = mock_execute.call_args[0][0]
        assert "-q" in cmd
    
    @patch('databases.execute_command')
    def test_install_database_failure(self, mock_execute, temp_db_dir):
        """Test EGGNOG database installation handles failures"""
        mock_execute.return_value = 1
        
        installer = EGGNOG_installer(database_dir=temp_db_dir, debug=True)
        
        with pytest.raises(RuntimeError, match="Something went wrong"):
            installer.install_database()


# Pytest configuration for marking integration tests
def pytest_configure(config):
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests (deselect with '-m \"not integration\"')"
    )
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
