import unittest
from unittest.mock import patch, call, MagicMock
from pathlib import Path
from context import cli
import pytest


class TestRunFunctions(unittest.TestCase):
    @patch("context.cli.SpeciesDemarcator")
    @patch("context.cli.validate_directory", return_value=True)
    def test_run_species_demarcation(
        self, mock_validate_directory, mock_SpeciesDemarcator
    ):
        args = {
            "in": "input_directory",
            "o": "output_directory",
            "fastani_cores": 4,
            "debug": True,
            "kmer": 16,
            "fraglen": 2000,
            "minfraction": 0.7,
            "inflation": 2.0,
            "mcl_cores": 4,
        }

        mock_instance = mock_SpeciesDemarcator.return_value

        cli.run_species_demarcation(args)

        mock_validate_directory.assert_called_once_with(Path("input_directory"))
        mock_SpeciesDemarcator.assert_called_once_with(
            in_dir=Path("input_directory"),
            out_dir=Path("output_directory"),
            fastani_cores=4,
            kmer=16,
            fraglen=2000,
            minfrac=0.7,
            inflation=2.0,
            mcl_cores=4,
            debug=True,
        )
        mock_instance.assign_species.assert_called_once()

    @patch("context.cli.Orthologues_identifier")
    def test_run_orthologues(self, mock_Orthologues_identifier):
        args = {
            "fasta_files": ["file1", "file2"],
            "out": "output_directory",
            "ref": "reference_strain",
            "ref_list": None,
            "input_type": "prot",
            "cores": 4,
            "evalue": 1e-10,
            "dmnd_sensitivity": "very_sensitive",
            "no_filter_orthologues": False,
            "no_concurrent": False,
        }

        mock_instance = mock_Orthologues_identifier.return_value

        cli.run_orthologues(args)

        mock_Orthologues_identifier.assert_called_once_with(
            fasta_files_list=["file1", "file2"],
            out_dir=Path("output_directory"),
            ref="reference_strain",
            ref_list=None,
            input_type="prot",
            available_cores=cli.software_defaults["system"]["usable_cores"],
            blast_cores=4,
            concurrent=True,
            evalue=1e-10,
            dmnd_sensitivity="very_sensitive",
            no_filter=False,
        )
        mock_instance.calculate_orthologues.assert_called_once()

    @patch("context.cli.read_excel")
    @patch("context.cli.Core_identifier")
    @patch(
        "context.cli.tqdm",
        return_value=[Path("ortholog_matrix_1"), Path("ortholog_matrix_2")],
    )
    def test_run_core(self, mock_tqdm, mock_Core_identifier, mock_read_excel):
        args = {
            "out": "output_directory",
            "species": "species_file.xlsx",
            "in": None,
            "ref_list": "ref_list.txt",
            "core_perc": 0.95,
        }

        mock_read_excel.return_value = MagicMock()
        mock_instance = mock_Core_identifier.return_value

        with patch(
            "builtins.open",
            unittest.mock.mock_open(read_data="ortholog_matrix_1\northolog_matrix_2"),
        ) as mock_file:
            cli.run_core(args)

        mock_read_excel.assert_called_once_with(Path("species_file.xlsx"), index_col=0)
        mock_Core_identifier.assert_has_calls(
            [
                call(
                    orthology_fin=Path("ortholog_matrix_1"),
                    out_dir=Path("output_directory"),
                    species_df=mock_read_excel.return_value,
                    core_perc=0.95,
                ),
                call().calculate_core(),
                call(
                    orthology_fin=Path("ortholog_matrix_2"),
                    out_dir=Path("output_directory"),
                    species_df=mock_read_excel.return_value,
                    core_perc=0.95,
                ),
                call().calculate_core(),
            ]
        )

    @patch("context.cli.Phylogenomic")
    def test_run_phylogenomic(self, mock_Phylogenomic):
        args = {
            "in": "orthologies_matrix.xlsx",
            "out": "output_directory",
            "cores": 8,
            "no_keep_fasta": True,
            "tree_model": "LG+G4",
            "method": "IQTree",
            "input_type": "prot",
            "debug": False,
            "fasta_files": ["file1.fasta", "file2.fasta"],
        }

        mock_instance = mock_Phylogenomic.return_value

        cli.run_phylogenomic(args)

        mock_Phylogenomic.assert_called_once_with(
            orthology_matrix_in=Path("orthologies_matrix.xlsx"),
            cores=8,
            fasta_files_list=["file1.fasta", "file2.fasta"],
            out_dir=Path("output_directory"),
            no_keep_fasta=True,
            tree_method="IQTree",
            iqtree_model="LG+G4",
            input_type="prot",
            debug=False,
        )
        mock_instance.run_phylogenomic.assert_called_once()

    @patch("context.cli.EggNOGRunner")
    @patch("context.cli.EggNOGParser")
    @patch("context.cli.DB_BASE_DIR", new_callable=Path)
    def test_run_eggnog(self, mock_DB_BASE_DIR, mock_EggNOGParser, mock_EggNOGRunner):
        args = {
            "debug": True,
            "out": "output_directory",
            "cores": 4,
            "pident": "0.0",
            "qcov": "0.0",
            "scov": "0.0",
            "input_type": "prot",
            "db": None,
            "fasta_files": ["file1.fasta", "file2.fasta"],
            "in": None,
            "in_list": None,
        }

        mock_runner_instance = mock_EggNOGRunner.return_value
        mock_parser_instance = mock_EggNOGParser.return_value

        cli.run_eggnog(args)

        mock_EggNOGRunner.assert_called_once_with(
            fasta_files_list=["file1.fasta", "file2.fasta"],
            out_dir=Path("output_directory"),
            cores=4,
            pident="0.0",
            qcov="0.0",
            scov="0.0",
            input_type="prot",
            database_dir=mock_DB_BASE_DIR / "site-packages" / "data",
            debug=True,
        )
        mock_runner_instance.execute_eggnog_mapper.assert_called_once()
        mock_EggNOGParser.assert_called_once_with(
            fasta_files_list=["file1.fasta", "file2.fasta"],
            core_proteins_table_f=None,
            out_dir=Path("output_directory"),
            debug=True,
        )
        mock_parser_instance.gather_eggnog_proteome_results.assert_called_once()
        self.assertEqual(
            mock_parser_instance.compare_proteome_to_core_and_fps.call_count, 0
        )

        # Test with core_proteins_file and core_protein_files_reflist
        args["in"] = "core_proteins_file.xlsx"
        args["in_list"] = None

        cli.run_eggnog(args)

        mock_EggNOGParser.assert_called_with(
            fasta_files_list=["file1.fasta", "file2.fasta"],
            core_proteins_table_f=[Path("core_proteins_file.xlsx")],
            out_dir=Path("output_directory"),
            debug=True,
        )
        mock_parser_instance.compare_proteome_to_core_and_fps.assert_called_once()

    @patch("context.cli.smBGCLocalRunner")
    @patch("context.cli.smBGCParser")
    def test_run_smbgc(self, mock_smBGCParser, mock_smBGCLocalRunner):
        args = {
            "debug": True,
            "out": "output_directory",
            "strictness": "strict",
            "genefinding_tool": "glimmer",
            "cores": 8,
            "db": None,
            "no_concurrent": False,
            "fasta_files": ["file1.fasta", "file2.fasta"],
        }

        mock_instance_runner = mock_smBGCLocalRunner.return_value
        mock_instance_parser = mock_smBGCParser.return_value

        cli.run_smbgc(args)

        mock_smBGCLocalRunner.assert_called_once_with(
            genome_fasta_files=["file1.fasta", "file2.fasta"],
            out_dir=Path("output_directory"),
            strictness="strict",
            genefinding_tool="glimmer",
            database_dir=None,
            cores=8,
            available_cores=cli.software_defaults["system"]["usable_cores"],
            concurrent=True,
            debug=True,
        )
        mock_instance_runner.analyze_genomes.assert_called_once()
        mock_smBGCParser.assert_called_once_with(Path("output_directory"), 8)
        mock_instance_parser.gather_results.assert_called_once()

    @patch("context.cli.VF_analyzer")
    def test_run_virulence(self, mock_VF_analyzer):
        args = {
            "out": "output_directory",
            "input_type": "prot",
            "cores": 8,
            "evalue": 1e-5,
            "dmnd_sensitivity": "very_sensitive",
            "debug": False,
            "db": None,
            "no_concurrent": False,
            "fasta_files": ["file1.fasta", "file2.fasta"],
        }

        mock_instance = mock_VF_analyzer.return_value

        cli.run_virulence(args)

        mock_VF_analyzer.assert_called_once_with(
            fasta_files_list=["file1.fasta", "file2.fasta"],
            database_dir=Path(cli.DB_BASE_DIR) / "Virulence",
            out_dir=Path("output_directory"),
            input_type="prot",
            available_cores=cli.software_defaults["system"]["usable_cores"],
            blast_cores=8,
            evalue=1e-5,
            dmnd_sensitivity="very_sensitive",
            concurrent=True,
            debug=False,
        )
        mock_instance.execute_homology_search.assert_called_once()
        mock_instance.parse_results.assert_called_once()

    @patch("context.cli.CAZY_analyzer")
    def test_run_cazy(self, mock_CAZY_analyzer):
        args = {
            "out": "output_directory",
            "input_type": "prot",
            "cores": 8,
            "evalue": 1e-5,
            "db": None,
            "no_concurrent": False,
            "debug": True,
            "fasta_files": ["file1.fasta", "file2.fasta"],
        }

        mock_instance = mock_CAZY_analyzer.return_value

        cli.run_cazy(args)

        mock_CAZY_analyzer.assert_called_once_with(
            fasta_files_list=["file1.fasta", "file2.fasta"],
            out_dir=Path("output_directory"),
            database_dir=Path(cli.DB_BASE_DIR) / "CAZY",
            input_type="prot",
            available_cores=cli.software_defaults["system"]["usable_cores"],
            cores=8,
            concurrent=True,
            evalue=1e-5,
            debug=True,
        )
        mock_instance.execute_cazy_search.assert_called_once()

    @patch("context.cli.AMR_analyzer")
    def test_run_amr(self, mock_AMR_analyzer):
        args = {
            "out": "output_directory",
            "input_type": "prot",
            "cores": 8,
            "db": None,
            "no_concurrent": False,
            "debug": True,
            "fasta_files": ["file1.fasta", "file2.fasta"],
        }

        mock_instance = mock_AMR_analyzer.return_value

        cli.run_amr(args)

        mock_AMR_analyzer.assert_called_once_with(
            fasta_files_list=["file1.fasta", "file2.fasta"],
            database_dir=Path(cli.DB_BASE_DIR) / "AMR",
            out_dir=Path("output_directory"),
            input_type="prot",
            available_cores=cli.software_defaults["system"]["usable_cores"],
            cores=8,
            concurrent=True,
            debug=True,
        )
        mock_instance.search_amr.assert_called_once()
        mock_instance.parse_results.assert_called_once()

    @patch("context.cli.GenomeDownloader")
    @patch("context.cli.validate_directory", return_value=True)
    @patch("context.cli.read_excel")
    @patch("context.cli.utils.join_pandas_dataframes")
    @patch("context.cli.utils.check_if_dir_exists", return_value=True)
    @patch("context.cli.utils.check_if_dir_is_empty", return_value=False)
    @patch("context.cli.SpeciesDemarcator")
    def test_run_download(
        self,
        mock_SpeciesDemarcator,
        mock_check_if_dir_exists,
        mock_check_if_dir_is_empty,
        mock_join_pandas_dataframes,
        mock_read_excel,
        mock_validate_directory,
        mock_GenomeDownloader,
    ):
        args = {
            "out": "output_directory",
            "taxon": "Bacteria",
            "level": "complete",
            "source": "RefSeq",
            "keep_plasmids": False,
            "keep_download": False,
            "perform_fastani": True,
            "debug": True,
        }

        downloader_instance = mock_GenomeDownloader.return_value
        demarcator_instance = mock_SpeciesDemarcator.return_value
        mock_join_pandas_dataframes.return_value = MagicMock()

        cli.run_download(args)

        mock_GenomeDownloader.assert_called_once_with(
            taxon="Bacteria",
            out_dir=Path("output_directory"),
            assembly_level="complete",
            assembly_source="RefSeq",
            keep_plasmids=False,
            debug=True,
        )
        downloader_instance.download_hydrated.assert_called_once()
        downloader_instance.extract_dataset_zip.assert_called_once()
        downloader_instance.create_output_directories.assert_called_once()
        downloader_instance.process_gbff_files.assert_called_once()
        downloader_instance.write_annotations.assert_called_once()
        downloader_instance.write_16S_fasta.assert_called_once()
        downloader_instance.remove_dataset_archive.assert_called_once()

        mock_SpeciesDemarcator.assert_called_once_with(
            in_dir=Path("output_directory") / "Genomic_fasta_files",
            out_dir=Path("output_directory"),
            fastani_cores=cli.software_defaults["species_demarcation"]["cores"],
            kmer=cli.software_defaults["species_demarcation"]["kmer"],
            fraglen=cli.software_defaults["species_demarcation"]["fraglen"],
            minfrac=cli.software_defaults["species_demarcation"]["minfrac"],
            inflation=cli.software_defaults["species_demarcation"]["inflation"],
            mcl_cores=cli.software_defaults["species_demarcation"]["cores"],
            debug=True,
        )
        demarcator_instance.assign_species.assert_called_once()
        mock_read_excel.assert_any_call(
            Path("output_directory") / "Genome_information.xlsx", index_col=0
        )
        mock_read_excel.assert_any_call(
            Path("output_directory")
            / "Species_demarcation"
            / "FastANI_species_clusters.xlsx",
            index_col=0,
        )
        mock_join_pandas_dataframes.assert_called_once()

    @patch("context.cli.SMBGC_installer")
    @patch("context.cli.EGGNOG_installer")
    @patch("context.cli.CAZY_installer")
    @patch("context.cli.AMR_installer")
    @patch("context.cli.VF_installer")
    def test_run_databases_invalid(
        self,
        mock_VF_installer,
        mock_AMR_installer,
        mock_CAZY_installer,
        mock_EGGNOG_installer,
        mock_SMBGC_installer,
    ):
        args = {
            "db": "all,smbgc",
            "db_dir": "custom_db_dir",
            "debug": True,
            "cores": 8,
        }

        # Test with invalid database
        with self.assertRaises(ValueError):
            cli.run_databases(args)

    @pytest.mark.skip(reason="Not implemented correctly")
    @patch("context.cli.SMBGC_installer")
    @patch("context.cli.EGGNOG_installer")
    def test_run_databases_valid_combination(
        self,
        mock_EGGNOG_installer,
        mock_SMBGC_installer,
    ):
        args = {
            "db": "eggnog,smbgc",
            "db_dir": "custom_db_dir",
            "debug": False,
            "cores": 8,
        }

        cli.run_databases(args)

        # Test with all
        mock_SMBGC_installer.assert_called_once_with(
            database_dir="custom_db_dir", debug=True
        )
        # mock_SMBGC_installer.return_value.install_database.assert_called_once()
        mock_EGGNOG_installer.assert_called_once_with(
            database_dir="custom_db_dir", debug=True
        )
        # mock_EGGNOG_installer.return_value.install_database.assert_called_once()

    @pytest.mark.skip(reason="Not implemented correctly")
    @patch("context.cli.SMBGC_installer")
    @patch("context.cli.EGGNOG_installer")
    @patch("context.cli.CAZY_installer")
    @patch("context.cli.AMR_installer")
    @patch("context.cli.VF_installer")
    def test_run_databases_all(
        self,
        mock_VF_installer,
        mock_AMR_installer,
        mock_CAZY_installer,
        mock_EGGNOG_installer,
        mock_SMBGC_installer,
    ):
        args = {
            "db": "all",
            "db_dir": "custom_db_dir",
            "debug": True,
            "cores": 8,
        }

        cli.run_databases(args)

        # Test with all
        mock_SMBGC_installer.assert_called_once_with(
            database_dir="custom_db_dir", debug=True
        )
        mock_SMBGC_installer.return_value.install_database.assert_called_once()
        mock_EGGNOG_installer.assert_called_once_with(
            database_dir="custom_db_dir", debug=True
        )
        mock_EGGNOG_installer.return_value.install_database.assert_called_once()
        mock_CAZY_installer.assert_called_once_with(
            database_dir="custom_db_dir", cores=8, debug=True
        )
        mock_CAZY_installer.return_value.install_database.assert_called_once()
        mock_AMR_installer.assert_called_once_with(
            database_dir="custom_db_dir", debug=True
        )
        mock_AMR_installer.return_value.install_database.assert_called_once()
        mock_VF_installer.assert_called_once_with(
            database_dir="custom_db_dir", debug=True
        )
        mock_VF_installer.return_value.install_database.assert_called_once()

    @patch("context.cli.WorkflowRunner")
    def test_run_workflow(self, mock_WorkflowRunner):
        args = {
            "in": "parameters_file.xlsx",
            "out": "output_directory",
            "debug": True,
            "db": "custom_db_dir",
        }

        workflow_instance = mock_WorkflowRunner.return_value

        cli.run_workflow(args)

        mock_WorkflowRunner.assert_called_once_with(
            out_dir=Path("output_directory"),
            param_file=Path("parameters_file.xlsx"),
            debug=True,
        )
        workflow_instance.read_param_file.assert_called_once()
        workflow_instance.validate_parameters.assert_called_once()
        workflow_instance.identify_tasks.assert_called_once()
        workflow_instance.gather_cds_or_protein_fasta_files.assert_called_once()
        workflow_instance.gather_genomic_fasta_files.assert_called_once()
        workflow_instance.get_per_group_representatives_files.assert_called_once()
