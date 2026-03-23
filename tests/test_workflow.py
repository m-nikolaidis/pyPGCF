import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pandas import DataFrame

from context import workflow


class TestWorkflowRunner(unittest.TestCase):
    def setUp(self):
        self.out_dir = Path("/path/to/output")
        self.param_file = Path("/path/to/params.xlsx")
        self.runner = workflow.WorkflowRunner(
            out_dir=self.out_dir, param_file=self.param_file, debug=True
        )

    @patch("context.workflow.read_excel")
    def test_read_param_file(self, mock_read_excel):
        # Mock the DataFrame returned by read_excel
        data = {
            "Option": ["CDS_or_proteins", "input_genomes_path"],
            "Value": ["CDS", "/path/to/genomes"],
        }
        df_mock = DataFrame(data).set_index("Option")
        mock_read_excel.return_value = df_mock

        self.runner.read_param_file()

        self.assertEqual(
            self.runner.params,
            {"CDS_or_proteins": "CDS", "input_genomes_path": "/path/to/genomes"},
        )

    @patch("context.workflow.Path")
    def test_validate_parameters(self, mock_path):
        self.runner.params = {
            "CDS_or_proteins": "CDS",
            "input_genomes_path": "/path/to/genomes",
            "input_cds_or_proteomes_path": "/path/to/proteomes",
            "input_genomes_or_proteomes_list": "/path/to/list",
        }

        mock_path.return_value = mock_path
        with self.assertRaises(RuntimeError):
            self.runner.validate_parameters()

        # Retry with complete parameters
        self.runner.params = {
            "CDS_or_proteins": "CDS",
            "input_genomes_path": "/path/to/genomes",
            "input_cds_or_proteomes_path": "/path/to/proteomes",
            "input_genomes_or_proteomes_list": "/path/to/list",
            "Calculate_orthologues": "Yes",
            "Calculate_cores": "Yes",
            "Calculate_fingerprints": "Yes",
            "Calculate_entire_phylogenomic_tree": "Fasttree",
            "Calculate_group_representatives_phylogenomic_tree": "IQTree",
            "Calculate_EGGNOG_on_Group_representatives": "Yes",
            "Calculate_EGGNOG_on_core/fingerprints": "Yes",
            "Calculate_SMBGCs_on_Group_representatives": "Yes",
            "Calculate_SMBGCs_on_entire_set": "Yes",
            "Calculate_CAZymes_on_Group_representatives": "Yes",
            "Calculate_CAZYmes_on_entire_set": "Yes",
            "Calculate_VFs_on_Group_representatives": "Yes",
            "Calculate_VFs_on_entire_set": "Yes",
            "Calculate_AMR_on_Group_representatives": "Yes",
            "Calculate_AMR_on_entire_set": "Yes",
        }
        self.runner.validate_parameters()

    @patch("context.workflow.read_csv")
    def test_load_input_genomes_or_proteomes_list(self, mock_read_csv):
        mock_df = DataFrame.from_dict(
            {
                "Proteome": ["Proteome1", "Proteome2"],
                "Group": ["Grp1", "Grp2"],
                "GroupRef": [1, 0],
                "WholeSet_Ref": [1, 0],
            },
        )
        mock_df = mock_df.set_index("Proteome")
        mock_read_csv.return_value = mock_df

        self.runner.params = {"input_genomes_or_proteomes_list": "/path/to/list"}
        self.runner.load_input_genomes_or_proteomes_list()

        # Test without WholeSetRef
        with self.assertRaises(ValueError):
            mock_df["WholeSet_Ref"] = [0, 0]
            mock_read_csv.return_value = mock_df
            self.runner.load_input_genomes_or_proteomes_list()

    def test_identify_tasks(self):
        self.runner.params = {
            "Calculate_orthologues": "Yes",
            "Calculate_cores": "No",
            "Calculate_fingerprints": "Yes",
            "Calculate_entire_phylogenomic_tree": "Fasttree",
        }

        self.runner.identify_tasks()

        expected_tasks = {
            "Calculate_orthologues": True,
            "Calculate_cores": False,
            "Calculate_fingerprints": True,
            "Calculate_entire_phylogenomic_tree": True,
            "Calculate_group_representatives_phylogenomic_tree": False,
            "Calculate_EGGNOG_on_Group_representatives": False,
            "Calculate_EGGNOG_on_core/fingerprints": False,
            "Calculate_SMBGCs_on_Group_representatives": False,
            "Calculate_SMBGCs_on_entire_set": False,
            "Calculate_CAZymes_on_Group_representatives": False,
            "Calculate_CAZYmes_on_entire_set": False,
            "Calculate_VFs_on_Group_representatives": False,
            "Calculate_VFs_on_entire_set": False,
            "Calculate_AMR_on_Group_representatives": False,
            "Calculate_AMR_on_entire_set": False,
        }

        self.assertEqual(self.runner.tasks_to_execute, expected_tasks)

    def test_validate_co_dependent_tasks(self):
        self.runner.tasks_to_execute = {
            "Calculate_orthologues": False,
            "Calculate_cores": True,
            "Calculate_fingerprints": True,
        }

        with self.assertRaises(ValueError) as context:
            self.runner.validate_co_dependent_tasks()
            self.assertTrue(
                "Cannot run 'Calculate_cores' if 'Calculate_orthologues'"
                in str(context.exception)
            )

    @patch("pathlib.Path.glob")
    def test_gather_cds_or_protein_fasta_files(self, mock_glob):
        self.runner.params = {"input_cds_or_proteomes_path": Path("/path/to/proteomes")}
        mock_file = MagicMock()
        mock_file.stem = "mock_gcf"
        mock_glob.return_value = [mock_file]

        mock_df = MagicMock()
        mock_df.index = ["mock_gcf"]
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.gather_cds_or_protein_fasta_files()

        self.assertEqual(self.runner.cds_or_protein_fasta_files, [mock_file])

    @patch("pathlib.Path.glob")
    def test_gather_genomic_fasta_files(self, mock_glob):
        self.runner.params = {"input_genomes_path": Path("/path/to/genomes")}
        mock_file = MagicMock()
        mock_file.stem = "mock_gcf"
        mock_glob.return_value = [mock_file]

        mock_df = MagicMock()
        mock_df.index = ["mock_gcf"]
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.gather_genomic_fasta_files()

        self.assertEqual(self.runner.genomic_fasta_files, [mock_file])

    def test_get_per_group_representatives_files(self):
        mock_df = DataFrame.from_dict(
            {
                "Proteome": ["Proteome1", "Proteome2"],
                "Group": ["Grp1", "Grp2"],
                "Group_Ref": [1, 0],
                "WholeSet_Ref": [1, 0],
            },
        )
        mock_df = mock_df.set_index("Proteome")
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.cds_or_protein_fasta_files = [
            Path("Proteome1.fa"),
            Path("Proteome2.fa"),
            Path("Proteome3.fa"),
        ]
        self.runner.genomic_fasta_files = [
            Path("Proteome1.fna"),
            Path("Proteome2.fna"),
            Path("Proteome3.fna"),
        ]

        self.runner.get_per_group_representatives_files()

        self.assertEqual(
            self.runner.cds_or_protein_fasta_files_representatives,
            [Path("Proteome1.fa")],
        )
        self.assertEqual(
            self.runner.genomic_fasta_files_representatives,
            [Path("Proteome1.fna")],
        )

    def test_create_orthologues_ref_list(self):
        mock_df = DataFrame.from_dict(
            {
                "Proteome": ["Proteome1", "Proteome2"],
                "Group": ["Grp1", "Grp2"],
                "Group_Ref": [1, 0],
                "WholeSet_Ref": [1, 0],
            },
        )
        mock_df = mock_df.set_index("Proteome")
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.create_orthologues_ref_list()

        self.assertEqual(self.runner.orthologues_ref_list, ["Proteome1"])

    @patch("pathlib.Path.exists")
    def test_create_core_ref_list(self, mock_exists):
        mock_df = DataFrame.from_dict(
            {
                "Proteome": ["Proteome1", "Proteome2"],
                "Group": ["Grp1", "Grp2"],
                "Group_Ref": [1, 0],
                "WholeSet_Ref": [1, 0],
            },
        )
        mock_df = mock_df.set_index("Proteome")
        self.runner.orthologues_ref_list = ["Proteome1"]
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.out_dir = Path("/mock/out/dir")
        self.runner.create_core_ref_list()

        expected_core_ref_list = [
            self.runner.out_dir / "Orthologues" / "Proteome1" / "OGMatrix.csv"
        ]
        expected_whole_set_ref_list = [
            self.runner.out_dir / "Orthologues" / "Proteome1" / "OGMatrix.csv"
        ]
        self.assertEqual(self.runner.core_ref_list, expected_core_ref_list)
        self.assertEqual(
            self.runner.core_whole_set_ref_list, expected_whole_set_ref_list
        )

    @patch("context.workflow.create_temporary_file")
    @patch("pathlib.Path.write_text")
    def test_write_species_file_for_core(self, mock_write_text, mock_create_temp_file):
        mock_create_temp_file.return_value = Path("/mock/temp/file")
        mock_df = MagicMock()
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.write_species_file_for_core()

        mock_create_temp_file.assert_called_once_with()
        mock_df[["Group"]].to_excel.assert_called_once_with(
            Path("/mock/temp/file.xlsx")
        )

    def test_get_phylogenomic_orthologue_matrix_file(self):
        mock_df = DataFrame.from_dict(
            {
                "Proteome": ["mock_ref", "Proteome2"],
                "Group": ["Grp1", "Grp2"],
                "Group_Ref": [1, 0],
                "WholeSet_Ref": [1, 0],
            },
        )
        mock_df = mock_df.set_index("Proteome")
        self.runner.input_genomes_or_proteomes_list_df = mock_df

        self.runner.out_dir = Path("/mock/out/dir")

        self.runner.get_phylogenomic_orthologue_matrix_file()

        self.assertEqual(
            self.runner.phylogenomic_og_matrix,
            Path("/mock/out/dir/Orthologues/mock_ref/OGmatrix.csv"),
        )

    @patch("context.workflow.Path.exists")
    def test_create_eggnog_core_protein_files_list(self, mock_path_exist):
        mock_path_exist.return_value = True
        self.runner.orthologues_ref_list = ["mock_ref1", "mock_ref2"]

        self.runner.out_dir = Path("/mock/out/dir")

        self.runner.create_eggnog_core_protein_files_list()

        expected_list = [
            self.runner.out_dir / "Core_and_fingerprints" / "mock_ref1_core.xlsx",
            self.runner.out_dir
            / "Core_and_fingerprints"
            / "mock_ref1_species_core.xlsx",
            self.runner.out_dir / "Core_and_fingerprints" / "mock_ref2_core.xlsx",
            self.runner.out_dir
            / "Core_and_fingerprints"
            / "mock_ref2_species_core.xlsx",
        ]

        self.assertEqual(self.runner.emapper_core_protein_files_reflist, expected_list)
