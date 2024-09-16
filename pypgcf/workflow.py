from pathlib import Path
from pandas import read_excel, read_csv
from typing import Dict, Any
import logging
from pypgcf.utils import create_temporary_file

OPTIONS = {
    "CDS_or_proteins": ["CDS", "Proteins"],
    "input_genomes_path": [Path],
    "input_cds_or_proteomes_path": [Path],
    "input_genomes_or_proteomes_list": [Path],
    "Calculate_orthologues": ["Yes", "No"],
    "Calculate_cores": ["Yes", "No"],
    "Calculate_fingerprints": ["Yes", "No"],
    "Calculate_entire_phylogenomic_tree": ["NJ", "IQTree", "Fasttree", "No"],
    "Calculate_group_representatives_phylogenomic_tree": [
        # "NJ",
        "IQTree",
        "Fasttree",
        "No",
    ],
    "Calculate_EGGNOG_on_Group_representatives": ["Yes", "No"],
    "Calculate_EGGNOG_on_core/fingerprints": ["Yes", "No"],
    "Calculate_SMBGCs_on_Group_representatives": ["Yes", "No"],
    "Calculate_SMBGCs_on_entire_set": ["Yes", "No"],
    "Calculate_CAZymes_on_Group_representatives": ["Yes", "No"],
    "Calculate_CAZYmes_on_entire_set": ["Yes", "No"],
    "Calculate_VFs_on_Group_representatives": ["Yes", "No"],
    "Calculate_VFs_on_entire_set": ["Yes", "No"],
    "Calculate_AMR_on_Group_representatives": ["Yes", "No"],
    "Calculate_AMR_on_entire_set": ["Yes", "No"],
}

NECESSARY_OPTIONS = [
    "CDS_or_proteins",
    "input_genomes_path",
    "input_cds_or_proteomes_path",
    "input_genomes_or_proteomes_list",
]

TASK_MAP = {  # Map tasks to module
    "Calculate_orthologues": "orthologues",
    "Calculate_cores": "core",
    "Calculate_fingerprints": "core",
    "Calculate_entire_phylogenomic_tree": "phylogenomic",
    "Calculate_group_representatives_phylogenomic_tree": "phylogenomic",
    "Calculate_EGGNOG_on_Group_representatives": "eggnog",
    "Calculate_SMBGCs_on_Group_representatives": "smbgc",
    "Calculate_SMBGCs_on_entire_set": "smbgc",
    "Calculate_CAZymes_on_Group_representatives": "cazy",
    "Calculate_CAZYmes_on_entire_set": "cazy",
    "Calculate_VFs_on_Group_representatives": "virulence",
    "Calculate_VFs_on_entire_set": "virulence",
    "Calculate_AMR_on_Group_representatives": "amr",
    "Calculate_AMR_on_entire_set": "amr",
}


class WorkflowRunner:
    def __init__(self, *, out_dir: Path, param_file: Path, debug: bool = False):
        self.out_dir = out_dir
        self.param_file = param_file
        self.debug = debug
        self.params: Dict[str, Any] = {}
        self.tasks: Dict[str, bool] = {}
        self.cds_or_protein_fasta_files = []
        self.genomic_fasta_files = []

        if self.debug:
            logging.debug(
                f"Initialized WorkflowRunner with param_file: {param_file} and debug: {debug}"
            )

    def read_param_file(self) -> None:
        if self.debug:
            logging.debug(f"Parsing Excel file: {self.param_file}")

        # Read the Excel file
        df = read_excel(self.param_file, index_col=0)

        if self.debug:
            logging.debug(
                f"Excel file read successfully with options: {df.index.tolist()}"
            )

        # Create the hashmap
        self.params = {key: value for key, value in df.itertuples(index=True)}

        if self.debug:
            logging.debug(f"Parameters hashmap created: {self.params}")

        return None

    def validate_parameters(self) -> None:
        if self.debug:
            logging.debug("Validating parameters")
        for input, expected_val in OPTIONS.items():
            val = self.params.get(input)
            if val is None:
                logging.error(f"{input} was ommitted from parameters file")
                raise RuntimeError(f"{input} was ommitted from parameters file")
            if "input_" in input:
                self.params[input] = Path(val)
                continue
            if val in expected_val:
                continue
        return None

    def load_input_genomes_or_proteomes_list(self) -> None:
        fin = self.params.get("input_genomes_or_proteomes_list")
        self.input_genomes_or_proteomes_list_df = read_csv(fin, sep="\t", index_col=0)

        mask = self.input_genomes_or_proteomes_list_df["WholeSet_Ref"] == 1
        whole_set_ref_list = self.input_genomes_or_proteomes_list_df.loc[
            mask
        ].index.tolist()
        if len(whole_set_ref_list) == 0:
            raise ValueError(f"WholeSet_Ref must be specified in {fin}")

        return None

    def identify_tasks(self) -> None:
        """
        Identify which tasks will be executed by the class instance based on the parameters
        """
        tasks = {
            "Calculate_orthologues": False,
            "Calculate_cores": False,
            "Calculate_fingerprints": False,
            "Calculate_entire_phylogenomic_tree": False,
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
        for task in TASK_MAP.keys():
            task_val = self.params.get(task, "No")
            if task_val != "No":
                tasks[task] = True
        self.tasks_to_execute = tasks

    def validate_co_dependent_tasks(self) -> None:
        """
        Certain tasks depend on the output of others in order to run.
        Identify if these tasks and those they depend upon are going to run.
        Such examples include the core and fingerprints that need orthologues to be run first
        """
        co_dependent_tasks = {
            "Calculate_cores": ["Calculate_orthologues"],
            "Calculate_fingerprints": ["Calculate_orthologues"],
            "Calculate_entire_phylogenomic_tree": ["Calculate_orthologues"],
            "Calculate_group_representatives_phylogenomic_tree": [
                "Calculate_orthologues"
            ],
            "Calculate_EGGNOG_on_core/fingerprints": [
                "Calculate_EGGNOG_on_Group_representatives",
                "Calculate_cores",
                "Calculate_fingerprints",
            ],
        }
        for dependent, dependee_list in co_dependent_tasks.items():
            dependent_status = self.tasks_to_execute.get(dependent)
            if dependent_status is False:
                continue
            dependee_status = set(
                [self.tasks_to_execute.get(dependee) for dependee in dependee_list]
            )
            if False in dependee_status or None in dependee_status:
                raise ValueError(
                    f"Cannot run '{dependent}' if '{','.join(dependee_list)}' are not set to 'Yes'"
                )

    def gather_cds_or_protein_fasta_files(self) -> None:
        """
        Gather all fasta files from the input directory and keep those found in the genomes/proteome list
        """
        fasta_dir = self.params.get("input_cds_or_proteomes_path")
        files = []
        wanted_gcfs = {i: True for i in self.input_genomes_or_proteomes_list_df.index}
        for f in fasta_dir.glob("*"):
            gcf = f.stem
            if gcf not in wanted_gcfs:
                continue
            files.append(f)
        self.cds_or_protein_fasta_files = files

        return None

    def gather_genomic_fasta_files(self) -> None:
        """
        Gather all fasta files from the input directory and keep those found in the genomes/proteome list
        """
        fasta_dir = self.params.get("input_genomes_path")
        files = []
        wanted_gcfs = {i: True for i in self.input_genomes_or_proteomes_list_df.index}
        for f in fasta_dir.glob("*"):
            gcf = f.stem
            if gcf not in wanted_gcfs:
                continue
            files.append(f)
        self.genomic_fasta_files = files

        return None

    def get_per_group_representatives_files(self) -> None:
        """
        Create sublists of genomic and cds/proteome fasta files based on the provided genomes/proteomes list
        """
        mask = self.input_genomes_or_proteomes_list_df["Group_Ref"] == 1
        per_group_representatives = set(
            self.input_genomes_or_proteomes_list_df.loc[mask].index.tolist()
        )
        cds_or_protein_fasta_files_representatives = []
        for file in self.cds_or_protein_fasta_files:
            stem = file.stem
            if stem not in per_group_representatives:
                continue
            cds_or_protein_fasta_files_representatives.append(file)

        genomic_fasta_files_representatives = []
        for file in self.genomic_fasta_files:
            stem = file.stem
            if stem not in per_group_representatives:
                continue
            genomic_fasta_files_representatives.append(file)

        self.cds_or_protein_fasta_files_representatives = (
            cds_or_protein_fasta_files_representatives
        )
        self.genomic_fasta_files_representatives = genomic_fasta_files_representatives
        return None

    def create_orthologues_ref_list(self) -> None:
        """
        Create the reference list need by the orthologues module.
        This list will be created based on the reference status of the organisms in the genomes/proteomes list.
        """
        mask = (self.input_genomes_or_proteomes_list_df["Group_Ref"] == 1) | (
            self.input_genomes_or_proteomes_list_df["WholeSet_Ref"] == 1
        )
        self.orthologues_ref_list = self.input_genomes_or_proteomes_list_df.loc[
            mask
        ].index.tolist()

        return None

    def create_core_ref_list(self) -> None:
        wholeset_ref = self.input_genomes_or_proteomes_list_df[
            self.input_genomes_or_proteomes_list_df["WholeSet_Ref"] == 1
        ].index.tolist()[0]
        ref_list = []

        wholeset_ref_file = None
        for ref in self.orthologues_ref_list:
            og_matrix_f = self.out_dir / "Orthologues" / ref / "OGMatrix.csv"
            ref_list.append(og_matrix_f)
            if ref == wholeset_ref:
                wholeset_ref_file = og_matrix_f

        if wholeset_ref_file is None:
            raise ValueError("OGMatrix based on WholeSet reference wasn't found")

        self.core_ref_list = ref_list
        self.core_whole_set_ref_list = [wholeset_ref_file]
        return None

    def write_species_file_for_core(self) -> None:
        tmp_species_file = create_temporary_file()
        tmp_species_file = tmp_species_file.with_suffix(".xlsx")
        self.input_genomes_or_proteomes_list_df[["Group"]].to_excel(tmp_species_file)
        self.core_species_file = tmp_species_file
        return None

    def get_phylogenomic_orthologue_matrix_file(self) -> None:
        """
        Get the orthologue matrix file that will be used by phylogenomic
        This file will be calculated using the WholeSet_Ref organism.
        """
        wholeset_ref = self.input_genomes_or_proteomes_list_df[
            self.input_genomes_or_proteomes_list_df["WholeSet_Ref"] == 1
        ].index.tolist()[0]
        file = self.out_dir / "Orthologues" / wholeset_ref / "OGmatrix.csv"
        self.phylogenomic_og_matrix = file
        return None

    def create_eggnog_core_protein_files_list(self) -> None:
        ref_list_core = []

        for ref in self.orthologues_ref_list:
            core = self.out_dir / "Core_and_fingerprints" / (ref + "_core.xlsx")
            species_core = (
                self.out_dir / "Core_and_fingerprints" / (ref + "_species_core.xlsx")
            )
            if core.exists():
                ref_list_core.append(core)
            if species_core.exists():
                ref_list_core.append(species_core)

        self.emapper_core_protein_files_reflist = ref_list_core

        return None
