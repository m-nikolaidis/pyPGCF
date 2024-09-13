"""
Module to query reference proteome to eggnog database
Used to add functional annotation to the proteins in the core genes excel files
"""

from datetime import datetime
from pathlib import Path
from typing import List, Union
import logging

from pandas import read_excel, read_csv, DataFrame, Series, ExcelWriter
from Bio.SeqIO import parse
from scipy.stats import hypergeom

from pypgcf.utils import multiprocess_dispatch, dict_to_dataframe

COG_CATEGORY_ANNOTATION = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: Prophages, transposons",
    "R": "General function, prediction only",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}

EGGNOG_HEADERS = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]

logging.basicConfig(level=logging.DEBUG)


class EggNOGRunner:
    def __init__(
        self,
        *,
        fasta_files_list: list[Path],
        out_dir: Path,
        cores: int,
        pident: float,
        qcov: float,
        scov: float,
        input_type: bool,
        database_dir: Path,
        debug: bool = False,
    ):
        self.fasta_files = fasta_files_list
        self.out_dir = out_dir / "eggNOG"
        self.emapper_resdir = self.out_dir / "eggnog_mapper"
        self.cores = cores
        self.pident = pident
        self.qcov = qcov
        self.scov = scov
        self.nucl = True
        self.database_dir = database_dir
        if input_type == "prot":
            self.nucl = False
        self.debug = debug

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)
        self.emapper_resdir.mkdir(exist_ok=True, parents=True)

    def create_eggnog_cmd(self):
        cmds = []
        for fasta_file in self.fasta_files:
            fout = self.emapper_resdir / (fasta_file.stem + "_eggNOG_results.csv")
            eggnog_cmd = " ".join(
                [
                    "emapper.py",
                    f"--cpu {self.cores}",
                    f"--pident {self.pident}",
                    f"--query_cov {self.qcov}",
                    f"--subject_cov {self.scov}",
                    f"-o {fout}",
                    f"-i {fasta_file}",
                    f"--data_dir {self.database_dir}",
                ]
            )
            if self.nucl:
                eggnog_cmd += " --itype CDS --translate"
            if not self.debug:
                eggnog_cmd += " > /dev/null"
            cmds.append(eggnog_cmd)
        self.eggnog_cmds = cmds

    def clean_unessecary_output(self):
        files_to_remove = []
        for fasta_file in self.fasta_files:
            tmpf = self.emapper_resdir / (fasta_file.stem + "_eggNOG_results.csv")
            files_to_remove.append(tmpf.with_suffix(".csv.emapper.hits"))
            files_to_remove.append(tmpf.with_suffix(".csv.emapper.seed_orthologs"))
            file_to_rename = tmpf.with_suffix(".csv.emapper.annotations")
            file_to_rename.rename(tmpf)

    def execute_eggnog_mapper(self):
        logging.info(
            f"Executing eggNOG-mapper: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.setup_directories()
        self.create_eggnog_cmd()
        _ = multiprocess_dispatch(
            "system",
            self.eggnog_cmds,
            num_procs=1,
            show_progress=True,
            description="Running eggNOG-mapper on representatives",
        )
        self.clean_unessecary_output()


class EggNOGParser:
    def __init__(
        self,
        *,
        fasta_files_list: list[Path],
        core_protein_table_file_list: Union[List[Path], None],
        out_dir: Path,
        debug: bool,
    ):
        self.fasta_files = fasta_files_list
        self.out_dir = out_dir / "eggNOG"

        if core_protein_table_file_list is None:
            self.core_protein_table_file_list = []
        else:
            self.core_protein_table_file_list = core_protein_table_file_list

        self.debug = debug
        self.proteome_cog_data = {}

        if self.debug:
            logging.info(
                f"Initialized EggNOGParser with {len(fasta_files_list)} fasta files."
            )

    # TODO: Need to figure out how to pre-gather the emapper results csv files and match the fasta files
    # then output a warning if the results for a certain fasta file are not found currently the error would stop the execution of the programm.
    # Should I leave it be?

    def load_all_reference_proteins(self, fasta_file: Path) -> dict:
        """
        Load reference proteins from a FASTA file, initializing their categories to ["S"].
        This way all proteins have a default value which will be updated later if the protein is reported in the eggNOG results.
        """
        record = parse(str(fasta_file), "fasta")
        return {protein.id: ["S"] for protein in record}

    def turn_category_to_list(self, category: str) -> list:
        """
        Convert a COG category string to a list. If "-", default to ["S"].
        """
        if category == "-":
            return ["S"]
        return list(category)

    def get_assigned_cog_categories(self, eggnog_results_df: DataFrame) -> dict:
        """
        Turn assigned COG categories a eggNOG results DataFrame to list and return a protein-based key dictionary
        """
        column = "COG_category"
        return eggnog_results_df[column].map(self.turn_category_to_list).to_dict()

    def parse_eggnog_raw_results(self):
        """
        Parse eggNOG raw results for each reference fasta file.
        """
        # Load reference proteins
        for fasta_file in self.fasta_files:
            reference = fasta_file.stem
            if self.debug:
                logging.info(f"Parsing eggNOG results for {reference}")
            ref_data = self.load_all_reference_proteins(fasta_file)

            # Load eggnog annotation results for current organism
            eggnog_results_f = (
                self.out_dir / "eggnog_mapper" / (reference + "_eggNOG_results.csv")
            )
            eggnog_results_df = self._load_eggnogmapper_csv(eggnog_results_f)

            # Get annotation dictionary
            eggnog_annotation = self.get_assigned_cog_categories(eggnog_results_df)

            # Update dictionary data
            ref_data.update(eggnog_annotation)
            self.proteome_cog_data[reference] = ref_data

    def _load_eggnogmapper_csv(self, file: Path) -> DataFrame:
        if not file.exists():
            logging.error(f"File {file} not found.")
            raise FileNotFoundError

        return read_csv(
            file,
            sep="\t",
            index_col=0,
            comment="#",
            names=EGGNOG_HEADERS,
        )

    def calculate_cog_category_counts_per_proteome(self, cog_data_dict: dict) -> dict:
        """
        Calculate the number of proteins in each COG category per proteome
        cog_data_dict example:
        {"Proteome1": {"p1": ["S"]}, "Proteome2": {"p1":["K"], "p2": ["K", "L"]}}
        """
        # Initialize counts
        counts = {
            p: {c: 0 for c in COG_CATEGORY_ANNOTATION.keys()}
            for p in cog_data_dict.keys()
        }

        for proteome, cog_data in cog_data_dict.items():
            counts[proteome]["Total"] = len(cog_data.keys())
            for cog_categories_list in cog_data.values():
                for category in cog_categories_list:
                    counts[proteome][category] += 1
        return counts

    def write_counts_dictionary_to_excel(self, counts_dict: dict, fout: Path) -> None:
        df = dict_to_dataframe(counts_dict)
        df_perc = df.div(df["Total"], axis=0) * 100
        with ExcelWriter(fout) as writer:
            df.to_excel(writer, sheet_name="counts")
            df_perc.to_excel(writer, sheet_name="percent")
        return None

    def read_core_protein_tables(self):
        """
        core_protein_data structure example:
        {
            "Proteome1": {
                 {"p1": ["L"], "p2": ["KT"]},
            },
            "Proteome2": {
                 {"p1": ["S"], "p2": ["KT"]},
            },
        }
        """
        core_protein_data = {}
        fingerprint_protein_data = {}
        for table_f in self.core_protein_table_file_list:
            core_protein_table = read_excel(table_f, index_col=0)
            ref = str(core_protein_table.index.name)
            core_proteins, fingerprints = self._get_protein_subsets_from_core_prot_df(
                core_protein_table
            )
            if ref not in self.proteome_cog_data:
                logging.error(f"emapper results for {ref} weren't parsed")
                core_protein_data[ref] = {}
                fingerprint_protein_data[ref] = {}
                continue

            # Initialize the dictionaries
            core_protein_data[ref] = {
                p: self.proteome_cog_data[ref].get(p) for p in core_proteins
            }
            fingerprint_protein_data[ref] = {
                p: self.proteome_cog_data[ref].get(p) for p in fingerprints
            }

        self.core_protein_cog_counts = self.calculate_cog_category_counts_per_proteome(
            core_protein_data
        )
        self.fingerprint_protein_cog_counts = (
            self.calculate_cog_category_counts_per_proteome(fingerprint_protein_data)
        )

    def _get_protein_subsets_from_core_prot_df(self, df: DataFrame) -> dict:
        """
        Get the core proteins and fingerprints (if availale) from core protein dataframe
        """
        core_proteins = []
        fingerprints = []
        for col in df.columns:
            mask = df[col] == 1
            if "fingerprint" in col:
                fingerprints = df.loc[mask].index.tolist()
            else:
                core_proteins = df.loc[mask].index.tolist()
        return core_proteins, fingerprints

    def perform_hypergeom_test(
        self,
        *,
        population_size: int,
        population_successes: int,
        sample_size: int,
        sample_successes: int,
    ) -> tuple:
        expected = (sample_size * population_successes) / population_size
        if sample_successes > sample_size:
            raise ValueError(
                "Hypergeometric test sample size cannot be less than sample sample_successes"
            )
        if sample_size > population_size:
            raise ValueError(
                "Hypergeometric test sample size cannot be more than population size"
            )

        # p-value
        pval = hypergeom.sf(
            sample_successes - 1, population_size, population_successes, sample_size
        )
        if sample_successes < expected or sample_successes == 0:
            pval = 1 - hypergeom.sf(
                sample_successes, population_size, population_successes, sample_size
            )
        # fold change
        ratio_population = population_successes / population_size
        ratio_sample = sample_successes / sample_size
        if ratio_population == 0 or ratio_sample == 0:
            return pval, 0
        fold_change = ratio_sample / ratio_population
        return pval, fold_change

    def compare_core_and_fingerprints_against_background(self) -> dict:
        hypergeom_dfs = []
        subsets = [self.core_protein_cog_counts, self.fingerprint_protein_cog_counts]
        for subset in subsets:
            if len(subset) == 0:  # Empty subset
                continue
            data = {}
            for proteome, protein_cog_counts in subset.items():
                if len(protein_cog_counts) == 0:
                    # No results were found by the read_core_protein_tables(self) function
                    continue
                data[proteome] = {}
                population_size = self.proteome_cog_counts[proteome]["Total"]
                sample_size = subset[proteome]["Total"]
                for category in protein_cog_counts.keys():
                    if category == "Total":
                        continue
                    population_successes = self.proteome_cog_counts[proteome][category]
                    sample_successes = subset[proteome][category]
                    pvalue, fold_change = self.perform_hypergeom_test(
                        population_size=population_size,
                        population_successes=population_successes,
                        sample_size=sample_size,
                        sample_successes=sample_successes,
                    )
                    data[proteome][f"{category}_pvalue"] = pvalue
                    data[proteome][f"{category}_fold_change"] = fold_change
            hypergeom_df = dict_to_dataframe(data)
            hypergeom_dfs.append(hypergeom_df)
        if len(hypergeom_dfs) == 1:
            return {"core": hypergeom_dfs[0]}
        return {"core": hypergeom_dfs[0], "fingerprint": hypergeom_dfs[1]}

    def highlight_pvalue_on_output(self, row: Series) -> List:
        over_rep_colour = "background-color:green"
        under_rep_colour = "background-color:red"
        return_value = ["" for _ in row.index]
        pvalue_idx = [idx for idx in row.index if "pvalue" in idx]
        fold_change_idx = [idx for idx in row.index if "fold_change" in idx]
        for pvalue_col, fold_change_col in zip(pvalue_idx, fold_change_idx):
            pvalue = row[pvalue_col]
            fold_change = row[fold_change_col]
            pvalue_col_idx = row.index.get_loc(pvalue_col)
            fold_change_col_idx = row.index.get_loc(fold_change_col)
            if type(pvalue) is str:
                return return_value
            if pvalue <= 0.05:
                if fold_change > 1:
                    return_value[pvalue_col_idx] = over_rep_colour
                    return_value[fold_change_col_idx] = over_rep_colour
                if fold_change < 1:
                    return_value[pvalue_col_idx] = under_rep_colour
                    return_value[fold_change_col_idx] = under_rep_colour
        return return_value

    def write_hypergeometric_dfs(self, fout: Path, hypergeom_results: dict) -> None:
        excel_writer = ExcelWriter(fout)
        for protein_set, df in hypergeom_results.items():
            df = df.style.apply(self.highlight_pvalue_on_output, axis=1)
            df.to_excel(excel_writer, sheet_name=protein_set)
        excel_writer.close()
        return None

    # Wrapper functions
    def gather_eggnog_proteome_results(self):
        logging.info(
            f"Parsing eggNOG results: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.parse_eggnog_raw_results()
        self.proteome_cog_counts = self.calculate_cog_category_counts_per_proteome(
            self.proteome_cog_data
        )
        fout = self.out_dir / "COG_categories_of_proteomes.xlsx"
        self.write_counts_dictionary_to_excel(self.proteome_cog_counts, fout)

    def compare_proteome_to_core_and_fps(self):
        logging.info(
            f"Parsing eggNOG results: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.parse_eggnog_raw_results()
        self.proteome_cog_counts = self.calculate_cog_category_counts_per_proteome(
            self.proteome_cog_data
        )
        self.read_core_protein_tables()
        hypergeom_results = self.compare_core_and_fingerprints_against_background()
        fout = self.out_dir / "eggNOG_hypergeometric.xlsx"
        self.write_hypergeometric_dfs(fout, hypergeom_results)

    # #
    # # def compare_proteome_to_core_and_fps(self):
    # #     print(f"Parsing results: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    # #     total_data = {"core": [], "fingerprint": []}
    # #     for core_protein_table_f in self.core_protein_table_file_list:
    # #         self.core_protein_table = read_excel(core_protein_table_f, index_col=0)
    # #         self.ref = str(self.core_protein_table.index.name)
    # #         self.parse_eggnog_raw_results()
    # #         self.get_protein_subsets()
    # #         self.calculate_cog_categories_for_proteins()
    # #         tmp_data = self.compare_sets_with_hypergeometric_test()
    # #         total_data["core"].append(tmp_data["core"])
    # #         if "fingerprint" in tmp_data:
    # #             total_data["fingerprint"].append(tmp_data["fingerprint"])
    # #     self.hypergeometric_dfs = [
    # #         concat(total_data["core"]),
    # #         concat(total_data["fingerprint"]),
    # #     ]
    # #     if len(total_data["fingerprint"]) != 0:
    # #         self.hypergeometric_dfs.append(concat(total_data["fingerprint"]))
    # #     print(f"Writing results: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    # #     self.write_hypergeometric_dfs()
    # #     print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")

    # def calculate_cog_categories_for_proteins(self):
    #     # Initialize data
    #     cog_categories = {}
    #     for category in COG_CATEGORY_ANNOTATION:
    #         cog_categories[category] = {"proteome": 0}
    #         for protein_subset, proteins in self.protein_subsets.items():
    #             cog_categories[category][protein_subset] = 0
    #
    #     # Add proteome data
    #     for protein, protein_categories in self.proteome_cog.items():
    #         for category in protein_categories:
    #             cog_categories[category]["proteome"] += 1
    #
    #     # Add data for each subset (core, fingerprints)
    #     for protein_subset, proteins in self.protein_subsets.items():
    #         for protein in proteins:
    #             protein_categories = self.proteome_cog[
    #                 protein
    #             ]  # WIll always exist as key
    #             for category in protein_categories:
    #                 cog_categories[category][protein_subset] += 1
    #
    #     self.cog_categories_data = cog_categories
    #
    #

    # def compare_sets_with_hypergeometric_test(self) -> dict:
    #     background_set = "proteome"
    #     population_size = len(
    #         self.proteome_cog.keys()
    #     )  # Number of proteins in proteome
    #     dfs = []
    #     for protein_set, protein_list in self.protein_subsets.items():
    #         data = {self.ref: {}}
    #         sample_size = len(protein_list)
    #         if sample_size == 0:
    #             print(f"No proteins are part of the {protein_set} set. Skipping...")
    #             continue
    #         for category in self.cog_categories_data:
    #             population_successes = self.cog_categories_data[category][
    #                 background_set
    #             ]
    #             sample_successes = self.cog_categories_data[category][protein_set]
    #             pvalue, fold_change = self.perform_hypergeom_test(
    #                 population_size, population_successes, sample_size, sample_successes
    #             )
    #             data[self.ref][f"{category}_pvalue"] = pvalue
    #             data[self.ref][f"{category}_fold_change"] = fold_change
    #             # data[self.ref][category] = {"Annotation": self.category_annotation.get(category, "X"), "p-value": pvalue, "Fold_change": fold_change}
    #         df = DataFrame.from_dict(data, orient="index")
    #         dfs.append(df)
    #     if len(dfs) == 1:
    #         return {"core": dfs[0]}
    #     return {"core": dfs[0], "fingerprint": dfs[1]}
