"""
Calculate the core proteins and fingerprints based on a given reference
"""

from pathlib import Path
from typing import Union

from pandas import DataFrame, read_csv
from pypgcf.utils import dict_to_dataframe


class Core_identifier:
    def __init__(
        self,
        *,
        orthology_fin: Path,
        out_dir: Path,
        species_df: Union[None, DataFrame],
        core_perc: float,
    ):
        """ """
        self.core_perc = core_perc
        self.out_dir = out_dir / "Core_and_fingerprints"
        self.orthology_fin = orthology_fin
        # Split the genomes into groups
        self.genus = True
        self.species_df = species_df
        if species_df is not None:
            self.genus = False

    def _turn_orthology_df_to_binary(self):
        return self.orthology_df.map(lambda x: False if x == "X" else True).astype(bool)

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)

    def split_genomes_into_groups(self):
        """
        Split the genomes into groups based on the species cluster (cluster_col)
        If genus is True, then the complete table (taxon) is used
        """
        total_genomes = set(self.orthology_df.columns.tolist())
        if self.genus:
            self.group_orgs = list(total_genomes)
            self.non_group_orgs = []
            return
        species_col = self.species_df.columns[0]
        ref_cluster = self.species_df.loc[self.ref, species_col]
        ref_sp_genomes = set(
            self.species_df[self.species_df[species_col] == ref_cluster].index.tolist()
        )
        non_group_orgs = set(
            self.species_df[self.species_df[species_col] != ref_cluster].index.tolist()
        )
        group_genomes = ref_sp_genomes.intersection(total_genomes)
        non_group_genomes = non_group_orgs.intersection(total_genomes)
        self.group_orgs = list(group_genomes)
        self.non_group_orgs = list(non_group_genomes)

    def calculate_protein_presence(self) -> DataFrame:
        tmpdf = self.orthology_df.copy().drop(self.orthology_df.columns, axis=1)
        tmpdf["Group orthologues"] = self.orthology_df[self.group_orgs].sum(axis=1)
        tmpdf["Group orthologues"] = tmpdf["Group orthologues"] + 1
        tmpdf["Group orthologues %"] = round(
            (tmpdf["Group orthologues"] / (len(self.group_orgs) + 1)) * 100, 2
        )
        tmpdf["Non group orthologues"] = self.orthology_df[self.non_group_orgs].sum(
            axis=1
        )
        tmpdf["Non group orthologues %"] = round(
            tmpdf["Non group orthologues"] / len(self.non_group_orgs) * 100, 2
        )
        return tmpdf

    def identify_core_proteins(self, df: DataFrame) -> dict:
        """"""
        core_proteins = df[df["Group orthologues %"] >= self.core_perc].index.tolist()
        fingerprints = df[
            (df["Group orthologues %"] == 100) & (df["Non group orthologues %"] == 0)
        ].index.tolist()
        data = {}
        core_perc_col = f"Core_{self.core_perc}%"
        for protein in df.index.tolist():
            data[protein] = {core_perc_col: 0, "Is fingerprint": 0}
            if protein in fingerprints:
                data[protein]["Is fingerprint"] = 1
                data[protein][core_perc_col] = 1
                continue
            if protein in core_proteins:
                data[protein][core_perc_col] = 1
        return data

    def load_orthology_matrix(self):
        if self.orthology_fin is not None:
            self.orthology_df = read_csv(self.orthology_fin, sep="\t", index_col=0)
            self.orthology_df = self._turn_orthology_df_to_binary()
            self.ref = self.orthology_df.index.name

    def calculate_core(self):
        self.setup_directories()
        self.load_orthology_matrix()
        self.split_genomes_into_groups()
        # Calculate the presence of each protein
        final_df = self.calculate_protein_presence()
        # Identify the core proteins
        core_protein_data = self.identify_core_proteins(final_df)
        core_protein_data_df = dict_to_dataframe(core_protein_data)
        core_protein_data_df.index.name = self.ref
        if self.genus:
            fout = self.out_dir / f"{self.ref}_core.xlsx"
            core_protein_data_df = core_protein_data_df.drop("Is fingerprint", axis=1)
        else:
            fout = self.out_dir / f"{self.ref}_species_core.xlsx"
        core_protein_data_df.to_excel(fout)
