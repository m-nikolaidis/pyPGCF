from pathlib import Path
import unittest
from pypgcf import core


class TestModule(unittest.TestCase):
    def test_split_genomes_into_groups(self):
        indir = Path("data")
        og_matrix = indir / "OGmatrix.csv"
        species_file = indir / "FastANI_species_clusters.xlsx"
        out_dir = Path("data")
        core_perc = 100.0
        instance = core.Core_identifier(og_matrix, out_dir, species_file, core_perc)
        instance.split_genomes_into_groups()
        self.assertEqual(len(instance.group_orgs), 1)
        self.assertNotIn("GCF_000009045.1", instance.group_orgs)
        self.assertIn("GCF_019704415.1", instance.group_orgs)
        self.assertEqual("GCF_000009045.1", instance.ref)
        self.assertEqual(len(instance.non_group_orgs), 1)

    def test_calculate_protein_presence(self):
        indir = Path("data")
        og_matrix = indir / "OGmatrix.csv"
        species_file = indir / "FastANI_species_clusters.xlsx"
        out_dir = Path("data")
        core_perc = 100
        instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
        instance.split_genomes_into_groups()
        presence_df = instance.calculate_protein_presence()
        num_core = len(
            presence_df[presence_df["Group orthologues %"] == 100].index.tolist()
        )
        self.assertEqual(num_core, 4073)
        specific_prot = "NP_387886.2"
        print(instance.orthology_df.loc[specific_prot])
        print(instance.group_orgs)
        specific_prot_presence = presence_df.loc[specific_prot, "Group orthologues %"]
        specific_prot_presence_outside = presence_df.loc[
            specific_prot, "Non group orthologues %"
        ]
        print(presence_df.loc[specific_prot])
        self.assertEqual(specific_prot_presence, 100)
        self.assertEqual(specific_prot_presence_outside, 0)

    def test_identify_core_proteins(self):
        indir = Path("data")
        og_matrix = indir / "OGmatrix.csv"
        species_file = indir / "FastANI_species_clusters.xlsx"
        out_dir = Path("data")
        core_perc = 100
        instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
        df, fout = instance.calculate_core()
        self.assertEqual(df["Core_100%"].sum(), 4073)
        self.assertEqual(df["Is fingerprint"].sum(), 1936)

    def test_split_genomes_into_groups_unequal_input(self):
        """
        Test if the function can handle species clusters with more genomes than the orthology matrix
        """
        pass

    # def test_split_genomes_into_groups_various_perc(self):
    #     indir = Path("data")
    #     og_matrix = indir / "OGmatrix.csv"
    #     species_file = indir / "FastANI_species_clusters.xlsx"
    #     out_dir = Path("data")
    #     core_percents = [100,98,95,90]
    #     for core_perc in core_percents:
    #         instance = Core_identifier(og_matrix, species_file, core_perc, out_dir)
    #         instance.calculate_core()
