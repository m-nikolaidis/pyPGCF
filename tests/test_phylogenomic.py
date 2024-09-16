from pathlib import Path
from Bio import SeqIO
import unittest

from context import phylogenomic


class TestModule(unittest.TestCase):
    def _cleanup_fasta(self):
        og_fasta_dir = Path("../data") / "Phylogenomic_tree" / "OGs_fasta"
        for fasta_file in og_fasta_dir.glob("*"):
            fasta_file.unlink()  # For cleaning
        og_fasta_dir.rmdir()

    def test_create_og_fasta(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_f,
            cores,
            fasta_dir,
            out_dir,
            no_keep_fasta,
            tree_model,
            debug=True,
        )

        phylo.setup_directories()
        phylo.create_og_fasta()
        og_fasta_dir = Path("../data") / "Phylogenomic_tree" / "OGs_fasta"
        fasta_files = list(og_fasta_dir.glob("*"))
        self.assertEqual(len(fasta_files), 2137)
        for fasta_file in fasta_files:
            fasta_file_parser = SeqIO.parse(str(fasta_file), "fasta")
            num_records = 0
            for rec in fasta_file_parser:
                self.assertGreater(len(rec.seq), 0)
                num_records += 1
            self.assertEqual(num_records, 3)
        # self._cleanup_fasta()

    def test_align_og_fasta(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_f,
            cores,
            fasta_dir,
            out_dir,
            no_keep_fasta,
            tree_model,
            debug=True,
        )
        phylo.align_og_fasta()
        og_fasta_dir = Path("../data") / "Phylogenomic_tree" / "OGs_fasta_aln"
        fasta_files = list(og_fasta_dir.glob("*"))
        self.assertEqual(len(fasta_files), 2137)
        for fasta_file in fasta_files:
            fasta_file_parser = SeqIO.parse(str(fasta_file), "fasta")
            num_records = 0
            for rec in fasta_file_parser:
                self.assertGreater(len(rec.seq), 0)
                num_records += 1
            self.assertEqual(num_records, 3)

    def test_create_supersequence_file(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_f,
            cores,
            fasta_dir,
            out_dir,
            no_keep_fasta,
            tree_model,
            debug=True,
        )
        phylo.create_supersequence_file()
        supersequence_file = Path("../data") / "Phylogenomic_tree" / "supersequence.fa"
        supersequence_parser = SeqIO.parse(str(supersequence_file), "fasta")
        num_records = 0
        for record in supersequence_parser:
            num_records += 1
            self.assertGreater(len(record.seq), 0)
        self.assertEqual(num_records, 3)

    def test_filter_supersequence_aln(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_f,
            cores,
            fasta_dir,
            out_dir,
            no_keep_fasta,
            tree_model,
            debug=True,
        )
        phylo.filter_supersequence_aln()
        supersequence_file = (
            Path("../data") / "Phylogenomic_tree" / "supersequence.fa-gb"
        )
        supersequence_parser = SeqIO.parse(str(supersequence_file), "fasta")
        num_records = 0
        for record in supersequence_parser:
            num_records += 1
            self.assertGreater(len(record.seq), 0)
        htm_supersequence_file = (
            Path("../data") / "Phylogenomic_tree" / "supersequence.fa-gb.htm"
        )
        self.assertEqual(htm_supersequence_file.exists(), False)

    def test_compute_tree(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_f,
            cores,
            fasta_dir,
            out_dir,
            no_keep_fasta,
            tree_model,
            debug=True,
        )
        phylo.compute_tree()
        phylo.move_iqtree_files()
        directory_to_check = Path("../data/Phylogenomic_tree/iqtree")
        files = directory_to_check.glob("*")
        filenames = [f.name for f in files]
        self.assertEqual(len(filenames), 7)
        self.assertIn("supersequence_IQTree2.nwk", filenames)

    def test_run_phylogenomic(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_f,
            cores,
            fasta_dir,
            out_dir,
            no_keep_fasta,
            tree_model,
            debug=True,
        )
        return_code = phylo.run_phylogenomic()
        self.assertEqual(return_code, 0)
