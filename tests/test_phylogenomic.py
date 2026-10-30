from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from context import phylogenomic


class TestCodonProjection(unittest.TestCase):
    def test_project_protein_alignment_to_codons_maps_gaps_and_preserves_cds(self):
        nucleotides = {
            "a": SeqRecord(Seq("ATGAAATTT"), id="a", description="sample a"),
            "b": SeqRecord(Seq("ATGTTT"), id="b", description="sample b"),
        }
        aligned = [
            SeqRecord(Seq("MKF"), id="a"),
            SeqRecord(Seq("M-F"), id="b"),
        ]

        result = phylogenomic.project_protein_alignment_to_codons(
            aligned, nucleotides
        )

        self.assertEqual(
            [str(record.seq) for record in result],
            ["ATGAAATTT", "ATG---TTT"],
        )
        self.assertEqual(str(result[1].seq).replace("-", ""), "ATGTTT")
        self.assertEqual(result[0].description, "sample a")

    def test_project_protein_alignment_to_codons_rejects_partial_cds(self):
        nucleotides = {"a": SeqRecord(Seq("ATGA"), id="a")}
        aligned = [SeqRecord(Seq("M"), id="a")]

        with self.assertRaisesRegex(ValueError, "a.*divisible by three"):
            phylogenomic.project_protein_alignment_to_codons(aligned, nucleotides)

    def test_project_protein_alignment_to_codons_requires_matching_identifiers(self):
        nucleotides = {"a": SeqRecord(Seq("ATG"), id="a")}
        aligned = [SeqRecord(Seq("M"), id="b")]

        with self.assertRaisesRegex(ValueError, "identifiers do not match"):
            phylogenomic.project_protein_alignment_to_codons(aligned, nucleotides)

    def test_project_protein_alignment_to_codons_preserves_stop_and_ambiguous_codons(
        self,
    ):
        nucleotides = {
            "a": SeqRecord(Seq("ATGNNNTAA"), id="a", description="ambiguous")
        }
        aligned = [SeqRecord(Seq("M-X*"), id="a")]

        result = phylogenomic.project_protein_alignment_to_codons(
            aligned, nucleotides
        )

        self.assertEqual(str(result[0].seq), "ATG---NNNTAA")
        self.assertEqual(str(result[0].seq).replace("-", ""), "ATGNNNTAA")


class TestCdsOgAlignment(unittest.TestCase):
    def _write_source(self, source):
        SeqIO.write(
            [
                SeqRecord(Seq("ATGAAATTT"), id="a", description="sample a"),
                SeqRecord(Seq("ATGTTT"), id="b", description="sample b"),
            ],
            source,
            "fasta",
        )

    def test_align_cds_orthologous_group_uses_temporary_protein_files(self):
        captured_temporary_paths = []
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            source = tmp_path / "OG1.fa"
            destination = tmp_path / "aligned" / "OG1.fa"
            self._write_source(source)

            def fake_muscle(command, debug=False):
                arguments = list(command)
                protein_input = Path(arguments[arguments.index("-in") + 1])
                protein_output = Path(arguments[arguments.index("-out") + 1])
                captured_temporary_paths.extend([protein_input, protein_output])
                with open(protein_input) as protein_handle:
                    translated = SeqIO.to_dict(SeqIO.parse(protein_handle, "fasta"))
                self.assertEqual(str(translated["a"].seq), "MKF")
                self.assertEqual(str(translated["b"].seq), "MF")
                SeqIO.write(
                    [
                        SeqRecord(Seq("MKF"), id="a"),
                        SeqRecord(Seq("M-F"), id="b"),
                    ],
                    protein_output,
                    "fasta",
                )
                return 0

            with patch.object(
                phylogenomic, "execute_command", side_effect=fake_muscle
            ):
                result = phylogenomic._align_cds_orthologous_group(
                    source, destination
                )

            self.assertEqual(result, 0)
            with open(destination) as destination_handle:
                aligned_sequences = [
                    str(record.seq)
                    for record in SeqIO.parse(destination_handle, "fasta")
                ]
            self.assertEqual(
                aligned_sequences,
                ["ATGAAATTT", "ATG---TTT"],
            )
            self.assertTrue(destination.exists())
            self.assertTrue(
                all(not path.exists() for path in captured_temporary_paths)
            )

    def test_align_cds_orthologous_group_cleans_up_after_muscle_failure(self):
        captured_temporary_paths = []
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            source = tmp_path / "OG1.fa"
            destination = tmp_path / "aligned" / "OG1.fa"
            self._write_source(source)

            def fail_muscle(command, debug=False):
                arguments = list(command)
                captured_temporary_paths.extend(
                    [
                        Path(arguments[arguments.index("-in") + 1]),
                        Path(arguments[arguments.index("-out") + 1]),
                    ]
                )
                return 1

            with patch.object(
                phylogenomic, "execute_command", side_effect=fail_muscle
            ):
                with self.assertRaisesRegex(RuntimeError, "MUSCLE failed.*OG1.fa"):
                    phylogenomic._align_cds_orthologous_group(source, destination)

            self.assertFalse(destination.exists())
            self.assertTrue(
                all(not path.exists() for path in captured_temporary_paths)
            )

    def test_align_cds_orthologous_group_requires_alignment_output(self):
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            source = tmp_path / "OG1.fa"
            destination = tmp_path / "aligned" / "OG1.fa"
            self._write_source(source)

            with patch.object(phylogenomic, "execute_command", return_value=0):
                with self.assertRaisesRegex(
                    RuntimeError, "did not create a protein alignment"
                ):
                    phylogenomic._align_cds_orthologous_group(source, destination)

            self.assertFalse(destination.exists())

    def test_align_cds_orthologous_group_rejects_duplicate_source_ids(self):
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            source = tmp_path / "OG1.fa"
            destination = tmp_path / "aligned" / "OG1.fa"
            SeqIO.write(
                [
                    SeqRecord(Seq("ATG"), id="duplicate"),
                    SeqRecord(Seq("TTT"), id="duplicate"),
                ],
                source,
                "fasta",
            )

            with self.assertRaisesRegex(ValueError, "Duplicate CDS identifier"):
                phylogenomic._align_cds_orthologous_group(source, destination)

            self.assertFalse(destination.exists())

    def test_align_cds_orthologous_group_rejects_unknown_aligned_ids(self):
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            source = tmp_path / "OG1.fa"
            destination = tmp_path / "aligned" / "OG1.fa"
            self._write_source(source)

            def write_unknown_alignment(command, debug=False):
                arguments = list(command)
                protein_output = Path(arguments[arguments.index("-out") + 1])
                SeqIO.write(
                    [SeqRecord(Seq("MKF"), id="unknown")],
                    protein_output,
                    "fasta",
                )
                return 0

            with patch.object(
                phylogenomic,
                "execute_command",
                side_effect=write_unknown_alignment,
            ):
                with self.assertRaisesRegex(ValueError, "identifiers do not match"):
                    phylogenomic._align_cds_orthologous_group(source, destination)

            self.assertFalse(destination.exists())


class TestPhylogenomicAlignmentRouting(unittest.TestCase):
    def test_routes_cds_alignment_jobs_to_codon_alignment_worker(self):
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            og_fasta_dir = tmp_path / "OGs_fasta"
            og_fasta_dir_aln = tmp_path / "OGs_fasta_aln"
            og_fasta_dir.mkdir()
            og_fasta_dir_aln.mkdir()
            sources = [og_fasta_dir / "OG1.fa", og_fasta_dir / "OG2.fa"]
            for source in sources:
                source.touch()

            runner = phylogenomic.Phylogenomic.__new__(phylogenomic.Phylogenomic)
            runner.input_type = "cds"
            runner.og_fasta_dir = og_fasta_dir
            runner.og_fasta_dir_aln = og_fasta_dir_aln
            runner.concurrent_jobs = 2
            runner.debug = False

            with patch.object(phylogenomic, "multiprocess_dispatch") as dispatch:
                runner.align_orthologous_groups_fasta()

            dispatched_function, jobs, processes = dispatch.call_args.args[:3]
            self.assertIs(
                dispatched_function,
                phylogenomic._align_cds_orthologous_group_job,
            )
            self.assertEqual(processes, 2)
            self.assertEqual(
                jobs,
                [
                    (source, og_fasta_dir_aln / source.name, False)
                    for source in sources
                ],
            )
            self.assertEqual(
                dispatch.call_args.kwargs["description"],
                "Aligning translated CDS orthologous groups",
            )

    def test_protein_alignment_keeps_direct_muscle_dispatch(self):
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            og_fasta_dir = tmp_path / "OGs_fasta"
            og_fasta_dir_aln = tmp_path / "OGs_fasta_aln"
            og_fasta_dir.mkdir()
            og_fasta_dir_aln.mkdir()
            source = og_fasta_dir / "OG1.fa"
            source.touch()

            runner = phylogenomic.Phylogenomic.__new__(phylogenomic.Phylogenomic)
            runner.input_type = "prot"
            runner.og_fasta_dir = og_fasta_dir
            runner.og_fasta_dir_aln = og_fasta_dir_aln
            runner.cores = 3
            runner.debug = False

            with patch.object(phylogenomic, "multiprocess_dispatch") as dispatch:
                runner.align_orthologous_groups_fasta()

            self.assertEqual(dispatch.call_args.args[0], "system")
            self.assertEqual(dispatch.call_args.args[2], 3)
            self.assertEqual(
                dispatch.call_args.args[1],
                [
                    [
                        "muscle",
                        "-in",
                        str(source),
                        "-out",
                        str(og_fasta_dir_aln / source.name),
                        "-quiet",
                    ]
                ],
            )

    def test_protein_alignment_failure_raises(self):
        with TemporaryDirectory() as tmp_dir:
            og_fasta_dir = Path(tmp_dir)
            (og_fasta_dir / "OG1.fa").touch()
            (og_fasta_dir / "OG2.fa").touch()
            runner = phylogenomic.Phylogenomic.__new__(phylogenomic.Phylogenomic)
            runner.input_type = "prot"
            runner.og_fasta_dir = og_fasta_dir
            runner.og_fasta_dir_aln = og_fasta_dir
            runner.cores, runner.debug = 1, False
            failed_one = patch.object(
                phylogenomic, "multiprocess_dispatch", return_value=[0, 1]
            )
            with failed_one, self.assertRaisesRegex(RuntimeError, "1 orthologous"):
                runner.align_orthologous_groups_fasta()

    def test_cds_alignment_end_to_end_writes_only_projected_nucleotides(self):
        captured_temporary_paths = []
        with TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            og_fasta_dir = tmp_path / "OGs_fasta"
            og_fasta_dir_aln = tmp_path / "OGs_fasta_aln"
            og_fasta_dir.mkdir()
            og_fasta_dir_aln.mkdir()
            source = og_fasta_dir / "OG1.fa"
            original_sequences = {"a": "ATGAAATTT", "b": "ATGTTT"}
            SeqIO.write(
                [
                    SeqRecord(Seq(sequence), id=record_id)
                    for record_id, sequence in original_sequences.items()
                ],
                source,
                "fasta",
            )

            runner = phylogenomic.Phylogenomic.__new__(phylogenomic.Phylogenomic)
            runner.input_type = "cds"
            runner.og_fasta_dir = og_fasta_dir
            runner.og_fasta_dir_aln = og_fasta_dir_aln
            runner.concurrent_jobs = 1
            runner.debug = False

            def fake_muscle(command, debug=False):
                arguments = list(command)
                protein_input = Path(arguments[arguments.index("-in") + 1])
                protein_output = Path(arguments[arguments.index("-out") + 1])
                captured_temporary_paths.extend([protein_input, protein_output])
                SeqIO.write(
                    [
                        SeqRecord(Seq("MKF"), id="a"),
                        SeqRecord(Seq("M-F"), id="b"),
                    ],
                    protein_output,
                    "fasta",
                )
                return 0

            def run_jobs(function, jobs, num_procs, **kwargs):
                return [function(job) for job in jobs]

            with (
                patch.object(
                    phylogenomic, "execute_command", side_effect=fake_muscle
                ),
                patch.object(
                    phylogenomic,
                    "multiprocess_dispatch",
                    side_effect=run_jobs,
                ),
            ):
                runner.align_orthologous_groups_fasta()

            destination = og_fasta_dir_aln / source.name
            with open(destination) as destination_handle:
                aligned = list(SeqIO.parse(destination_handle, "fasta"))
            self.assertEqual(
                [str(record.seq) for record in aligned],
                ["ATGAAATTT", "ATG---TTT"],
            )
            for record in aligned:
                self.assertEqual(len(record.seq) % 3, 0)
                self.assertEqual(
                    str(record.seq).replace("-", ""),
                    original_sequences[record.id],
                )
            self.assertTrue(destination.exists())
            self.assertTrue(
                all(not path.exists() for path in captured_temporary_paths)
            )


class TestModule(unittest.TestCase):
    def _cleanup_fasta(self):
        og_fasta_dir = Path("../data") / "Phylogenomic_tree" / "OGs_fasta"
        for fasta_file in og_fasta_dir.glob("*"):
            fasta_file.unlink()  # For cleaning
        og_fasta_dir.rmdir()

    def _iqtree_command(self, fast):
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=Path("matrix.csv"),
            cores=2,
            available_cores=2,
            fasta_files_list=[],
            out_dir=Path("out"),
            no_keep_fasta=False,
            tree_method="IQTree",
            iqtree_model="TEST",
            input_type="prot",
            fast=fast,
            debug=True,
        )
        with patch("pypgcf.phylogenomic.execute_command", return_value=0) as run:
            phylo.compute_tree_iqtree()
        return run.call_args[0][0]

    def test_filter_superalignment_keeps_conserved_columns_for_few_genomes(self):
        with TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            protein = Seq("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK")
            records = [SeqRecord(protein, id=g) for g in ["G1", "G2", "G3"]]
            SeqIO.write(records, tmp / "superalignment.fa", "fasta")
            phylo = phylogenomic.Phylogenomic.__new__(phylogenomic.Phylogenomic)
            phylo.out_dir = tmp
            phylo.filter_superalignment()
            filtered = list(SeqIO.parse(tmp / "superalignment.fa-gb", "fasta"))
        # Fixed thresholds for about 16 genomes removed every column here
        self.assertEqual(len(filtered[0].seq), len(records[0].seq))

    def test_iqtree_command_uses_fast_only_when_requested(self):
        self.assertNotIn("--fast", self._iqtree_command(fast=False))
        self.assertIn("--fast", self._iqtree_command(fast=True))

    def test_iqtree_failure_raises(self):
        phylo = phylogenomic.Phylogenomic.__new__(phylogenomic.Phylogenomic)
        phylo.out_dir, phylo.tree_model, phylo.cores = Path("out"), "TEST", 1
        phylo.fast, phylo.debug = False, False
        with patch("pypgcf.phylogenomic.execute_command", return_value=2):
            with self.assertRaises(RuntimeError):
                phylo.compute_tree_iqtree()

    @unittest.skip("needs ../data (OGmatrix and proteome FASTAs), not included in the repository")
    def test_create_og_fasta(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=orthology_matrix_f,
            cores=cores,
            available_cores=cores,
            fasta_files_list=list(fasta_dir.glob("*")),
            out_dir=out_dir,
            no_keep_fasta=no_keep_fasta,
            tree_method="iqtree",
            iqtree_model=tree_model,
            input_type="prot",
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

    @unittest.skip("needs ../data (OGmatrix and proteome FASTAs), not included in the repository")
    def test_align_og_fasta(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=orthology_matrix_f,
            cores=cores,
            available_cores=cores,
            fasta_files_list=list(fasta_dir.glob("*")),
            out_dir=out_dir,
            no_keep_fasta=no_keep_fasta,
            tree_method="iqtree",
            iqtree_model=tree_model,
            input_type="prot",
            debug=True,
        )
        phylo.align_orthologous_groups_fasta()
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

    @unittest.skip("needs ../data (OGmatrix and proteome FASTAs), not included in the repository")
    def test_create_supersequence_file(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=orthology_matrix_f,
            cores=cores,
            available_cores=cores,
            fasta_files_list=list(fasta_dir.glob("*")),
            out_dir=out_dir,
            no_keep_fasta=no_keep_fasta,
            tree_method="iqtree",
            iqtree_model=tree_model,
            input_type="prot",
            debug=True,
        )
        phylo.create_superalignment_file()
        supersequence_file = Path("../data") / "Phylogenomic_tree" / "supersequence.fa"
        supersequence_parser = SeqIO.parse(str(supersequence_file), "fasta")
        num_records = 0
        for record in supersequence_parser:
            num_records += 1
            self.assertGreater(len(record.seq), 0)
        self.assertEqual(num_records, 3)

    @unittest.skip("needs ../data (OGmatrix and proteome FASTAs), not included in the repository")
    def test_filter_supersequence_aln(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=orthology_matrix_f,
            cores=cores,
            available_cores=cores,
            fasta_files_list=list(fasta_dir.glob("*")),
            out_dir=out_dir,
            no_keep_fasta=no_keep_fasta,
            tree_method="iqtree",
            iqtree_model=tree_model,
            input_type="prot",
            debug=True,
        )
        phylo.filter_superalignment()
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

    @unittest.skip("needs ../data (OGmatrix and proteome FASTAs), not included in the repository")
    def test_compute_tree(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=orthology_matrix_f,
            cores=cores,
            available_cores=cores,
            fasta_files_list=list(fasta_dir.glob("*")),
            out_dir=out_dir,
            no_keep_fasta=no_keep_fasta,
            tree_method="iqtree",
            iqtree_model=tree_model,
            input_type="prot",
            debug=True,
        )
        phylo.compute_tree_iqtree()
        phylo.move_iqtree_files()
        directory_to_check = Path("../data/Phylogenomic_tree/iqtree")
        files = directory_to_check.glob("*")
        filenames = [f.name for f in files]
        self.assertEqual(len(filenames), 7)
        self.assertIn("supersequence_IQTree2.nwk", filenames)

    @unittest.skip("needs ../data (OGmatrix and proteome FASTAs), not included in the repository")
    def test_run_phylogenomic(self):
        fasta_dir = Path("../data/Protein_fasta_files")
        orthology_matrix_f = Path("../data/OGmatrix.csv")
        cores = 4
        out_dir = Path("../data")
        no_keep_fasta = False
        tree_model = None
        phylo = phylogenomic.Phylogenomic(
            orthology_matrix_in=orthology_matrix_f,
            cores=cores,
            available_cores=cores,
            fasta_files_list=list(fasta_dir.glob("*")),
            out_dir=out_dir,
            no_keep_fasta=no_keep_fasta,
            tree_method="iqtree",
            iqtree_model=tree_model,
            input_type="prot",
            debug=True,
        )
        return_code = phylo.run_phylogenomic()
        self.assertEqual(return_code, 0)
