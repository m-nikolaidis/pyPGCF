from datetime import datetime
from pathlib import Path
from typing import Union

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pypgcf._vendor.pygblocks import Options, compute_mask, trim_sequence
from pandas import read_csv
from tqdm import tqdm

from pypgcf.utils import (
    calc_avail_dispatchers,
    execute_command,
    multiprocess_dispatch,
    recursive_unlink,
)

# from itaxotools.pygblocks import Options, compute_mask, trim_sequence
# Revert this when the itaxotools recipe is created


class Phylogenomic:
    def __init__(
        self,
        *,
        orthology_matrix_in: Path,
        cores: int,
        available_cores: int,
        fasta_files_list: list[Path],
        out_dir: Path,
        no_keep_fasta: bool,
        tree_method: str,
        iqtree_model: Union[str, None],
        input_type: str,
        debug: bool = False,
    ):
        self.orthology_matrix_in = orthology_matrix_in
        self.out_dir = out_dir / "Phylogenomic_tree"
        self.og_fasta_dir = self.out_dir / "OGs_fasta"
        self.og_fasta_dir_aln = self.out_dir / "OGs_fasta_aln"
        self.iqtree_results_dir = self.out_dir / "IQTree"
        self.cores = cores
        self.no_keep_fasta = no_keep_fasta
        self.debug = debug
        self.tree_method = tree_method
        self.input_type = input_type
        self.fasta_files = fasta_files_list

        if iqtree_model is None:
            self.tree_model = "TEST"
        else:
            self.tree_model = iqtree_model

            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, cores, avoid_throttle=True
            )
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1

    def load_orthology_matrix(self) -> None:
        orthology_matrix = read_csv(self.orthology_matrix_in, sep="\t", index_col=0)
        columns = set(orthology_matrix.columns.tolist())
        genomes_to_keep: list = []
        self.ref = orthology_matrix.index.name

        # Filter unecessary genomes from orthology matrix
        for file in self.fasta_files:
            genome = file.stem
            if genome == self.ref:
                continue
            if genome not in columns:
                continue
            genomes_to_keep.append(genome)

        if self.debug:
            print(
                f"Num of genomes to keep: {len(genomes_to_keep)}. Num fasta files: {len(self.fasta_files)}"
            )

        if len(genomes_to_keep) == 0:
            raise RuntimeError(
                "Phylogenomic: The orthology matrix doesn't contain orthologues from any provided organisms. Please check the names of the fasta files or the input matrix "
            )
        if len(genomes_to_keep) + 1 != len(self.fasta_files):
            raise RuntimeError(
                "Phylogenomic: The orthology matrix doesn't contain orthologues from certain organisms"
            )

        self.orthology_matrix = orthology_matrix[genomes_to_keep]
        self.genomes = genomes_to_keep

        return None

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)
        self.og_fasta_dir.mkdir(exist_ok=True, parents=True)
        self.og_fasta_dir_aln.mkdir(exist_ok=True, parents=True)
        if self.tree_method == "IQTree":
            self.iqtree_results_dir.mkdir(exist_ok=True, parents=True)

    def _replace_empty_with_na_in_orthology_matrix(self):
        self.orthology_matrix = self.orthology_matrix.map(
            lambda x: None if x == "X" else x
        )

    def create_og_fasta(self):
        """
        Create a fasta file for each cluster of orthologous genes
        :return: None
        """
        self._replace_empty_with_na_in_orthology_matrix()
        self.orthology_matrix = self.orthology_matrix.dropna()
        if self.orthology_matrix.shape[0] == 0:
            raise RuntimeError(
                "Couldn't identify orthologous groups with orthologues from all genomes"
            )
        organisms = [self.orthology_matrix.index.name]
        organisms.extend(self.orthology_matrix.columns.tolist())
        fasta_files_dict = {f.stem: f for f in self.fasta_files}
        for org_count, organism in enumerate(organisms):
            fasta_file = fasta_files_dict[organism]
            protein_records = SeqIO.to_dict(
                SeqIO.parse(open(fasta_file, "r"), format="fasta")
            )
            if org_count == 0:
                genes = list(self.orthology_matrix.index)
            else:
                genes = list(self.orthology_matrix[organism].values)
            for x, gene in enumerate(genes):
                info = protein_records[gene]
                info.description = organism
                info.name = ""
                info.id = info.id
                fout = self.og_fasta_dir / ("OG" + str(x) + ".fa")
                fout_handle = open(fout, "a")
                SeqIO.write(info, fout_handle, "fasta")
                fout_handle.close()
            # TODO: If input type is CDS we need to translate before feeding into muscle

    def align_orthologous_groups_fasta(self) -> None:
        """
        Use muscle to align each file of orthologous group
        """
        orthologous_groups_files = list(self.og_fasta_dir.glob("*"))
        commands = []
        for file in orthologous_groups_files:
            cmd = " ".join([
                "muscle",
                "-in",
                str(file),
                "-out",
                str(self.og_fasta_dir_aln / file.name),
                "-quiet",
            ])
            commands.append(cmd)
        _ = multiprocess_dispatch(
            "system",
            commands,
            self.cores,
            show_progress=True,
            description="Aligning OG fasta files",
        )
        return None

    def create_superalignment_file(self) -> None:
        """
        Join the aligned orthologous groups into a superalignment
        :return: None
        """

        def init_superalignment_file(ref, genomes, superseq_file):
            """
            Initialise the superalignment file with empty sequences for each organism
            """
            seqrecord_list = []
            organisms = [ref]
            organisms.extend(genomes[:])
            for org in organisms:
                record = SeqRecord(Seq(""), id=org, name=org, description="")
                seqrecord_list.append(record)
            superseq_file_handle_out = open(superseq_file, "w")
            SeqIO.write(seqrecord_list, superseq_file_handle_out, "fasta")

        superseq_file = self.out_dir / "superalignment.fa"
        init_superalignment_file(self.ref, self.genomes, superseq_file)
        superseq_records = SeqIO.to_dict(AlignIO.read(str(superseq_file), "fasta"))
        aln_files = self.og_fasta_dir_aln.glob("*")

        for aln_file in tqdm(
            aln_files, desc="Joining orthologous groups", ascii=True, leave=True
        ):
            parser = SeqIO.parse(str(aln_file), "fasta")
            aln_records = {
                record.description.replace(record.name + " ", ""): record
                for record in parser
            }  # Rename the keys, need to use the org name
            aln_organisms = list(aln_records.keys())
            for organism in superseq_records:
                if organism in aln_organisms:
                    superseq_records[organism].seq += aln_records[organism].seq
        superseq_seqrecords = list(superseq_records.values())
        superseq_file_handle_out = open(superseq_file, "w")
        SeqIO.write(superseq_seqrecords, superseq_file_handle_out, "fasta")
        return None

    def filter_superalignment(self):
        """Filter the superalignment using Gblocks with default parameters"""
        aln = AlignIO.read(self.out_dir / "superalignment.fa", "fasta")
        options = Options(
            IS=9,  # Minimum Number Of Sequences For A Conserved Position
            FS=14,  # Minimum Number Of Sequences For A Flank Position
            CP=8,  # Maximum Number Of Contiguous Nonconserved Positions
            BL1=10,  # Minimum Length Of A Block, 1st iteration
            BL2=10,  # Minimum Length Of A Block, 2nd iteration
            GT=0,  # Maximum Number of Allowed Gaps For Any Position
            GC="-",  # Definition of Gap Characters
        )
        mask = compute_mask((record.seq for record in aln), options)
        for record in aln:
            record.seq = Seq(trim_sequence(record.seq, mask))
        with open(self.out_dir / "superalignment.fa-gb", "w") as wf:
            AlignIO.write(aln, wf, "fasta")
        return None

    def compute_tree_bionj(self):
        """
        Compute the phylogenomic tree using bionj
        """
        raise NotImplementedError("Not implemented yet")

    def compute_tree_fasttree(self):
        """
        Compute the phylogenomic tree using FastTree
        """
        superalignment_file = self.out_dir / "superalignment.fa-gb"
        if self.input_type == "nt" or self.input_type == "cds":
            cmd = f"fasttree -nt -gtr {superalignment_file}"
        else:  # prot
            cmd = f"fasttree -lg {superalignment_file}"
        if not self.debug:
            cmd += " -quiet"
        fout = self.out_dir / "superalignment_Fasttree.nwk"
        cmd += f" > {fout}"
        _ = execute_command(cmd)

    def compute_tree_iqtree(self):
        """
        Compute the phylogenomic tree using IQtree2
        """
        superalignment_file = self.out_dir / "superalignment.fa-gb"
        cmd = "iqtree2 -m {} -merit AIC --alrt 1000 -T {} -s {}".format(
            self.tree_model, self.cores, superalignment_file
        )
        if not self.debug:
            cmd += " --quiet"
        _ = execute_command(cmd)

    def move_iqtree_files(self):
        files = list(self.out_dir.glob("superalignment.fa-gb.*"))
        for f in files:
            renamed = self.iqtree_results_dir / f.name
            if f.suffix == ".treefile":
                renamed = self.out_dir / "superalignment_IQTree2.nwk"
            f.rename(renamed)

    def clean_fasta_files(self):
        directories_to_clean = [self.og_fasta_dir, self.og_fasta_dir_aln]
        for directory in directories_to_clean:
            recursive_unlink(directory)

    def run_phylogenomic(self) -> Union[None, int]:
        print(
            f"Loading orthology matrix: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.load_orthology_matrix()
        self.setup_directories()
        print(
            f"Creating fasta files of each orthologous group: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.create_og_fasta()
        print(
            f"Aligning fasta files of each orthologous group: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.align_orthologous_groups_fasta()
        print(
            f"Creating super-alignment file: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.create_superalignment_file()
        print(
            f"Filtering super-alignment file: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.filter_superalignment()
        print(
            f"Computing phylogenomic tree: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        if self.tree_method == "Fasttree":
            self.compute_tree_fasttree()
        if self.tree_method == "NJ":
            self.compute_tree_bionj()
        if self.tree_method == "IQTree":
            self.compute_tree_iqtree()
            self.move_iqtree_files()
        print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        if self.no_keep_fasta:
            self.clean_fasta_files()
        if self.debug:
            return 0
