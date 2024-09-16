from datetime import datetime
from pathlib import Path
from typing import Iterable, Union

from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

from pypgcf.utils import (
    execute_command,
    multiprocess_dispatch,
    calc_avail_dispatchers,
    recursive_unlink,
)


class Orthologues_identifier:
    def __init__(
        self,
        *,
        fasta_files_list: list[Path],
        out_dir: Path,
        ref: Union[str, None],
        ref_list: Union[str, None],
        input_type: str,
        available_cores: int,
        blast_cores: int,
        concurrent: bool,
        evalue: float,
        dmnd_sensitivity: str,
        no_filter: bool,
    ):
        self.ref = ref
        self.ref_list = ref_list
        self.out_dir = out_dir / "Orthologues"
        self.blast_cores = blast_cores
        self.blast_evalue = evalue
        self.dmnd_sensitivity = dmnd_sensitivity
        self.fasta_files = fasta_files_list
        self.input_type = input_type
        self.no_filter = no_filter
        self.blast_db_dir = out_dir / "Blast_DB"

        if concurrent:
            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, blast_cores, avoid_throttle=True
            )
        else:
            self.concurrent_jobs = 0
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1

    def _create_directories(self, ref: str) -> None:
        directories = [
            self.out_dir / ref / "Blast_results",
            self.out_dir / ref / "Best_reciprocal_hits",
            self.blast_db_dir,
        ]
        for directory in directories:
            directory.mkdir(exist_ok=True, parents=True)

    def _get_blast_binaries(self) -> None:
        diamond_bin = "diamond"
        blastn_bin = "blastn"
        makeblastdb_bin = "makeblastdb"

        self.blast_bin = diamond_bin
        if self.input_type == "CDS":
            self.blast_bin = blastn_bin
            self.blast_db_bin = makeblastdb_bin

    def create_blast_db(self) -> None:
        """
        Create diamond/blast database for a fasta file
        """
        for fasta_file in tqdm(
            self.fasta_files, desc="Preparing DIAMOND/BLAST database", ascii=True
        ):
            database_f = self.blast_db_dir / fasta_file.stem
            cmd = f"{self.blast_bin} makedb --in {fasta_file} --quiet --db {database_f} --threads {self.blast_cores}"  # DIAMOND
            if self.input_type == "CDS":
                cmd = f"{self.blast_db_bin} -in {fasta_file} -dbtype nucl -out {database_f}"
            execute_command(cmd)

    def _create_blast_cmd(
        self, fasta_file: Path, database_f: Path, out_file: Path
    ) -> str:
        # DIAMOND case first
        added_sensitivity = " --very-sensitive"
        if self.dmnd_sensitivity == "sensitive":
            added_sensitivity = " --sensitive"
        if self.dmnd_sensitivity == "mid_sensitive":
            added_sensitivity = " --mid-sensitive"
        if self.dmnd_sensitivity == "more_sensitive":
            added_sensitivity = " --more-sensitive"
        if self.dmnd_sensitivity == "ultra_sensitive":
            added_sensitivity = " --ultra-sensitive"
        cmd = f"{self.blast_bin} blastp --query {fasta_file} --quiet --db {database_f} --outfmt 6 --out {out_file} --evalue {self.blast_evalue} --threads {self.blast_cores} "
        cmd += added_sensitivity
        if self.input_type == "nucl":
            cmd = f"{self.blast_bin} -query {fasta_file} -db {database_f} -outfmt 6 -out {out_file} -evalue {self.blast_evalue} -num_threads {self.blast_cores} "
        return cmd

    def setup(self) -> None:
        if self.ref is not None:
            self._create_directories(self.ref)
        if self.ref_list is not None:
            with open(self.ref_list, "r") as rf:
                for ref in rf.readlines():
                    ref = ref.rstrip()
                    self._create_directories(ref)
        self._get_blast_binaries()

    def perform_reciprocal_blast(self, ref: str):
        """ """
        # Check if the ref file exists
        ref_fasta = None
        for fasta_file in self.fasta_files:
            if ref in fasta_file.name:
                ref_fasta = fasta_file
        if ref_fasta is None:
            raise FileNotFoundError(f"{ref} is not in input fasta file names")

        ref_db = self.blast_db_dir / ref

        cmds = []
        for fasta_file in self.fasta_files:
            if ref == fasta_file.stem:
                continue
            database_f = self.blast_db_dir / fasta_file.stem
            fout_f = (
                self.out_dir
                / ref
                / "Blast_results"
                / (fasta_file.name + "_vs_reference.txt")
            )
            fout_r = (
                self.out_dir
                / ref
                / "Blast_results"
                / (fasta_file.name + "_vs_reference_reverse.txt")
            )
            forward_blast = self._create_blast_cmd(ref_fasta, database_f, fout_f)
            reverse_blast = self._create_blast_cmd(fasta_file, ref_db, fout_r)
            cmds.append(forward_blast)
            cmds.append(reverse_blast)

        _ = multiprocess_dispatch(
            "system",
            cmds,
            self.concurrent_jobs,
            show_progress=True,
            description="Performing reciprocal homology search",
        )
        return ref_fasta

    def _get_best_subject(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Helper function to keep the best subject hit per query
        Input: Pandas Data frame
        Return: Filtered Pandas data frame
        """
        indeces_to_keep = []
        first_col = df.columns[0]
        for _, grp_df in df.groupby(first_col):
            indeces_to_keep.append(grp_df["Bitscore"].idxmax())
        return df.loc[indeces_to_keep].reset_index(drop=True)

    def _create_reciprocal_matrix(self, ref_vs_query_df, query_vs_ref_df):
        """
        Helper function which compare the two dataframes given line by line
        to identify reciprocal hits
        Input: The one way and the reverse way dataframe
        Return a new data frame which contains only the reciprocal hits
        """
        tmp_df = pd.merge(ref_vs_query_df, query_vs_ref_df, on="RefSeq")
        df = tmp_df.drop(tmp_df[tmp_df.QuerySeq_x != tmp_df.QuerySeq_y].index)
        case1 = df[df["Pident_x"] >= df["Pident_y"]]
        case1 = case1[["RefSeq", "QuerySeq_x", "Pident_x", "Evalue_x", "Bitscore_x"]]
        case1 = case1.rename(
            columns={
                "QuerySeq_x": "QuerySeq",
                "Pident_x": "Pident",
                "Evalue_x": "Evalue",
                "Bitscore_x": "Bitscore",
            }
        )
        case2 = df[df["Pident_x"] < df["Pident_y"]]
        case2 = case2.rename(
            columns={
                "QuerySeq_y": "QuerySeq",
                "Pident_y": "Pident",
                "Evalue_y": "Evalue",
                "Bitscore_y": "Bitscore",
            }
        )
        df = pd.concat([case1, case2])
        return df[["RefSeq", "QuerySeq", "Pident"]]

    def _orthologue_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Helper function to create the distribution of percent identities for all the genes in the data frame
        Filter the hits that have perc. identity lower than -2 SD
        """
        pid_std = df.Pident.std(axis=0, numeric_only=True)
        pid_mean = df.Pident.mean(axis=0, numeric_only=True)
        drop_condition = pid_mean - (2 * pid_std)
        return df.drop(df[df.Pident < drop_condition].index)

    def parse_blast_results(self, ref: str):
        """
        Reads the txt blast output and create the reciprocal table
        First read the files and filter the redundant information for each protein hit, keep the best hit
        Create the reciprocal blast matrices and write into files
        Input: The ref and fasta files
        Output: Best reciprocal hits per strain
        Return: None
        """

        fasta_files = self.fasta_files
        for fasta_file in tqdm(
            fasta_files, ascii=True, leave=True, desc="Parsing BLAST output"
        ):
            if ref in fasta_file.stem:
                continue
            ref_vs_query_file = (
                self.out_dir
                / ref
                / "Blast_results"
                / (fasta_file.name + "_vs_reference.txt")
            )
            query_vs_ref_file = (
                self.out_dir
                / ref
                / "Blast_results"
                / (fasta_file.name + "_vs_reference_reverse.txt")
            )

            c1 = ["RefSeq", "QuerySeq", "Pident", "Evalue", "Bitscore"]
            c2 = ["QuerySeq", "RefSeq", "Pident", "Evalue", "Bitscore"]

            ref_vs_query_df = pd.read_csv(
                ref_vs_query_file, sep="\t", names=c1, usecols=[0, 1, 2, 10, 11]
            )
            query_vs_ref_df = pd.read_csv(
                query_vs_ref_file, sep="\t", names=c2, usecols=[0, 1, 2, 10, 11]
            )

            # Blast output
            # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
            ref_vs_query_df = self._get_best_subject(ref_vs_query_df)
            query_vs_ref_df = self._get_best_subject(query_vs_ref_df)

            reciprocal_df = self._create_reciprocal_matrix(
                ref_vs_query_df, query_vs_ref_df
            )
            if not self.no_filter:  # No_filtering is false
                reciprocal_df = self._orthologue_filter(reciprocal_df)
            # Create reciprocal files with headers
            rec_fout = (
                self.out_dir
                / ref
                / "Best_reciprocal_hits"
                / (fasta_file.stem + "_BRH.txt")
            )
            reciprocal_df.to_csv(rec_fout, sep="\t", index=False)

    def _init_orthology_matrix(
        self, refgenes: Iterable, reciprocal_files: Iterable
    ) -> dict:
        orthology_matrix_dict = {gene.id: {} for gene in refgenes}
        for reciprocal_file in reciprocal_files:
            query_genome = reciprocal_file.stem.replace("_BRH", "")
            for gene in orthology_matrix_dict:
                orthology_matrix_dict[gene][query_genome] = "X"
        return orthology_matrix_dict

    def create_orthology_matrix(self, ref: str, ref_fasta: Path) -> None:
        """
        Create the initial cog matrix with empty (NaN) vectors and replace with the true values
        :return: None
        """
        # Pass all genes of the Refseq into a list

        refgenes = SeqIO.parse(str(ref_fasta), "fasta")
        reciprocal_f_dir = self.out_dir / ref / "Best_reciprocal_hits"
        reciprocal_files = list(reciprocal_f_dir.glob("*"))

        init_orthology_matrix_dict = self._init_orthology_matrix(
            refgenes, reciprocal_files
        )

        for reciprocal_file in reciprocal_files:
            query_genome = reciprocal_file.stem.replace("_BRH", "")
            header = True
            for lines in open(reciprocal_file, mode="r"):
                if header:
                    header = False
                    continue
                refseq, queryseq, *_ = lines.rstrip().split("\t")
                init_orthology_matrix_dict[refseq][query_genome] = queryseq
        orthology_matrix_df = pd.DataFrame.from_dict(
            init_orthology_matrix_dict, orient="index"
        )
        orthology_matrix_df.index.name = ref
        orthology_matrix_f = self.out_dir / ref / "OGmatrix.csv"
        orthology_matrix_df.to_csv(orthology_matrix_f, sep="\t")

    def calculate_orthologues(self):
        refs = []
        self.setup()
        if self.ref is not None:
            ref = self.ref
            refs.append(ref)
        if self.ref_list is not None:  # Duplicate code in setup function
            list_in = open(self.ref_list, "r")
            for ref in list_in:
                ref = ref.rstrip()
                refs.append(ref)
            list_in.close()
        self.create_blast_db()
        for idx, ref in enumerate(refs):
            if idx > 0:
                print("-" * 200)
            print(
                f"Calulating orthologues with reference {ref}: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
            )
            ref_fasta = self.perform_reciprocal_blast(ref)
            self.parse_blast_results(ref)
            print(
                f"Creating orthology matrix: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
            )
            self.create_orthology_matrix(ref, ref_fasta)
            print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        recursive_unlink(self.blast_db_dir)
