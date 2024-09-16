from pathlib import Path

import pandas as pd
from tqdm import tqdm

from pypgcf.utils import (
    calc_avail_dispatchers,
    create_blastn_cmd,
    create_diamond_blastp_cmd,
    keep_best_homology_hit,
    multiprocess_dispatch,
)


class VF_analyzer:
    def __init__(
        self,
        *,
        fasta_files_list: list[Path],
        database_dir: Path,
        out_dir: Path,
        input_type: bool,
        available_cores: int,
        blast_cores: int,
        evalue: float,
        dmnd_sensitivity: str,
        concurrent: bool,
        debug: bool,
    ):
        self.cores = blast_cores
        self.protein = True
        if input_type != "prot":  # either prot or nucl
            self.protein = False
        self.evalue = evalue
        # The asigned column names will be the same everywhere, need to fix the homology search function
        self.column_names = [
            "qseqid",
            "sseqid",
            "qcovhsp",
            "pident",
            "evalue",
            "bitscore",
        ]
        self.database_dir = database_dir
        if self.protein:
            self.database_file = database_dir / "VFDB_setA_pro"
            self.outfmt = "6 qseqid sseqid qcovhsp pident evalue bitscore"
            self.dmnd_sensitivity = dmnd_sensitivity
        else:
            self.database_file = database_dir / "VFDB_setA_nt"
            self.outfmt = "'6 qaccver saccver qcovhsp pident evalue bitsc"
        self.results_dir = out_dir / "Virulence"
        self.fasta_files = fasta_files_list
        self.debug = debug

        if concurrent:
            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, blast_cores, avoid_throttle=True
            )
        else:
            self.concurrent_jobs = 0
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1

        # TODO: If the input CDS the input should be first translated and then run diamond

    def _read_vf_desc_filt(self) -> dict:
        vf_desc = {}
        with open(self.database_dir / "vfdb_desc.tsv", "r") as rf:
            for line in rf:
                line = line.rstrip()
                if line.startswith("#"):  # Is comment
                    continue
                vf, category, origin, desc = line.split("\t")
                vf_desc[vf] = {
                    "VF_Category": category,
                    "VF_Origin": origin,
                    "VF_Desc": desc,
                }
        return vf_desc

    def execute_homology_search(self) -> None:
        cmds = []
        outdir = self.results_dir / "homology_search"
        outdir.mkdir(exist_ok=True)
        for query_file in self.fasta_files:
            outfile = outdir / (query_file.stem + ".txt")
            if self.protein:
                cmd = create_diamond_blastp_cmd(
                    query_file,
                    self.database_file,
                    outfile,
                    self.dmnd_sensitivity,
                    self.evalue,
                    self.cores,
                    self.outfmt,
                )
            else:
                cmd = create_blastn_cmd(
                    query_file,
                    self.database_file,
                    outfile,
                    self.evalue,
                    self.cores,
                    self.outfmt,
                )
            cmds.append(cmd)

        _ = multiprocess_dispatch(
            "system",
            cmds,
            self.concurrent_jobs,
            show_progress=True,
            description="Scanning for virulence factors",
        )
        return None

    def parse_results(self):
        # Get the results form previous analysis
        outdir = self.results_dir / "homology_search"
        files = list(outdir.glob("*"))
        dfs = [None] * len(files)
        for idx, f in tqdm(enumerate(files), desc="Parsing VF results"):
            df = pd.read_csv(f, sep="\t", names=self.column_names)
            df = keep_best_homology_hit(df)
            mask = (df["qcovhsp"] >= 50) & (
                df["pident"] >= 50
            )  # Get hits 50% identical over 50% of length
            df_filt = df.loc[mask]
            df_filt["QueryGenome"] = f.stem
            dfs[idx] = df_filt
        total_df = pd.concat(dfs)

        # Add VF annotation
        vf_desc = self._read_vf_desc_filt()
        total_df["VF_Category"] = total_df["qseqid"].map(
            lambda x: vf_desc[x]["VF_Category"]
        )
        total_df["VF_Origin"] = total_df["qseqid"].map(
            lambda x: vf_desc[x]["VF_Origin"]
        )
        total_df["VF_Desc"] = total_df["qseqid"].map(lambda x: vf_desc[x]["VF_Desc"])

        # Write to excel
        total_df.to_excel(self.results_dir / "VF_results_50pident50qcov.xlsx")
