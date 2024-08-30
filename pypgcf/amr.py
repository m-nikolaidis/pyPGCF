from pathlib import Path

import pandas as pd

from pypgcf.utils import (
    calc_avail_dispatchers,
    create_blastn_cmd,
    create_diamond_blastp_cmd,
    multiprocess_dispatch,
)


class AMR_analyzer:
    def __init__(
        self,
        *,
        fasta_dir: Path,
        database_dir: Path,
        results_dir: Path,
        input_type: bool,
        available_cores: int,
        blast_cores: int,
        evalue: float,
        dmnd_sensitivity: str,
        concurrent: bool,
        debug: bool,
    ):
        self.cpus = blast_cores
        self.protein = True
        if input_type != "prot":  # either prot or nucl
            self.protein = False
        self.results_dir = results_dir
        self.fasta_dir = fasta_dir
        self.debug = debug
        if concurrent:
            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, blast_cores, avoid_throttle=True
            )
        else:
            self.concurrent_jobs = 0
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1

        # TODO: If the input CDS the input should be first translated first?

    def search_amr(self) -> None:
        query_files = list(self.fasta_dir.glob("*"))
        cmds = []
        outdir = self.results_dir / "AMR"
        outdir.mkdir(exist_ok=True)
        for query_file in query_files:
            outfile = outdir / query_file.stem + ".txt"
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
        outdir = self.results_dir / "AMR"
        files = list(outdir.glob("*"))
        dfs = [pd.DataFrame()] * len(files)
        for idx, f in enumerate(files):
            genome = f.stem
            df = pd.read_csv(f, sep="\t", index_col=0)
            df = df[df["Element type"] == "AMR"]
            df["Genome"] = genome
            dfs[idx] = df

        total_df = pd.concat(dfs)
        fout = self.results_dir / "AMR" / "AMRfinder_results.xlsx"

        # Write to excel
        total_df.to_excel(fout)
