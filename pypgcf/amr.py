from pathlib import Path
from typing import Union

from pandas import DataFrame, read_csv, concat

from pypgcf.utils import (
    calc_avail_dispatchers,
    multiprocess_dispatch,
)


class AMR_analyzer:
    def __init__(
        self,
        *,
        fasta_files_list: list[Path],
        database_dir: Union[Path, None],
        out_dir: Path,
        input_type: bool,
        available_cores: int,
        blast_cores: int,
        concurrent: bool,
        debug: bool,
    ):
        if input_type != "prot":  # either prot or nucl
            self.protein = False
        else:
            self.protein = True
        self.out_dir = out_dir / "AMR"
        self.fasta_files = fasta_files_list
        self.debug = debug
        self.blast_cores = blast_cores
        if concurrent:
            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, blast_cores, avoid_throttle=True
            )
        else:
            self.concurrent_jobs = 0
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1
        self.database_dir = database_dir

    def search_amr(self) -> None:
        cmds = []
        outdir = self.out_dir / "amrfinder"
        outdir.mkdir(exist_ok=True, parents=True)
        for query_file in self.fasta_files:
            outfile = outdir / (query_file.stem + ".txt")

            cmd = ["amrfinder", "--plus", "-i", "-1", "-o", str(outfile)]
            cmd += ["--threads", str(self.blast_cores)]
            if self.database_dir is not None:
                cmd += ["-d", str(self.database_dir)]
            cmd += ["-p" if self.protein else "-n", str(query_file)]
            cmd.append("--debug" if self.debug else "--quiet")
            cmds.append(cmd)

        _ = multiprocess_dispatch(
            "system",
            cmds,
            self.concurrent_jobs,
            show_progress=True,
            description="Scanning for AMR genes",
            debug=self.debug,
        )
        return None

    def parse_results(self):
        # Get the results form previous analysis
        outdir = self.out_dir / "amrfinder"
        files = list(outdir.glob("*"))
        dfs = [DataFrame()] * len(files)
        for idx, f in enumerate(files):
            genome = f.stem
            df = read_csv(f, sep="\t", index_col=0)
            df = df[df["Element type"] == "AMR"]
            df["Genome"] = genome
            dfs[idx] = df

        total_df = concat(dfs)
        fout = self.out_dir / "AMRfinder_results.xlsx"

        # Write to excel
        total_df.to_excel(fout)
