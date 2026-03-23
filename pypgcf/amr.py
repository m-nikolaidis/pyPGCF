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
        if database_dir is None:
            self.database_dir = None
        else:
            self.database_dir = database_dir / "AMR"

    def search_amr(self) -> None:
        cmds = []
        outdir = self.out_dir / "amrfinder"
        outdir.mkdir(exist_ok=True, parents=True)
        for query_file in self.fasta_files:
            outfile = outdir / (query_file.stem + ".txt")

            if self.database_dir is None:
                cmd = (
                    f"amrfinder --plus -i -1 -o {outfile} --threads {self.blast_cores}"
                )
            else:
                cmd = f"amrfinder --plus -i -1 -o {outfile} --threads {self.blast_cores} -d {self.database_dir}"
            if self.protein:
                cmd += f" -p {query_file}"
            else:
                cmd += f" -n {query_file}"
            if not self.debug:
                cmd += " --quiet"
            else:
                cmd += " --debug"
            cmds.append(cmd)

        _ = multiprocess_dispatch(
            "system",
            cmds,
            self.concurrent_jobs,
            show_progress=True,
            description="Scanning for AMR genes",
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
