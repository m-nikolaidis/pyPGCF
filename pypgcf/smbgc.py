import json
from pathlib import Path
from typing import Union

import pandas as pd
from tqdm import tqdm

from pypgcf.utils import (
    check_if_file_exists,
    execute_command,
    multiprocess_dispatch,
    calc_avail_dispatchers,
)


class smBGCLocalRunner:
    def __init__(
        self,
        *,
        genome_fasta_files: Path,
        out_dir: Path,
        strictness: str,
        genefinding_tool: str,
        database_dir: Union[Path, None],
        cores: int,
        available_cores: int,
        concurrent: bool,
        debug: bool,
    ):
        self.antismash_raw_results_dir = out_dir / "smBGC" / "antismash_raw_results"
        self.cores = cores
        self.strictness = strictness
        self.genefinding_tool = genefinding_tool
        self.database_dir = database_dir
        self.debug = debug
        self.genome_files = genome_fasta_files

        if concurrent:
            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, cores, avoid_throttle=True
            )
        else:
            self.concurrent_jobs = 0
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1

    def create_antismash_cmd(self, genome_fasta_file: Path, genome_out_dir: Path):
        cmd = " ".join(
            [
                "antismash",
                f"--output-dir {genome_out_dir}",
                f"--cpu {self.cores}",
                f"--hmmdetection-strictness {self.strictness}",
                f"--genefinding-tool {self.genefinding_tool}",
                "--cb-knownclusters",
                f"{genome_fasta_file.absolute()}",
            ]
        )
        if self.database_dir is not None:
            cmd += f" --databases {self.database_dir}"

        if self.debug:
            cmd += " -v -d"
        return cmd

    def analyze_genomes(self):
        genome_out_dirs = []
        for genome_file in self.genome_files:
            genome_out_dir = self.antismash_raw_results_dir / genome_file.stem
            genome_out_dirs.append(genome_out_dir)
            genome_out_dir.mkdir(exist_ok=True, parents=True)
        commands = [
            self.create_antismash_cmd(gf, go)
            for gf, go in zip(self.genome_files, genome_out_dirs)
        ]
        if self.debug:
            for cmd in tqdm(
                commands, ascii=True, leave=True, desc="Running antiSMASH in debug mode"
            ):
                execute_command(cmd)

        else:
            _ = multiprocess_dispatch(
                "system",
                commands,
                self.concurrent_jobs,
                show_progress=True,
                description="Running antiSMASH",
            )


class smBGCParser:
    def __init__(self, out_dir: Path, cores: int):
        self.out_dir = out_dir / "smBGC"
        self.antismash_raw_result_dir = self.out_dir / "antismash_raw_results"
        self.antismash_subdirs = list(self.antismash_raw_result_dir.glob("*"))
        self.parsing_cores = cores

    def load_json(self, f: str) -> dict:
        json_str = ""
        with open(f, "r") as rf:
            for line in rf:
                line = line.rstrip()
                json_str += line
        return json.loads(json_str)

    def get_known_cluster_results(self, seq_obj) -> str:
        clusterblast_res = seq_obj["modules"]["antismash.modules.clusterblast"][
            "knowncluster"
        ]["results"]
        return clusterblast_res

    def has_identified_smbgcs(self, identified_smbgcs: list) -> bool:
        if len(identified_smbgcs) == 0:
            return False
        return True

    def parse_strain_results(self, strain_subdir: Path) -> pd.DataFrame:
        genome = strain_subdir.name
        json_f = strain_subdir / (genome + ".json")
        if not check_if_file_exists(json_f):
            return pd.DataFrame()
        json_str = self.load_json(str(json_f))
        data = []
        for seq_obj in json_str["records"]:
            identified_smbgcs = seq_obj["areas"]
            has_identified_smbgcs = self.has_identified_smbgcs(identified_smbgcs)
            if not has_identified_smbgcs:
                continue
            known_cluster_results = self.get_known_cluster_results(seq_obj)
            region_id = seq_obj["id"]
            for idx, area in enumerate(identified_smbgcs):
                # smbgcs
                region_s = area["start"]
                region_e = area["end"]
                region_p = ";".join(area["products"])
                # known cluster blast results
                if (
                    known_cluster_results is None
                    or len(known_cluster_results[idx]["ranking"]) == 0
                ):
                    num_hits = 0
                    description = "X"
                    database_total_genes = 0
                    perc_identity = 0
                else:
                    num_hits = known_cluster_results[idx]["ranking"][0][1]["hits"]
                    database_total_genes = len(
                        known_cluster_results[idx]["ranking"][0][0]["proteins"]
                    )
                    perc_identity = round(num_hits / database_total_genes * 100, 0)
                    if (
                        perc_identity > 100
                    ):  # There may be a gene duplicate in the query smBGC
                        perc_identity = 100
                    description = known_cluster_results[idx]["ranking"][0][0][
                        "description"
                    ]
                tmp_data = {
                    "Genome": genome,
                    "Region": region_id,
                    "Region_start": region_s,
                    "Region_end": region_e,
                    "Products": region_p,
                    "Most_similar_known_cluster": description,
                    "Perc._similarity": perc_identity,
                }
                data.append(tmp_data)
        return pd.DataFrame(data)

    def write_antismash_results(self):
        excel_filename = "antiSMASH_results.xlsx"
        fout = self.out_dir / excel_filename
        self.final_df.to_excel(fout, index=False)

    def gather_results(self):
        total_data = multiprocess_dispatch(
            self.parse_strain_results,
            self.antismash_subdirs,
            self.parsing_cores,
            show_progress=True,
            description="Parsing antiSMASH results",
        )
        # total_data = []
        # with ProcessPoolExecutor(self.parsing_cores) as executor:
        #     for dataframe in tqdm(
        #         executor.map(self.parse_strain_results, self.antismash_subdirs),
        #         total=len(self.antismash_subdirs),
        #         desc="Parsing antiSMASH results",
        #         ascii=True,
        #         leave=True,
        #     ):
        #         total_data.append(dataframe)
        self.final_df = pd.concat(total_data)
        self.write_antismash_results()
