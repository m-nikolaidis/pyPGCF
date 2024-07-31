from multiprocessing import cpu_count
from pathlib import Path
import re

from tqdm import tqdm
import logging

from pypgcf.utils import (
    calc_avail_dispatchers,
    execute_command,
    multiprocess_dispatch,
    create_temporary_file,
    translate_fasta_records,
    seqrecords_to_fasta,
    dict_to_dataframe,
)


class CAZY_builder:
    def __init__(self, database_dir: Path, cores: int, verbose: bool):
        self.database_dir = database_dir / "CAZY"
        self.cores = cores

    def setup(self) -> int:
        build_cmd = (
            f"dbcan_build --cores {self.cores} --db-dir {self.database_dir} --clean"
        )
        retval = execute_command(build_cmd)
        return retval

    def validate_if_database_exists(self): ...


class CAZY_analyzer:
    def __init__(
        self,
        cores: int,
        evalue: float,
        database_dir: Path,
        results_dir: Path,
        fasta_dir: Path,
        dmnd_sensitivity: str,
        input_type: str,
        verbose: bool,
    ):
        self.cores = cores
        self.evalue = evalue
        # The asigned column names will be the same, need to fix the homology search function
        self.database_dir = database_dir / "CAZY"
        self.results_dir = results_dir / "CAZY"
        self.search_res_dir = results_dir / "CAZY" / "CAZY_search"
        self.fasta_dir = fasta_dir
        self.verbose = verbose
        self.concurrent_jobs = calc_avail_dispatchers(
            cpu_count(), cores, avoid_throttle=True
        )
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1
        self.protein = True
        if input_type != "prot":  # either prot or nucl
            self.protein = False

        self.HMM_COORDS_RE = re.compile(r"\(\d+-\d+\)")
        self.FAMILY_MAPPING = {
            "AA": "Auxiliary activities",
            "CBM": "Carbohydrate-Binding Module Family",
            "CE": "Carbohydrate Esterase Family",
            "GH": "Glycoside Hydrolase Family",
            "GT": "Glycosyltransferase Family",
            "PL": "Polysaccharide Lyase Family",
        }

    def execute_cazy_search(self) -> None:
        query_files = list(self.fasta_dir.glob("*"))
        if len(query_files) == 0:
            raise FileNotFoundError
        cmds = []
        if not self.search_res_dir.exists():
            self.search_res_dir.mkdir(parents=True)
        for query_file in query_files:
            outdir = self.search_res_dir / query_file.stem
            if not self.protein:
                protein_f = create_temporary_file()
                records = translate_fasta_records(query_file)
                seqrecords_to_fasta(records, protein_f)
            else:
                protein_f = query_file
            cmd = " ".join(
                [
                    f"run_dbcan {protein_f} protein",
                    f"--out_dir {outdir}",
                    f"--db_dir {self.database_dir}",
                    f"--dbcan_thread {self.cores}",
                    f"--dia_cpu {self.cores}",
                    f"--dia_eval {self.evalue}",
                    f"--hmm_cpu {self.cores}",
                    f"--hmm_eval {self.evalue}",
                ]
            )
            cmds.append(cmd)

        _ = multiprocess_dispatch(
            "system",
            cmds,
            self.concurrent_jobs,
            show_progress=True,
            description="Scanning for CAZYmes",
        )
        return None

    def _clean_hmm_output(self, cazy: str) -> str:
        cazy_domains = cazy.split("+")  # Split the various results
        for idx, cazy_domain in enumerate(cazy_domains):
            # Remove coordinates
            cazy_domain = re.sub(self.HMM_COORDS_RE, "", cazy_domain)
            # Remove subfamily
            cazy_domain = cazy_domain.split("_")[0]
            cazy_domain = "".join([s for s in cazy_domain if not s.isdigit()])
            cazy_domains[idx] = cazy_domain
        cazy_domains = set(cazy_domains)
        return ";".join(cazy_domains)

    def _clean_dmnd_dbcansub_output(self, cazy: str) -> str:
        cazy_domains = cazy.split("+")  # Split the various results
        for idx, cazy_domain in enumerate(cazy_domains):
            # Remove subfamily
            cazy_domain = cazy_domain.split("_")[0]
            cazy_domain = "".join([s for s in cazy_domain if not s.isdigit()])
            cazy_domains[idx] = cazy_domain
        cazy_domains = set(cazy_domains)
        return ";".join(cazy_domains)

    def parse_results(self) -> None:
        tool_cutoff = 3  # Number of tools for CAZyome annotation (3 == all)
        directories = list(self.search_res_dir.glob("*"))
        if len(directories) == 0:
            logging.info("No output for CAZY search was created")
        data = {}
        for directory in tqdm(directories, ascii=True):
            genome = directory.name
            data[genome] = {f: 0 for f in self.FAMILY_MAPPING.keys()}
            f = directory / "overview.txt"
            with open(f, "r") as rf:
                next(rf)  # Advance headers to remove the need of the check
                for line in rf:
                    line = line.rstrip()
                    prot, ec, _, _, dmnd_r, num_tools = line.split("\t")
                    if int(num_tools) < tool_cutoff:
                        continue
                    # hmmer_r = self._clean_hmm_output(hmmer_r)
                    dmnd_r = self._clean_dmnd_dbcansub_output(dmnd_r)
                    # dbcansub_r = self._clean_dmnd_dbcansub_output(dbcansub_r)
                    for family in dmnd_r.split(";"):
                        data[genome][family] += 1
        fout = self.results_dir / "CAZY_families.xlsx"
        df = dict_to_dataframe(data)
        df.to_excel(fout)
        return None
