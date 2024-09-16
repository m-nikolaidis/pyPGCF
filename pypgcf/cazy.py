from pathlib import Path
import re

from tqdm import tqdm
import logging

from pypgcf.utils import (
    calc_avail_dispatchers,
    multiprocess_dispatch,
    create_temporary_file,
    translate_fasta_records,
    seqrecords_to_fasta,
    dict_to_dataframe,
)


class CAZY_analyzer:
    def __init__(
        self,
        *,
        fasta_files_list: list[Path],
        database_dir: Path,
        out_dir: Path,
        input_type: str,
        available_cores: int,
        cores: int,
        concurrent: bool,
        evalue: float,
        verbose: bool,
    ):
        self.cores = cores
        self.evalue = evalue
        # The asigned column names will be the same, need to fix the homology search function
        self.database_dir = database_dir / "CAZY"
        self.out_dir = out_dir / "CAZY"
        self.search_res_dir = out_dir / "CAZY" / "CAZY_search"
        self.fasta_files = fasta_files_list
        self.verbose = verbose
        if concurrent:
            self.concurrent_jobs = calc_avail_dispatchers(
                available_cores, cores, avoid_throttle=True
            )
        else:
            self.concurrent_jobs = 0
        if self.concurrent_jobs == 0:
            self.concurrent_jobs = 1
        if input_type != "prot":  # either prot or CDS
            self.protein = False
        else:
            self.protein = True

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
        cmds = []
        if not self.search_res_dir.exists():
            self.search_res_dir.mkdir(parents=True, exist_ok=True)
        for query_file in self.fasta_files:
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
        fout = self.out_dir / "CAZY_families.xlsx"
        df = dict_to_dataframe(data)
        df.to_excel(fout)
        return None
