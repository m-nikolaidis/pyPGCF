import csv
from datetime import datetime
from pathlib import Path
from typing import Generator, List, Union

from pandas import read_csv

from pypgcf.utils import check_if_file_exists, execute_command, dict_to_dataframe


class SpeciesDemarcator:
    def __init__(
        self,
        *,
        in_dir: Path,
        out_dir: Path,
        fastani_cores: int,
        kmer: int,
        fraglen: int,
        minfrac: float,
        inflation: float,
        mcl_cores: int,
        debug: bool = False,
    ):
        self.in_dir = in_dir
        self.out_dir = out_dir / "Species_demarcation"
        self.fastani_cores = fastani_cores
        self.kmer = kmer
        self.fraglen = fraglen
        self.minfrac = minfrac
        self.inflation = inflation
        self.mcl_cores = mcl_cores
        self.debug = debug

    def create_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)

    def create_input_for_fastani(
        self, files_for_fastani: Union[List, Generator], tmp_file_for_fastani: Path
    ) -> None:
        print("Preparing FastANI input")
        with open(str(tmp_file_for_fastani), "w") as f:
            for file in files_for_fastani:
                f.write(str(file) + "\n")
        return None

    def perform_fastani(self, org_list: Path, fout: Path) -> None:
        print(f"Performing FastANI: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        cmd = "fastANI --ql {} --rl {} -t {} -k {} --fragLen {} --minFraction {} -o {}".format(
            org_list,
            org_list,
            self.fastani_cores,
            self.kmer,
            self.fraglen,
            self.minfrac,
            fout,
        )
        if not self.debug:
            cmd += " > /dev/null 2>&1"
        ret = execute_command(cmd)
        if ret != 0:
            raise RuntimeError("fastANI command was not successful")
        return None

    def prepare_input_for_mcl(self, input_file: Path) -> Path:
        """
        This script prepares the input file for fastANI.
        """
        headers = ["query", "target", "ANI", "query_length", "target_length"]
        df = read_csv(input_file, sep="\t", index_col=0, names=headers)
        df = df[df["ANI"] >= 95]
        df = df.drop(columns=["query_length", "target_length"])
        fout = input_file.parent / "fastANI_for_mcl.txt"
        df.to_csv(fout, sep="\t", header=False)
        return fout

    def clean_mcl(self, outdir: Path) -> None:
        to_remove = [
            outdir / "fastANI_for_mcl.txt",
            outdir / "fastANI_mcx_mtrx.txt",
            outdir / "fastANI_annot.tab",
            outdir / "fastANI_mcl_out.txt",
            outdir / "FastANI_input.txt",
        ]
        for f in to_remove:
            f.unlink(missing_ok=True)

        to_rename = outdir / "fastANI_mcx_dump.txt"
        new_name = outdir / "fastANI_clusters.tsv"
        to_rename.rename(new_name)

    def run_mcl(self, fastani_for_mcl: Path) -> Path:
        print(
            f"Running MCL clustering: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        outdir = fastani_for_mcl.parent
        mcxload_cmd = f"mcxload -abc {fastani_for_mcl} -o {outdir}/fastANI_mcx_mtrx.txt -write-tab {outdir}/fastANI_annot.tab"
        mcl_cmd = f"mcl {outdir}/fastANI_mcx_mtrx.txt -te {self.mcl_cores} -I {self.inflation} -o {outdir}/fastANI_mcl_out.txt"
        mcxdump_cmd = f"mcxdump -icl {outdir}/fastANI_mcl_out.txt -tabr {outdir}/fastANI_annot.tab -o {outdir}/fastANI_mcx_dump.txt"
        cmds = [mcxload_cmd, mcl_cmd, mcxdump_cmd]

        for cmd in cmds:
            if not self.debug:
                cmd += " > /dev/null 2>&1"
            ret = execute_command(cmd)
            if ret != 0:
                raise RuntimeError("Something went wrong with MCL")
        self.clean_mcl(outdir)
        return outdir / "fastANI_clusters.tsv"

    def parse_mcx_output(self, fastani_from_mcl: Path) -> None:
        print(f"Parsing MCL output: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        outdir = fastani_from_mcl.parent
        results = {}
        clust_num = 0
        for lines in csv.reader(open(str(fastani_from_mcl), "r"), delimiter="\t"):
            for line in lines:
                results[line] = clust_num
            clust_num += 1
        df = dict_to_dataframe(results)
        df.columns = ["ClustNum"]
        df["FastANI_species"] = df["ClustNum"].apply(lambda x: "C" + str(x))
        df = df.drop("ClustNum", axis=1)
        df.index.name = "Genome"
        df.index = [idx.split("/")[-1] for idx in df.index]
        df.index = [".".join(idx.split(".")[:-1]) for idx in df.index]
        fastani_from_mcl.unlink()
        fout = outdir / "FastANI_species_clusters.xlsx"
        df.to_excel(fout)

    def assign_species(self):
        # Check if input is file or directory
        self.create_directories()

        files_for_fastani = self.in_dir.glob("*")
        tmp_file_for_fastani = self.out_dir / "FastANI_input.txt"
        self.create_input_for_fastani(files_for_fastani, tmp_file_for_fastani)
        if not check_if_file_exists(tmp_file_for_fastani):
            raise FileNotFoundError(f"{tmp_file_for_fastani} was not created")

        fastani_out = self.out_dir / "FastANI.tsv"
        self.perform_fastani(tmp_file_for_fastani, fastani_out)
        if not check_if_file_exists(fastani_out):
            raise FileNotFoundError(f"{fastani_out} was not created")

        fastani_for_mcl = self.prepare_input_for_mcl(fastani_out)
        if not check_if_file_exists(fastani_for_mcl):
            raise FileNotFoundError(f"{fastani_for_mcl} was not created")

        fastani_from_mcl = self.run_mcl(fastani_for_mcl)
        if not check_if_file_exists(fastani_from_mcl):
            raise FileNotFoundError(f"{fastani_from_mcl} was not created")

        self.parse_mcx_output(fastani_from_mcl)
        print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
