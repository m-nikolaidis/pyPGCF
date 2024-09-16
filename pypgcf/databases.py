from datetime import datetime
from os import linesep
from pathlib import Path
from re import search as regex_search
from typing import Union
import logging

from Bio import SeqIO

from pypgcf.utils import (
    execute_command,
    download_file,
    get_remote_file_size,
    unzip_file,
    create_blastdb_cmd,
)
from pypgcf.config import DB_BASE_DIR


class CAZY_installer:
    def __init__(self, *, database_dir: Union[None, Path], cores: int, debug: bool):
        if database_dir is None:
            self.database_dir = DB_BASE_DIR / "CAZY"
        else:
            self.database_dir = database_dir / "CAZY"
        self.database_dir.mkdir(exist_ok=True, parents=True)
        self.cores = cores
        self.debug = debug

    def install_database(self) -> int:
        build_cmd = (
            f"dbcan_build --cores {self.cores} --db-dir {self.database_dir} --clean"
        )
        retval = execute_command(build_cmd)

        build_date = datetime.now().strftime("%D")
        with open(self.database_dir / "log.txt", "w") as wf:
            wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return retval


class VF_installer:
    def __init__(self, *, database_dir: Union[None, Path], debug: bool):
        if database_dir is None:
            self.database_dir = DB_BASE_DIR / "Virulence"
        else:
            self.database_dir = database_dir / "Virulence"
        self.database_dir.mkdir(exist_ok=True, parents=True)
        self.urls = [
            "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz",
            "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",
        ]
        self.mol_types = ["Prot", "Nucl"]
        if debug:
            print(f"{self.database_dir} was created")

    def _split_VF_desc(self, string: str) -> tuple:
        category_pattern = r"\[.+\) - (.+) \(.+\] \[.+$"
        origin_pattern = r"\] \[(.+)\]$"
        description_pattern = r"(.+) \[.+ \["

        match1 = regex_search(category_pattern, string)
        category = match1.group(1) if match1 else None

        match2 = regex_search(origin_pattern, string)
        origin = match2.group(1) if match2 else None

        match3 = regex_search(description_pattern, string)
        desc = match3.group(1) if match3 else None
        return category, origin, desc

    def install_database(self) -> None:
        # Get the fasta annotations
        for url, mol_type in zip(self.urls, self.mol_types):
            filename = url.split("/")[-1]
            database_fasta = self.database_dir / filename
            dl_size = download_file(url, database_fasta)
            exp_size = get_remote_file_size(url)
            if dl_size != exp_size:
                raise ConnectionError(
                    "Virulence database was not downloaded, please retry"
                )

            unzip_file(database_fasta, "gzip")
            dbfout = self.database_dir / filename.split(".")[0]
            parser = SeqIO.parse(database_fasta, "fasta")
            vf_desc = {}
            for record in parser:
                vf = record.id
                desc = record.description
                category, origin, desc = self._split_VF_desc(desc)
                vf_desc[vf] = [category, origin, desc]
            desc_fout = self.database_dir / "vfdb_desc.tsv"
            with open(desc_fout, "w") as wf:
                build_date = datetime.now().strftime("%D")
                wf.write(f"# Database downloaded and built: {build_date}{linesep}")
                for vf, items in vf_desc.items():
                    str_to_write = vf + "\t" + "\t".join(items) + linesep
                    wf.write(str_to_write)

            if mol_type == "Prot":
                cmd = f"diamond makedb --db {dbfout} --in {database_fasta} --threads 4"
                if not self.debug:
                    cmd += " --quiet"
            else:
                cmd = f"makeblastdb -dbtype nucl -in {database_fasta} -out {dbfout}"
                if not self.debug:
                    cmd += " > /dev/null 2> /dev/null"

            res = execute_command(cmd)
            if res != 0:
                logging.error(
                    "Something went wrong during the build process of VF database"
                )
            return None


class AMR_installer:
    def __init__(self, *, database_dir: Union[None, Path], debug: bool):
        if database_dir is None:
            self.database_dir = DB_BASE_DIR / "AMR"
        else:
            self.database_dir = database_dir / "AMR"
        self.database_dir.mkdir(exist_ok=True, parents=True)

        self.url = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/"
        self.files = [
            "AMR.LIB",
            "AMRProt",
            "AMRProt-mutation.tab",
            "AMRProt-suppress",
            "AMRProt-susceptible.tab",
            "AMR_CDS",
            "AMR_DNA-Acinetobacter_baumannii",
            "AMR_DNA-Acinetobacter_baumannii.tab",
            "AMR_DNA-Campylobacter",
            "AMR_DNA-Campylobacter.tab",
            "AMR_DNA-Clostridioides_difficile",
            "AMR_DNA-Clostridioides_difficile.tab",
            "AMR_DNA-Enterococcus_faecalis",
            "AMR_DNA-Enterococcus_faecalis.tab",
            "AMR_DNA-Enterococcus_faecium",
            "AMR_DNA-Enterococcus_faecium.tab",
            "AMR_DNA-Escherichia",
            "AMR_DNA-Escherichia.tab",
            "AMR_DNA-Klebsiella_oxytoca",
            "AMR_DNA-Klebsiella_oxytoca.tab",
            "AMR_DNA-Klebsiella_pneumoniae",
            "AMR_DNA-Klebsiella_pneumoniae.tab",
            "AMR_DNA-Neisseria_gonorrhoeae",
            "AMR_DNA-Neisseria_gonorrhoeae.tab",
            "AMR_DNA-Salmonella",
            "AMR_DNA-Salmonella.tab",
            "AMR_DNA-Staphylococcus_aureus",
            "AMR_DNA-Staphylococcus_aureus.tab",
            "AMR_DNA-Streptococcus_pneumoniae",
            "AMR_DNA-Streptococcus_pneumoniae.tab",
            "ReferenceGeneCatalog.txt",
            "ReferenceGeneHierarchy.txt",
            "amr_targets.fa",
            "changelog.txt",
            "changes.txt",
            "database_format_version.txt",
            "fam.tab",
            "mapgenelist.txt",
            "taxgroup.tab",
            "version.txt",
            "AMR.LIB",
            "AMRProt",
            "AMRProt-mutation.tab",
            "AMRProt-suppress",
            "AMRProt-susceptible.tab",
            "AMR_CDS",
            "database_format_version.txt",
            "fam.tab",
            "taxgroup.tab",
        ]
        self.debug = debug

    def _index_database(self) -> None:
        for file in self.files:
            if (
                file.endswith(".tab")
                or file.endswith(".txt")
                or file.endswith(".LIB")
                or file.endswith("suppress")
            ):
                continue
            filetype = "nucl"
            if "AMRProt" in file:
                filetype = "prot"
            outfile = self.database_dir / file
            cmd = create_blastdb_cmd(outfile, filetype)
            _ = execute_command(cmd)

        return None

    def install_database(self) -> None:
        for file in self.files:
            url = self.url + "/" + file
            fout = self.database_file / file
            obs_size = download_file(url, fout)
            exp_size = get_remote_file_size(url)
            if obs_size != exp_size:
                logging.error(
                    f"Something went wrong with the download of {file} (AMR database)"
                )
        self._index_database()
        build_date = datetime.now().strftime("%D")
        with open(self.database_dir / "log.txt", "w") as wf:
            wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return None


class SMBGC_installer:
    def __init__(self, *, database_dir: Union[None, Path], debug: bool = False):
        if database_dir is not None:
            self.database_dir = database_dir / "smBGC"
        else:
            self.database_dir = None
        self.database_dir.mkdir(exist_ok=True, parents=True)
        self.debug = debug

    def install_database(self) -> None:
        print("Downloading databases of antiSMASH")
        cmd = "download-antismash-databases"
        if self.database_dir is not None:
            cmd += f" --database_dir {self.database_dir}"
        ret = execute_command(cmd)
        if self.debug:
            if ret == 0:
                print("Installed antiSMASH database successfully")
            else:
                raise RuntimeError(
                    "Something went wrong with the download of antiSMASH database"
                )
        build_date = datetime.now().strftime("%D")
        with open(self.database_dir / "log.txt", "w") as wf:
            wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return None


class EGGNOG_installer:
    def __init__(self, *, database_dir: Union[None, Path] = None, debug: bool = False):
        self.debug = debug
        if database_dir is None:
            self.database_dir = DB_BASE_DIR / "site-packages" / "data"
        else:
            self.database_dir = database_dir / "eggNOG"
        self.database_dir.mkdir(exist_ok=True, parents=True)

    def install_database(self) -> None:
        cmd = f"download_eggnog_data.py --data_dir {self.database_dir} -y"
        if not self.debug:
            cmd += " -q"
        ret = execute_command(cmd)
        if self.debug:
            if ret == 0:
                print("Installed eggNOG database successfully")
            else:
                raise RuntimeError(
                    "Something went wrong with the download of eggNOG database"
                )
        build_date = datetime.now().strftime("%D")
        with open(self.database_dir / "log.txt", "w") as wf:
            wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return None
