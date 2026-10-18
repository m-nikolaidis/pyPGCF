import logging
from datetime import datetime
from os import linesep
from pathlib import Path
from re import search as regex_search
from typing import Union

from Bio import SeqIO

from pypgcf.config import database_dir as resolve_database_dir
from pypgcf.utils import (
    create_blastdb_cmd,
    createdmnddb_cmd,
    download_file,
    execute_command,
    get_remote_file_size,
    unzip_file,
    sanitize_vfdb_file,
)


class CAZY_installer:
    def __init__(self, *, database_dir: Union[None, Path], cores: int, debug: bool):
        self.database_dir = resolve_database_dir("cazy", database_dir)
        self.database_dir.mkdir(exist_ok=True, parents=True)
        self.cores = cores
        self.debug = debug

    def install_database(self) -> int:
        build_cmd = (
            f"dbcan_build --cpus {self.cores} --db-dir {self.database_dir} --clean"
        )
        retval = (
            execute_command(build_cmd, debug=True)
            if self.debug
            else execute_command(build_cmd)
        )

        build_date = datetime.now().strftime("%D")
        with open(self.database_dir / "log.txt", "w") as wf:
            wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return retval


class VF_installer:
    def __init__(self, *, database_dir: Union[None, Path], debug: bool):
        self.database_dir = resolve_database_dir("virulence", database_dir)
        self.database_dir.mkdir(exist_ok=True, parents=True)
        self.urls = [
            "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz",
            "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",
        ]
        self.mol_types = ["Prot", "Nucl"]
        self.debug = debug

        if self.debug:
            logging.debug(f"VFInstaller: {self.database_dir} was created")

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

    def create_vf_description_file(self, dbfout: Path) -> None:
        parser = SeqIO.parse(str(dbfout), "fasta")
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
        return None

    def install_database(self) -> None:
        for url, mol_type in zip(self.urls, self.mol_types):
            filename = url.split("/")[-1]
            database_fasta = self.database_dir / filename
            dl_size = download_file(url, database_fasta)
            exp_size = get_remote_file_size(url)
            if exp_size == 0:
                logging.warning("The remote file size could not be established")
            if dl_size != exp_size:
                logging.warning(
                    f"Expected and downloaded file sizes differ {dl_size}/{exp_size}"
                )
            #     raise ConnectionError(
            #         "Virulence database was not downloaded, please retry"
            #     )

            unzip_file(database_fasta, "gzip")
            dbfout = self.database_dir / filename.replace(".gz", "")
            sanitize_vfdb_file(dbfout)
            self.create_vf_description_file(dbfout)
            fout = self.database_dir / dbfout.stem
            if mol_type == "Prot":
                cmd = createdmnddb_cmd(fin=dbfout, fout=fout, debug=self.debug)
            else:
                cmd = create_blastdb_cmd(
                    fin=dbfout, fout=fout, input_type="nucl", debug=self.debug
                )

            res = execute_command(cmd, debug=self.debug)
            if res != 0:
                logging.error(
                    "Something went wrong during the build process of VF database"
                )
        return None


class AMR_installer:
    def __init__(self, *, database_dir: Union[None, Path], debug: bool):
        self.database_dir = resolve_database_dir("amr", database_dir)
        self.database_dir.mkdir(exist_ok=True, parents=True)
        self.debug = debug

    def install_database(self) -> None:
        """
        Update the amrfinder database from NCBI
        Return: None
        """
        cmd = ["amrfinder_update", "-d", str(self.database_dir)]
        if not self.debug:
            cmd.append("--quiet")
        if execute_command(cmd, debug=self.debug) != 0:
            logging.error("amrfinder_update failed to install the AMR database")
        return None


class SMBGC_installer:
    def __init__(self, *, database_dir: Union[None, Path], debug: bool = False):
        self.database_dir = resolve_database_dir("smbgc", database_dir)
        if self.database_dir is not None:
            self.database_dir.mkdir(exist_ok=True, parents=True)
        self.debug = debug

    def install_database(self) -> None:
        logging.debug("Downloading databases of antiSMASH")
        cmd = ["download-antismash-databases"]
        if self.database_dir is not None:
            cmd += ["--database-dir", str(self.database_dir)]
        ret = execute_command(cmd, debug=self.debug)
        if ret != 0:
            raise RuntimeError(
                "Something went wrong with the download of antiSMASH database"
            )
        logging.debug("Installed antiSMASH database successfully")
        if self.database_dir is not None:
            build_date = datetime.now().strftime("%D")
            with open(self.database_dir / "log.txt", "w") as wf:
                wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return None


class EGGNOG_installer:
    def __init__(self, *, database_dir: Union[None, Path] = None, debug: bool = False):
        self.debug = debug
        self.database_dir = resolve_database_dir("eggnog", database_dir)
        self.database_dir.mkdir(exist_ok=True, parents=True)

    def install_database(self) -> None:
        cmd = ["download_eggnog_data.py", "--data_dir", str(self.database_dir), "-y"]
        if not self.debug:
            cmd.append("-q")
        ret = execute_command(cmd, debug=self.debug)
        if ret != 0:
            raise RuntimeError(
                "Something went wrong with the download of eggNOG database"
            )
        logging.debug("Installed eggNOG database successfully")
        build_date = datetime.now().strftime("%D")
        with open(self.database_dir / "log.txt", "w") as wf:
            wf.write(f"# Database downloaded and built: {build_date}{linesep}")

        return None
