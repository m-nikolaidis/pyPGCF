"""
Module that contains utlity functions used throughout the software
"""

from concurrent.futures import ProcessPoolExecutor
from gzip import GzipFile
from math import floor
from os import fork, system
from pathlib import Path
from subprocess import Popen
from typing import Callable, List, Union
from zipfile import ZipFile
from tempfile import mkstemp

from Bio import Seq, SeqIO, SeqRecord
from pandas import DataFrame
from tqdm import tqdm

from pypgcf import config


def recursive_unlink(input: Path) -> None:
    """
    Delete a directory and all its contents
    If the input is a file it just deletes it
    """
    if input.is_file():
        input.unlink()
        return None

    items = input.glob("*")
    for item in items:
        if item.is_dir():
            recursive_unlink(item)
        else:
            item.unlink()
    input.rmdir()
    return None


def calc_avail_dispatchers(
    system_cores: int,
    cores_per_job: int,
    avoid_throttle: bool = True,
    usemax: bool = False,
):
    if cores_per_job > system_cores:
        raise ValueError("You have set to use more cores than the available resources")
    # Old value calculation = floor((cpu_count() - 2) / cores)
    if usemax:
        avoid_throttle = False
    val = floor(system_cores / cores_per_job)
    if val > 1 and avoid_throttle:
        return val - 1
    return val


def join_pandas_dataframes(df1: DataFrame, df2: DataFrame) -> DataFrame:
    """
    Join two pandas DataFrames on their indices using an outer join.

    If there are columns with the same name in both DataFrames, the columns from `df1` will have the suffix '_l'
    and the columns from `df2` will have the suffix '_r' appended to their names.

    Parameters:
    df1 (DataFrame): The first DataFrame to join.
    df2 (DataFrame): The second DataFrame to join.

    Returns:
    DataFrame: The joined DataFrame with combined columns from both `df1` and `df2`.

    Example:
    >>> df1 = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}, index=[0, 1, 2])
    >>> df2 = pd.DataFrame({'C': [7, 8, 9], 'D': [10, 11, 12]}, index=[0, 1, 2])
    >>> join_pandas_dataframes(df1, df2)
       A  B  C   D
    0  1  4  7  10
    1  2  5  8  11
    2  3  6  9  12

    >>> df1 = pd.DataFrame({'A': [1, 2], 'B': [3, 4]}, index=[0, 1])
    >>> df2 = pd.DataFrame({'C': [5, 6], 'D': [7, 8]}, index=[1, 2])
    >>> join_pandas_dataframes(df1, df2)
         A    B    C    D
    0  1.0  3.0  NaN  NaN
    1  2.0  4.0  5.0  7.0
    2  NaN  NaN  6.0  8.0
    """
    return df1.join(df2, rsuffix="_r", lsuffix="_l", how="outer")


def create_diamond_blastp_cmd(
    fasta_file: Path,
    database_f: Path,
    out_file: Path,
    dmnd_sensitivity: str,
    blast_evalue: float,
    blast_cores: int,
    outfmt: str,
) -> str:
    # DIAMOND case first
    dmnd_sensitivity_values = [
        "very-sensitive",
        "sensitive",
        "mid-sensitive",
        "more-sensitive",
        "ultra-sensitive",
    ]
    if dmnd_sensitivity not in dmnd_sensitivity_values:
        return ""
    cmd = f"diamond blastp --query {fasta_file} --quiet --db {database_f} --outfmt {outfmt} --out {out_file} --evalue {blast_evalue} --threads {blast_cores} --{dmnd_sensitivity}"
    return cmd


def create_blastn_cmd(
    fasta_file: Path,
    database_f: Path,
    out_file: Path,
    blast_evalue: float,
    blast_cores: int,
    outfmt: str,
) -> str:
    cmd = f"blastn -query {fasta_file} -db {database_f} -outfmt {outfmt} -out {out_file} -evalue {blast_evalue} -num_threads {blast_cores}"
    return cmd


def perform_homology_search(self): ...


def keep_best_homology_hit(
    df: DataFrame, group_col: str, score_col: str = "bitscore"
) -> DataFrame:
    keep = []
    for qseqid, qseqid_df in df.groupby(group_col):
        qseqid_df = qseqid_df.sort_values(by=score_col, ascending=False)
        keep.append(qseqid_df.index[0])
    return df.loc[keep]


def check_if_dir_is_empty(directory: Path) -> bool:
    """
    Check if a directory is empty
    """
    items = directory.glob("*")
    if not any(items):
        return True
    return False


def check_if_dir_exists(directory: Path) -> bool:
    """
    Check if a directory exists
    """
    if directory.exists():
        return True
    return False


def check_if_file_exists(file: Path) -> bool:
    """
    Check if a file exists
    """
    if file.exists():
        return True
    return False


def is_valid_antismash_strict(strictness: str):
    """
    Check if the strictness level is valid
    """

    if strictness in config.smbgc_valid_strictness:
        return True
    print(
        f"Invalid strictness level, please use one of the following: {','.join(config.smbgc_valid_strictness)}"
    )
    return False


def is_valid_genefinding_tool(tool: str) -> bool:
    """
    Check if the gene finding tool is valid
    """

    if tool in config.smbgc_genefinding_tools:
        return True
    print(
        f"Invalid gene finding tool, please use one of the following: {','.join(config.smbgc_genefinding_tools)}"
    )
    return False


def is_valid_assembly_source(source: str) -> bool:
    """
    Check if assembly source is valid.
    Needed for the download_genomes module
    """
    if source == "RefSeq" or source == "Genbank":
        return True
    return False


def multiprocess_dispatch(
    f: Union[Callable, str],
    args: list,
    num_procs: int,
    show_progress: bool,
    description: str = "",
) -> List:
    """
    Map a list of executable to a ProcessPoolExecutor
    Return: list
    """
    if isinstance(f, str):  # If it is a string then just call the system
        f = system
    with ProcessPoolExecutor(num_procs) as executor:
        if show_progress:
            results = list(
                tqdm(
                    executor.map(f, args),
                    total=len(args),
                    desc=description,
                    ascii=True,
                    leave=True,
                )
            )
        else:
            results = list(executor.map(f, args))
    return results


def execute_command(cmd: str) -> int:
    ret = system(cmd)
    return ret


def unzip_file(file: Path, method: str) -> None:
    if not check_if_file_exists(file):
        raise FileNotFoundError(f"{file} doesn't exist")
    parent_dir = file.absolute().parent
    if method != "zip" and method != "gzip":
        raise ValueError(f"Invalid compress method: {method}")
    if method == "zip":
        with ZipFile(file, "r") as zf:
            zf.extractall(parent_dir)
    if method == "gzip":
        fout = parent_dir / file.name.replace(".gz", "")
        with GzipFile(file) as rf:
            with open(fout, "wb") as wf:
                for line in rf.readlines():
                    wf.write(line)
    return None


def translate_fasta_records(file: Path) -> list:
    records = []
    for item in SeqIO.parse(file, "fasta"):
        transl_item = item.translate()
        transl_item.id = item.id
        transl_item.description = item.description
        transl_item.name = item.name
        if transl_item.seq[-1] == "*":  # Remove stop codon translation
            transl_item.seq = transl_item.seq[:-1]
        records.append(transl_item)

    return records


def seqrecords_to_fasta(records: list, fout: Path) -> None:
    SeqIO.write(records, fout, "fasta")


def create_temporary_file(dir: Union[Path, None] = None) -> Path:
    return Path(mkstemp(dir=dir)[1])


def dict_to_dataframe(d: dict) -> DataFrame:
    return DataFrame.from_dict(d, orient="index")


# TODO::     # def verify_databases_exist(database_directory: Path) -> List[bool]:
# Implement: Check if the database files from VFs, CAZY, eggNOG, antiSMASH exist
# Shouldn't throw an error if something doesn't exist, the user might want to use only specific databases
