"""
Perform necessary checks during environment setup
"""
from pathlib import Path

from pypgcf import config


def check_if_dir_is_empty(directory: Path):
    """
    Check if a directory is empty
    """
    items = directory.glob("*")
    if not any(items):
        return True
    return False


def check_if_dir_exists(directory: Path):
    """
    Check if a directory exists
    """
    if directory.exists():
        return True
    return False


def check_if_file_exists(file: Path):
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


def is_valid_genefinding_tool(tool: str):
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
    if source == "RefSeq" or source == "Genbank":
        return True
    return False
