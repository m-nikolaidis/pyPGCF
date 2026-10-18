# Default configurations for the various modules
# By M. Nikolaidis

import logging
from multiprocessing import cpu_count
import os
from sys import platform 
import sysconfig
from pathlib import Path
from typing import Optional, Union

DB_BASE_DIR = Path(sysconfig.get_config_var("BINLIBDEST")) / "site-packages"
DB_DIR_ENV = "PYPGCF_DB_DIR"
# DIAMOND sensitivity modes, passed to DIAMOND as --<mode>
DMND_SENSITIVITIES = [
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
]
DB_SUBDIRS = {
    "cazy": "CAZY",
    "virulence": "Virulence",
    "amr": "AMR",
    "smbgc": "smBGC",
    "eggnog": "eggNOG",
}


def database_dir(
    module: str, cli_value: Union[str, Path, None] = None
) -> Optional[Path]:
    """
    Database directory for a module. The --db option wins, then the
    PYPGCF_DB_DIR environment variable, then the default location.
    """
    root = cli_value or os.environ.get(DB_DIR_ENV)
    if root:
        path = Path(root) / DB_SUBDIRS[module]
    elif module == "smbgc":
        path = None  # antiSMASH uses its own default database location
    elif module == "eggnog":
        path = DB_BASE_DIR / "site-packages" / "data"  # used by earlier releases
    else:
        path = DB_BASE_DIR / DB_SUBDIRS[module]
    logging.info(f"{module} database directory: {path}")
    return path
TEST_DATA_DIR = Path(__file__).parent / ".." / "test_data/"

system_cores = cpu_count()


def usable_core_count(system_cores: int) -> int:
    """Leave two cores free for the system, but always use at least one."""
    return max(1, system_cores - 2)


usable_cores = usable_core_count(system_cores)
default_no_concurrent = False
default_debug = False

software_defaults = {
    "system": {
        "platform": platform,
        "system_cores": system_cores,
        "usable_cores": usable_cores,
    },
    "general": {
        "input_type": "prot",
        "valid_input_types": ["prot", "nucl", "cds"],
    },
    "species_demarcation": {
        "cores": usable_cores,
        "kmer": 16,
        "fraglen": 3000,
        "minfrac": 0.2,
        "inflation": 2,
        "prefix":"C",
        "debug": default_debug,
    },
    "orthologues": {
        "cores": min(6, usable_cores),
        "evalue": 1e-5,
        "dmnd_sensitivity": "very-sensitive",
        "no_concurrent": default_no_concurrent,
        "no_filter_orthologues": False,
        "debug": default_debug,
    },
    "core": {
        "core_perc": 100,
        "debug": default_debug,
    },
    "phylogenomic": {
        "cores": usable_cores,
        "valid_tree_methods": ["NJ", "IQTree"],
        "method": "IQTree",
        "tree_model": "TEST",
        "fast": False,
        "input_type": "prot",
        "valid_input_types": ["prot", "cds"],
        "no_keep_fasta": False,
        "debug": default_debug,
    },
    "eggnog": {
        "cores": usable_cores,
        "pident": 40,
        "scov": 20,
        "qcov": 20,
        "debug": default_debug,
        "input_type": "prot",
        "valid_input_types": ["prot", "cds"],
    },
    "smbgc": {
        "cores": min(6, usable_cores),
        "strictness": "strict",
        "genefinding_tool": "prodigal",
        "valid_strictness": ["loose", "relaxed", "strict"],
        "genefinding_tools": ["prodigal"],
        "no_concurrent": default_no_concurrent,
        "debug": default_debug,
    },
    "virulence": {
        "blast_cores": min(6, usable_cores),
        "blast_evalue": 1e-5,
        "dmnd_sensitivity": "very-sensitive",
        "no_concurrent": default_no_concurrent,
        "debug": default_debug,
        "input_type": "prot",
        "valid_input_types": ["prot", "cds"],
    },
    "cazy": {
        "cores": min(6, usable_cores),
        "evalue": 1e-5,
        "no_concurrent": default_no_concurrent,
        "debug": default_debug,
        "input_type": "prot",
        "valid_input_types": ["prot", "cds"],
    },
    "amr": {
        "cores": min(6, usable_cores),
        "evalue": 1e-5,
        "no_concurrent": default_no_concurrent,
        "debug": default_debug,
        "input_type": "prot",
        "valid_input_types": ["prot", "cds"],
    },
    "database": {
        "db": "all",
        "db_dir": None,
        "cores": min(4, usable_cores),
        "debug": default_debug,
    },
    "download": {
        "valid_sources": ["RefSeq", "GenBank"],
        "assembly_level": "chromosome,complete",
    },
}

# Conda package (conda-forge or bioconda) that provides each external binary
TOOL_PACKAGES = {
    "fastANI": "fastani",
    "mcxload": "mcl",
    "mcl": "mcl",
    "mcxdump": "mcl",
    "diamond": "diamond",
    "blastn": "blast",
    "makeblastdb": "blast",
    "muscle": "muscle",
    "iqtree": "iqtree",
    "emapper.py": "eggnog-mapper",
    "antismash": "antismash",
    "amrfinder": "ncbi-amrfinderplus",
    "run_dbcan": "dbcan",
    "datasets": "ncbi-datasets-cli",
}

# External binaries each module calls
MODULE_TOOLS = {
    "species_demarcation": ["fastANI", "mcxload", "mcl", "mcxdump"],
    "orthologues": ["diamond"],
    "orthologues_cds": ["blastn", "makeblastdb"],
    "core": [],
    "phylogenomic": ["muscle", "iqtree"],
    "eggnog": ["emapper.py"],
    "smbgc": ["antismash"],
    "amr": ["amrfinder"],
    "cazy": ["run_dbcan"],
    "virulence": ["diamond"],
    "download": ["datasets"],
}
