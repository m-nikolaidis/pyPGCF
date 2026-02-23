# Default configurations for the various modules
# By M. Nikolaidis

from multiprocessing import cpu_count
from sys import platform 
import sysconfig
from pathlib import Path

DB_BASE_DIR = Path(sysconfig.get_config_var("BINLIBDEST")) / "site-packages"
TEST_DATA_DIR = Path(__file__).parent / ".." / "test_data/"

system_cores = cpu_count()
usable_cores = system_cores - 2
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
        "valid_tree_methods": ["NJ", "IQTree", "Fasttree"],
        "method": "IQTree",
        "tree_model": "TEST",
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
