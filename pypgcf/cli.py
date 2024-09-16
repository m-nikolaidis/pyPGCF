"""
Author: Marios Nikolaidis
Git: https://github.com/m-nikolaidis
email: marionik23@gmail.com
"""

import argparse
import logging
from datetime import datetime
from pathlib import Path

from pandas import read_excel
from tqdm import tqdm

from pypgcf import utils
from pypgcf.amr import AMR_analyzer
from pypgcf.cazy import CAZY_analyzer
from pypgcf.config import software_defaults, DB_BASE_DIR
from pypgcf.core import Core_identifier
from pypgcf.databases import (
    AMR_installer,
    CAZY_installer,
    EGGNOG_installer,
    SMBGC_installer,
    VF_installer,
)
from pypgcf.download_genomes import GenomeDownloader
from pypgcf.eggnog import EggNOGParser, EggNOGRunner
from pypgcf.orthologues import Orthologues_identifier
from pypgcf.phylogenomic import Phylogenomic
from pypgcf.smbgc import smBGCLocalRunner, smBGCParser
from pypgcf.species_demarcation import SpeciesDemarcator
from pypgcf.virulence import VF_analyzer
from pypgcf.workflow import WorkflowRunner

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# TODO: For CLI, in order to install the databases the user should invoke
# pyPGCF databases --install all --databse_dir ./ for all
# or
# pyPGCF databases --install cazy,smbgc for specific --databse_dir ./
# if --all is passed then the other options should be ignored


def validate_directory(path: Path) -> bool:
    if not utils.check_if_dir_exists(path):
        logging.error(f"{path} does not exist")
        return False
    if utils.check_if_dir_is_empty(path):
        logging.error(f"{path} is empty")
        return False
    return True


def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pyPGCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="pyPGCF: PhyloGenomic, Core and Fingerprint analysis software",
        allow_abbrev=False,
    )
    subparsers = parser.add_subparsers(
        help="The various modules of the program", dest="module", required=True
    )

    # Add subparsers for each module
    add_species_demarcation_subparser(subparsers)
    add_orthologues_subparser(subparsers)
    add_core_subparser(subparsers)
    add_phylogenomic_subparser(subparsers)
    add_eggnog_subparser(subparsers)
    add_smbgc_subparser(subparsers)
    add_download_subparser(subparsers)
    add_virulence_subparser(subparsers)
    add_cazy_subparser(subparsers)
    add_amr_subparser(subparsers)
    add_workflow_subparser(subparsers)
    add_databases_subparser(subparsers)

    return parser


def add_species_demarcation_subparser(subparsers):
    species_demarcation = subparsers.add_parser(
        "species_demarcation",
        help="species_demarcation module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Assign the input genomes to species clusters using FastANI and MCL",
    )
    species_demarcation_basic = species_demarcation.add_argument_group("Basic options")
    species_demarcation_basic.add_argument(
        "-in", metavar="input", help="Genome fasta directory", required=True
    )
    species_demarcation_basic.add_argument(
        "-o", metavar="out", help="Output directory", required=True
    )
    species_demarcation_basic.add_argument(
        "--debug", help="Used for debugging purposes", action="store_true"
    )

    species_demarcation_fastani = species_demarcation.add_argument_group(
        "FastANI options"
    )
    species_demarcation_fastani.add_argument(
        "--fastani_cores",
        help="Number of cores for FastANI",
        default=software_defaults["species_demarcation"]["cores"],
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--kmer",
        metavar="N",
        help="kmer size (<= 16)",
        default=software_defaults["species_demarcation"]["kmer"],
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--fraglen",
        metavar="N",
        help="FastANI Fragment length",
        default=software_defaults["species_demarcation"]["fraglen"],
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--minfraction",
        metavar="N",
        help="FastANI Minimum fraction",
        default=software_defaults["species_demarcation"]["minfrac"],
        type=float,
    )

    species_demarcation_mcl = species_demarcation.add_argument_group("MCL options")
    species_demarcation_mcl.add_argument(
        "--inflation",
        help="Inflation parameter for MCL",
        default=software_defaults["species_demarcation"]["inflation"],
    )
    species_demarcation_mcl.add_argument(
        "--mcl_cores",
        help="Number of cores for MCL",
        default=software_defaults["species_demarcation"]["cores"],
        type=int,
    )


def add_orthologues_subparser(subparsers):
    orthologues = subparsers.add_parser(
        "orthologues",
        help="orthologues module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Calculate the orthologues of the dataset using one reference strain. By default this module will run multiple searches in parallel to accelerate the analysis",
    )
    orthologues_basic = orthologues.add_argument_group("Basic options")
    orthologues_basic.add_argument(
        "-in", metavar="in", help="Fasta directory", required=True
    )
    orthologues_basic.add_argument(
        "-out", metavar="out", help="Output directory", required=True
    )
    ref_exclusive = orthologues_basic.add_mutually_exclusive_group(required=True)
    ref_exclusive.add_argument("-ref", help="Reference strain")
    ref_exclusive.add_argument("-ref_list", help="List of reference strains to use")
    orthologues_basic.add_argument(
        "--input_type",
        help="Type of input [proteins (DIAMOND) or CDS (BLASTN)]",
        default=software_defaults["general"]["input_type"],
    )
    orthologues_blast = orthologues.add_argument_group("DIAMOND/BLASTN options")
    orthologues_blast.add_argument(
        "--cores",
        help="Number of BLAST cores",
        default=software_defaults["orthologues"]["cores"],
        type=int,
    )
    orthologues_blast.add_argument(
        "--evalue",
        help="E-value cut-off",
        default=software_defaults["orthologues"]["evalue"],
    )
    orthologues_blast.add_argument(
        "--dmnd_sensitivity",
        help="Sensitivity settings for DIAMOND (i.e. very_sensitive). Not used when input type is CDS/nucl",
        default=software_defaults["orthologues"]["dmnd_sensitivity"],
    )
    orthologues_blast.add_argument(
        "--no_concurrent",
        help="Do not run multiple searches in parallel",
        action="store_true",
    )
    orthologues_blast.add_argument(
        "--no_filter_orthologues", help="Do not filter orthologues", action="store_true"
    )


def add_core_subparser(subparsers):
    core = subparsers.add_parser(
        "core",
        help="core module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    core_input_group = core.add_mutually_exclusive_group(required=True)
    core_input_group.add_argument(
        "-in", help="Input orthology matrix (from 'orthologues' module)"
    )
    core_input_group.add_argument(
        "-ref_list",
        help="List of orthology matrices to use for multiple analyses",
    )
    core.add_argument("-out", help="Output directory", required=True)
    core.add_argument(
        "--species",
        help="Input species assignment matrix (excel format; from 'species_demarcation' module)",
    )
    core.add_argument(
        "--core_perc",
        help="Percent presence of a protein/gene in a cluster to be considered core",
        default=software_defaults["core"]["core_perc"],
    )


def add_phylogenomic_subparser(subparsers):
    phylogenomic = subparsers.add_parser(
        "phylogenomic",
        help="phylogenomic module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    phylogenomic.add_argument("-fasta_dir", help="Input fasta directory", required=True)
    phylogenomic.add_argument(
        "-in", help="Input orthology matrix (from 'orthologues' module)", required=True
    )
    phylogenomic.add_argument("-out", help="Output directory", required=True)
    phylogenomic.add_argument(
        "--input_type",
        help="Input type [prot or cds]",
        default=software_defaults["phylogenomic"]["input_type"],
        choices=software_defaults["phylogenomic"]["valid_input_types"],
    )
    phylogenomic.add_argument(
        "--cores",
        help="Number of cores",
        default=software_defaults["phylogenomic"]["cores"],
        type=int,
    )
    phylogenomic.add_argument(
        "--method",
        help="Tree method (IQTree, Fasttree, NJ)",
        default=software_defaults["phylogenomic"]["method"],
        choices=software_defaults["phylogenomic"]["valid_tree_methods"],
    )

    phylogenomic.add_argument(
        "--debug",
        help="Used for debugging purposes",
        action="store_true",
    )
    phylogenomic_alns = phylogenomic.add_argument_group(
        "Intermediate fasta files options"
    )
    phylogenomic_alns.add_argument(
        "--no_keep_fasta",
        help="Remove the orthologous groups fasta files after completion",
        action="store_true",
    )
    # TODO:  Need to add other tree calculation methods and make these options
    # mutually exclusive when the others are chosen
    phylogenomic_tree = phylogenomic.add_argument_group("IQTree2 options")
    phylogenomic_tree.add_argument(
        "--tree_model",
        help="Specific evolutionary model for tree calculation with IQTree2. Ignored by other methods",
        default=software_defaults["phylogenomic"]["tree_model"],
    )


def add_eggnog_subparser(subparsers):
    eggnog = subparsers.add_parser(
        "eggnog",
        help="The eggnog module. Use eggnogmapper to annotate the COG categories of the provided proteomes.\nIf excel files with core/fingerprint genes/proteins are provided the module will also perform an enrichment analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    eggnog_input = eggnog.add_mutually_exclusive_group()
    eggnog_input.add_argument(
        "-in",
        metavar="in",
        help="Excel file with core genes/proteins (from 'core' module)",
    )
    eggnog_input.add_argument(
        "-in_list",
        metavar="in",
        help="List of excel files with core genes/proteins (for multiple reference strains)",
    )
    # TODO: -in and -in_list should be mutually exclusive and either one should be required
    # If i do this I can get rid of the if statement in run_eggnog function
    eggnog.add_argument(
        "-fasta_dir", metavar="fasta_dir", help="Input fasta directory", required=True
    )
    eggnog.add_argument("-out", metavar="out", help="Output directory", required=True)
    eggnog.add_argument("--debug", help="Print debug information", action="store_true")
    eggnog_mapper = eggnog.add_argument_group("eggNOG mapper options")
    eggnog_mapper.add_argument(
        "--cores",
        help="Number of cores",
        default=software_defaults["eggnog"]["cores"],
        type=int,
    )
    eggnog_mapper.add_argument(
        "--pident",
        metavar="N",
        help="Percent identity",
        default=software_defaults["eggnog"]["pident"],
    )
    eggnog_mapper.add_argument(
        "--qcov",
        metavar="N",
        help="Query coverage",
        default=software_defaults["eggnog"]["qcov"],
    )
    eggnog_mapper.add_argument(
        "--scov",
        metavar="N",
        help="Suject coverage",
        default=software_defaults["eggnog"]["scov"],
    )
    eggnog_mapper.add_argument(
        "--input_type",
        help="Type of input [prot or cds]",
        default=software_defaults["eggnog"]["input_type"],
        choices=software_defaults["eggnog"]["valid_input_types"],
    )
    eggnog_mapper.add_argument(
        "--db",
        help="Non default database directory. This directory is expected to contain a 'eggNOG' subdirectory",
    )


def add_smbgc_subparser(subparsers):
    smbgc = subparsers.add_parser(
        "smbgc",
        help="smbgc module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    smbgc.add_argument(
        "-fasta_dir",
        metavar="fasta_dir",
        help="Directory with genomic fasta files",
        required=True,
    )
    smbgc.add_argument(
        "-out", metavar="output_dir", help="Output directory", required=True
    )
    smbgc.add_argument(
        "--strictness",
        metavar="strictness",
        help="antiSMASH strictness settings (loose, relaxed, strict)",
        default=software_defaults["smbgc"]["strictness"],
        choices=software_defaults["smbgc"]["valid_strictness"],
        type=str,
    )
    smbgc.add_argument(
        "--genefinding_tool",
        metavar="tool",
        help="Gene finding tool used by antiSMASH",
        default=software_defaults["smbgc"]["genefinding_tool"],
        choices=software_defaults["smbgc"]["genefinding_tools"],
        type=str,
    )
    smbgc.add_argument(
        "--db",
        help="Non default database directory. This directory is expected to contain a 'smBGC' subdirectory",
    )
    smbgc.add_argument(
        "--cores",
        help="Number of cores",
        default=software_defaults["smbgc"]["cores"],
        type=int,
    )
    smbgc.add_argument(
        "--no_concurrent",
        help="Do not run multiple searches in parallel",
        action="store_true",
    )
    smbgc.add_argument("--debug", help="Print debug information", action="store_true")


def add_download_subparser(subparsers):
    download = subparsers.add_parser(
        "download",
        help="download module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Used to download genomes from NCBI datasets",
    )
    download.add_argument("-out", metavar="out", help="Output directory", required=True)
    download.add_argument(
        "-taxon", metavar="taxon", help="Taxon to search", required=True
    )
    download.add_argument(
        "--source",
        help="Assembly source [either: RefSeq or GenBank]",
        default="RefSeq",
        choices=software_defaults["download"]["valid_sources"],
    )
    download.add_argument(
        "--level",
        help="Assembly level of genomes [available options: contig, scaffold, chromosome, complete]",
        default=software_defaults["download"]["assembly_level"],
    )
    download.add_argument(
        "--keep_plasmids",
        help="Keep plasmids for analyses. Removed by default",
        action="store_true",
    )
    download.add_argument(
        "--keep_download",
        help="Keep downloaded archive",
        action="store_true",
    )
    download.add_argument(
        "--perform_fastani",
        help="Perform fastANI once the genomes are downloaded",
        action="store_true",
    )
    download.add_argument(
        "--debug", help="Print debug information", action="store_true"
    )


def add_virulence_subparser(subparsers):
    virulence = subparsers.add_parser(
        "virulence",
        help="Virulence identification module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Used to identify virulence factors with VFDB",
    )
    virulence.add_argument(
        "-fasta_dir", metavar="in", help="Fasta directory", required=True
    )
    virulence.add_argument(
        "-out", metavar="out", help="Output directory", required=True
    )
    virulence.add_argument(
        "--input_type",
        help="Input type [either: prot or cds]",
        default=software_defaults["virulence"]["input_type"],
        choices=software_defaults["virulence"]["valid_input_types"],
    )
    virulence.add_argument(
        "--db",
        help="Non default database directory. This directory is expected to contain a 'Virulence' subdirectory",
    )
    virulence.add_argument(
        "--cores",
        help="Number of cores",
        default=software_defaults["virulence"]["blast_cores"],
        type=int,
    )
    virulence.add_argument(
        "--evalue",
        help="E-value cut-off",
        default=software_defaults["virulence"]["blast_evalue"],
    )
    virulence.add_argument(
        "--dmnd_sensitivity",
        help="Sensitivity settings for DIAMOND (i.e. very_sensitive)",  # TODO: Is this example good? Should it be _ or -?
        default=software_defaults["virulence"]["dmnd_sensitivity"],
    )
    virulence.add_argument(
        "--no_concurrent",
        help="Do not run parallel search jobs",
        action="store_true",
    )

    virulence.add_argument(
        "--debug", help="Print debug information", action="store_true"
    )


def add_cazy_subparser(subparsers):
    cazy = subparsers.add_parser(
        "cazy",
        help="CAZYme identification module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Used to identify CAZYmes with dbcan",
    )
    cazy.add_argument("-fasta_dir", metavar="in", help="Fasta directory", required=True)
    cazy.add_argument("-out", metavar="out", help="Output directory", required=True)
    cazy.add_argument(
        "--input_type",
        help="Input type [either: prot, nucl or cds]",
        default=software_defaults["cazy"]["input_type"],
        choices=software_defaults["cazy"]["valid_input_types"],
    )
    cazy.add_argument(
        "--db",
        help="Non default database directory. This directory is expected to contain a 'CAZY' subdirectory",
    )
    cazy.add_argument(
        "--cores",
        help="Number of cores",
        default=software_defaults["cazy"]["cores"],
        type=int,
    )
    cazy.add_argument(
        "--evalue",
        help="E-value cut-off",
        default=software_defaults["cazy"]["evalue"],
    )

    cazy.add_argument("--debug", help="Print debug information", action="store_true")


def add_amr_subparser(subparsers):
    amr = subparsers.add_parser(
        "amr",
        help="Antimicrobial resistance (AMR) genes identification module. This module will run multiple search jobs by default.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Used to identify AMR genes with amrfinder",
    )
    amr.add_argument("-fasta_dir", metavar="in", help="Fasta directory")
    amr.add_argument("-out", metavar="out", help="Output directory")
    amr.add_argument(
        "--input_type",
        help="Input type [either: prot or cds]",
        default=software_defaults["amr"]["input_type"],
        choices=software_defaults["amr"]["valid_input_types"],
    )
    amr.add_argument(
        "--db",
        help="Non default database directory. This directory is expected to contain a 'AMR' subdirectory",
    )
    amr.add_argument(
        "--cores",
        help="Number of cores",
        default=software_defaults["amr"]["cores"],
        type=int,
    )
    amr.add_argument(
        "--no_concurrent",
        help="Do not run parallel search jobs",
        action="store_true",
    )
    amr.add_argument("--debug", help="Print debug information", action="store_true")


def add_workflow_subparser(subparsers):
    workflow = subparsers.add_parser(
        "workflow",
        help="workflow module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Wrapper module for automated analyses",
    )
    workflow.add_argument("-in", metavar="in", help="Parameters file", required=True)
    workflow.add_argument("-out", metavar="out", help="Output directory", required=True)
    workflow.add_argument("-db", metavar="db", help="Database directory")
    workflow.add_argument(
        "--debug", help="Print debug information", action="store_true"
    )


def add_databases_subparser(subparsers):
    databases = subparsers.add_parser(
        "databases",
        help="databases module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Download the required databases (smbgc, eggnog, cazy, virulence, amr)",
    )
    databases.add_argument(
        "--install",
        metavar="db",
        help="Which database to install: Either 'all' or comma separated values i.e. 'smbgc,eggnog,cazy'",
        required=True,
        default=software_defaults["database"]["db"],
    )
    databases.add_argument(
        "--database_dir",
        metavar="db_dir",
        help="Specific download directory",
        default=software_defaults["database"]["db_dir"],
    )
    databases.add_argument(
        "--cores",
        help="Number of cores for parallel downloads. Applicable to certain databases. Use a low value (i.e. 4) to avoid throttling",
        default=software_defaults["database"]["cores"],
        type=int,
    )
    databases.add_argument(
        "--debug", help="Print debug information", action="store_true"
    )


def run_species_demarcation(args):
    in_dir = Path(args["in"])
    if not validate_directory(in_dir):
        logging.error(in_dir)
        return
    out_dir = Path(args["o"])
    fastani_cores = args["fastani_cores"]
    debug = args["debug"]
    kmer = int(args["kmer"])
    fraglen = int(args["fraglen"])
    minfraction = float(args["minfraction"])
    inflation = float(args["inflation"])
    mcl_cores = int(args["mcl_cores"])
    demarcator = SpeciesDemarcator(
        in_dir=in_dir,
        out_dir=out_dir,
        fastani_cores=fastani_cores,
        kmer=kmer,
        fraglen=fraglen,
        minfrac=minfraction,
        inflation=inflation,
        mcl_cores=mcl_cores,
        debug=debug,
    )
    demarcator.assign_species()


def run_orthologues(args):
    # fasta_in_dir = Path(args["in"])
    fasta_files = args["fasta_files"]
    out_dir = Path(args["out"])
    ref = args["ref"]
    ref_list = args["ref_list"]
    input_type = args["input_type"]
    cores = int(args["cores"])
    evalue = float(args["evalue"])
    dmnd_sensitivity = args["dmnd_sensitivity"]
    no_filter_orthologues = args["no_filter_orthologues"]
    concurrent = True
    if args["no_concurrent"]:
        concurrent = False
    orthologues_identifier = Orthologues_identifier(
        fasta_files_list=fasta_files,
        out_dir=out_dir,
        ref=ref,
        ref_list=ref_list,
        input_type=input_type,
        available_cores=software_defaults["system"]["usable_cores"],
        blast_cores=cores,
        concurrent=concurrent,
        evalue=evalue,
        dmnd_sensitivity=dmnd_sensitivity,
        no_filter=no_filter_orthologues,
    )
    orthologues_identifier.calculate_orthologues()


def run_core(args):
    out_dir = Path(args["out"])
    species_file = args.get("species")
    if species_file:
        species_file = Path(species_file)
        species_df = read_excel(species_file, index_col=0)
    else:
        species_df = None
    og_matrix_in = args.get("in")
    og_matrix_list_f = args.get("ref_list")
    if og_matrix_list_f is None:
        og_matrix_list = [og_matrix_in]
    else:
        og_matrix_list = []
        with open(og_matrix_list_f, "r") as rf:
            for line in rf:
                line = line.rstrip()
                og_matrix_list.append(line)
    og_matrix_list = [Path(og_matrix_in) for og_matrix_in in og_matrix_list]
    core_perc = float(args["core_perc"])
    for og_matrix_in in tqdm(og_matrix_list, ascii=True, desc="Calculating core"):
        core_identifier = Core_identifier(
            orthology_fin=og_matrix_in,
            out_dir=out_dir,
            species_df=species_df,
            core_perc=core_perc,
        )
        core_identifier.calculate_core()


def run_phylogenomic(args):
    og_matrix_in = Path(args["in"])
    out_dir = Path(args["out"])
    cores = int(args["cores"])
    no_keep_fasta = args["no_keep_fasta"]
    tree_model = args["tree_model"]
    tree_method = args["method"]
    input_type = args["input_type"]
    debug = args["debug"]
    fasta_files = args["fasta_files"]
    phylogenomic = Phylogenomic(
        orthology_matrix_in=og_matrix_in,
        cores=cores,
        fasta_files_list=fasta_files,
        out_dir=out_dir,
        no_keep_fasta=no_keep_fasta,
        tree_method=tree_method,
        iqtree_model=tree_model,
        input_type=input_type,
        debug=debug,
    )
    phylogenomic.run_phylogenomic()


def run_eggnog(args):
    debug = args["debug"]
    out_dir = Path(args["out"])
    cores = int(args["cores"])
    pident = args["pident"]
    qcov = args.get("qcov")
    scov = args.get("scov")
    input_type = args.get("input_type")

    if input_type is None:
        raise RuntimeError("eggNOG: Input type cannot be None")

    if args.get("db") is None:
        database_dir = Path(DB_BASE_DIR) / "site-packages" / "data"
    else:
        database_dir = Path(args["db"]) / "eggNOG"

    fasta_files = args["fasta_files"]

    run_inference = args.get("execute_eggnog_mapper", True)
    # When calling the module from the cli run_inference will be None
    # Thus default it to True
    if run_inference:
        # First run the references
        runner = EggNOGRunner(
            fasta_files_list=fasta_files,
            out_dir=out_dir,
            cores=cores,
            pident=pident,
            qcov=qcov,
            scov=scov,
            input_type=input_type,
            database_dir=database_dir,
            debug=debug,
        )
        runner.execute_eggnog_mapper()

    core_proteins_file = args.get("in")
    core_protein_files_reflist = args.get("in_list")
    if core_proteins_file is not None or core_protein_files_reflist is not None:
        core_proteins_file_list = (
            [Path(core_proteins_file)]
            if core_proteins_file
            else [Path(line.rstrip()) for line in open(core_protein_files_reflist, "r")]
        )
    else:
        core_proteins_file_list = None

    parser = EggNOGParser(
        fasta_files_list=fasta_files,
        core_proteins_table_f=core_proteins_file_list,  # Accepts List[Path] or None
        out_dir=out_dir,
        debug=debug,
    )
    parser.gather_eggnog_proteome_results()
    if core_proteins_file_list is not None:
        parser.compare_proteome_to_core_and_fps()


def run_smbgc(args):
    logging.info(
        f"Starting antiSMASH run: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
    )
    debug = args["debug"]
    out_dir = Path(args.get("out"))
    strictness = args.get("strictness")
    genefinding_tool = args.get("genefinding_tool")
    cores = int(args.get("cores", 0))

    database_dir = args.get("db", None)
    if database_dir is not None:
        database_dir = Path(database_dir) / "smBGC"
    # antismash gets the default database if nothing is supplied

    concurrent = True
    if args["no_concurrent"]:
        concurrent = False

    fasta_files = args["fasta_files"]

    # Initialize runner and analyze the genomes
    local_runner = smBGCLocalRunner(
        genome_fasta_files=fasta_files,
        out_dir=out_dir,
        strictness=strictness,
        genefinding_tool=genefinding_tool,
        database_dir=database_dir,
        cores=cores,
        available_cores=software_defaults["system"]["usable_cores"],
        concurrent=concurrent,
        debug=debug,
    )
    local_runner.analyze_genomes()
    parser = smBGCParser(out_dir, cores)
    parser.gather_results()
    logging.info(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")


def run_virulence(args):
    out_dir = Path(args["out"])
    input_type = args["input_type"]
    cores = int(args["cores"])
    evalue = float(args["evalue"])
    dmnd_sensitivity = args["dmnd_sensitivity"]
    debug = args["debug"]
    if args.get("db") is None:
        database_dir = Path(DB_BASE_DIR)
    else:
        database_dir = Path(args["db"])
    database_dir = database_dir / "Virulence"
    concurrent = True
    if args["no_concurrent"]:
        concurrent = False
    fasta_files = args["fasta_files"]
    vf_analyzer = VF_analyzer(
        fasta_files_list=fasta_files,
        database_dir=database_dir,
        out_dir=out_dir,
        input_type=input_type,
        available_cores=software_defaults["system"]["usable_cores"],
        blast_cores=cores,
        evalue=evalue,
        dmnd_sensitivity=dmnd_sensitivity,
        concurrent=concurrent,
        debug=debug,
    )
    vf_analyzer.execute_homology_search()
    vf_analyzer.parse_results()


def run_cazy(args):
    out_dir = Path(args["out"])
    input_type = args["input_type"]
    cores = int(args["cores"])
    evalue = float(args["evalue"])
    if args.get("db") is None:
        database_dir = Path(DB_BASE_DIR)
    else:
        database_dir = Path(args["db"])
    database_dir = database_dir / "CAZY"
    concurrent = True
    if args["no_concurrent"]:
        concurrent = False
    debug = args["debug"]
    fasta_files = args["fasta_files"]

    cazy_analyzer = CAZY_analyzer(
        fasta_files_list=fasta_files,
        out_dir=out_dir,
        database_dir=database_dir,
        input_type=input_type,
        available_cores=software_defaults["system"]["usable_cores"],
        cores=cores,
        concurrent=concurrent,
        evalue=evalue,
        debug=debug,
    )
    cazy_analyzer.execute_cazy_search()


def run_amr(args):
    out_dir = Path(args["out"])
    input_type = args["input_type"]
    cores = int(args["cores"])
    if args.get("db") is None:
        database_dir = Path(DB_BASE_DIR) / "AMR"
    else:
        database_dir = Path(args["db"]) / "AMR"
    concurrent = True
    if args["no_concurrent"]:
        concurrent = False
    debug = args["debug"]
    fasta_files = args["fasta_files"]
    amr_analyzer = AMR_analyzer(
        fasta_files_list=fasta_files,
        database_dir=database_dir,
        out_dir=out_dir,
        input_type=input_type,
        available_cores=software_defaults["system"]["usable_cores"],
        cores=cores,
        concurrent=concurrent,
        debug=debug,
    )
    amr_analyzer.search_amr()
    amr_analyzer.parse_results()


def run_download(args):
    out_dir = Path(args["out"])
    if not validate_directory(out_dir):
        return
    taxon = args["taxon"]
    debug = args["debug"]
    downloader = GenomeDownloader(
        taxon=taxon,
        out_dir=out_dir,
        assembly_level=args["level"],
        assembly_source=args["source"],
        keep_plasmids=args["keep_plasmids"],
        debug=debug,
    )
    logging.info(f"Starting download: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
    downloader.download_hydrated()
    downloader.extract_dataset_zip()
    downloader.create_output_directories()
    downloader.process_gbff_files()
    downloader.write_annotations()
    downloader.write_16S_fasta()
    if not args["keep_download"]:
        downloader.remove_dataset_archive()

    if args["perform_fastani"]:
        logging.info(
            f"Finished downloading, initiating FastANI: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        in_dir = out_dir / "Genomic_fasta_files"
        if not validate_directory(in_dir):
            logging.error(f"{in_dir} is invalid")
            return

        demarcator = SpeciesDemarcator(
            in_dir=in_dir,
            out_dir=out_dir,
            fastani_cores=software_defaults["species_demarcation"]["cores"],
            kmer=software_defaults["species_demarcation"]["kmer"],
            fraglen=software_defaults["species_demarcation"]["fraglen"],
            minfrac=software_defaults["species_demarcation"]["minfrac"],
            inflation=software_defaults["species_demarcation"]["inflation"],
            mcl_cores=software_defaults["species_demarcation"]["cores"],
            debug=debug,
        )
        demarcator.assign_species()
        # Join excel files
        f1 = out_dir / "Genome_information.xlsx"
        f2 = out_dir / "Species_demarcation" / "FastANI_species_clusters.xlsx"
        df1 = read_excel(f1, index_col=0)
        df2 = read_excel(f2, index_col=0)
        df1 = utils.join_pandas_dataframes(df1, df2)
        df1.to_excel(f1)
    logging.info(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")


def run_databases(args):
    db_value = args["db"]
    database_dir = args["db_dir"]
    debug = args["debug"]
    cores = args["cores"]
    if "all" in db_value and db_value != "all":
        raise ValueError(
            "Cannot use the value 'all' with other databases during installation"
        )
    if db_value == "all":
        databases_to_install = ["smbgc", "eggnog", "cazy", "virulence", "amr"]
    else:
        databases_to_install = db_value.split(",")

    print(databases_to_install)
    for database in databases_to_install:
        if database == "smbgc":
            smbgc = SMBGC_installer(database_dir=database_dir, debug=debug)
            smbgc.install_database()
        if database == "eggnog":
            eggnog = EGGNOG_installer(database_dir=database_dir, debug=debug)
            eggnog.install_database()
        if database == "cazy":
            cazy = CAZY_installer(database_dir=database_dir, cores=cores, debug=debug)
            cazy.install_database()
        if database == "amr":
            amr = AMR_installer(database_dir=database_dir, debug=debug)
            amr.install_database()
        if database == "virulence":
            virulence = VF_installer(database_dir=database_dir, debug=debug)
            virulence.install_database()
        return None


def run_workflow(args):
    arguments = software_defaults.copy()
    fin = Path(args["in"])
    out = Path(args["out"])
    debug = args.get("debug", False)
    db_dir = args.get("db", None)

    # Update the default arguments with the CLI parameters
    for key in arguments:
        arguments[key]["out"] = out
        arguments[key]["db"] = db_dir
        arguments[key]["debug"] = debug

    # Initialize the workflow runner and read the parameters excel file
    workflow_runner = WorkflowRunner(out_dir=out, param_file=fin, debug=debug)
    workflow_runner.read_param_file()
    workflow_runner.validate_parameters()

    # Identify which tasks the user wants to run
    workflow_runner.identify_tasks()

    # Update the parameters of each task using the parameters file
    input_type = workflow_runner.params.get("CDS_or_proteins")

    # Gather the fasta files that will be used
    workflow_runner.gather_cds_or_protein_fasta_files()
    # workflow_runner.cds_or_protein_fasta_files attribute contains the whole set
    workflow_runner.gather_genomic_fasta_files()
    # workflow_runner.genomic_fasta_files attribute contains the whole set
    workflow_runner.get_per_group_representatives_files()
    # workflow_runner.cds_or_protein_fasta_files_representatives
    # workflow_runner.genomic_fasta_files_representatives

    for task in workflow_runner.tasks:
        if task == "Calculate_orthologues":
            module = "orthologues"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files
            )
            arguments[module]["ref"] = None
            arguments[module]["ref_list"] = workflow_runner.orthologues_ref_list
            arguments[module]["input_type"] = input_type
            run_orthologues(arguments[module])

        if task == "Calculate_cores" or task == "Calculate_fingerprints":
            module = "core"
            # Calculate core of whole set and of individual groups

            # The orthologous matrix list
            workflow_runner.create_core_ref_list()

            # Create a temporary species file
            workflow_runner.write_species_file_for_core()

            for runtype in ["whole", "groups"]:
                if runtype == "whole":
                    arguments[module]["ref_list"] = (
                        workflow_runner.core_whole_set_ref_list
                    )
                    arguments[module]["species"] = None
                else:
                    arguments[module]["ref_list"] = workflow_runner.core_ref_list
                    arguments[module]["species"] = workflow_runner.core_species_file
                run_core(arguments[module])

        if task == "Calculate_entire_phylogenomic_tree":
            module = "phylogenomic"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files
            )
            arguments[module]["input_type"] = input_type
            arguments[module]["method"] = workflow_runner.params.get(
                "Calculate_entire_phylogenomic_tree", "IQTree"
            )
            workflow_runner.get_phylogenomic_orthologue_matrix_file()
            arguments[module]["in"] = workflow_runner.phylogenomic_og_matrix
            arguments[module]["out"] = arguments[module]["out"] / "WholeSet"
            run_phylogenomic(arguments[module])

        if task == "Calculate_group_representatives_phylogenomic_tree":
            module = "phylogenomic"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files_representatives
            )
            arguments[module]["input_type"] = input_type
            arguments[module]["method"] = workflow_runner.params.get(
                "Calculate_group_representatives_phylogenomic_tree", "IQTree"
            )
            workflow_runner.get_phylogenomic_orthologue_matrix_file()
            arguments[module]["in"] = workflow_runner.phylogenomic_og_matrix
            arguments[module]["out"] = arguments[module]["out"] / "Representatives"
            arguments[module]["execute_eggnog_mapper"] = True
            run_phylogenomic(arguments[module])

        if task == "Calculate_EGGNOG_on_Group_representatives":
            module = "eggnog"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files_representatives
            )
            arguments[module]["input_type"] = input_type
            if arguments[module].get("execute_eggnog_mapper") is True:
                arguments[module]["execute_eggnog_mapper"] = False
                # In order to skip running eggnog mapper if
                # it was already run from the previous task
                # This option will not run during the cli call of the module
                # in run_eggnog function the return of arg[module].get("execute..") is default True
            run_eggnog(arguments[module])

        if task == "Calculate_EGGNOG_on_core/fingerprints":
            module = "eggnog"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files_representatives
            )
            arguments[module]["input_type"] = input_type

            # Get the core protein and fingerprint excel files
            workflow_runner.create_eggnog_core_protein_files_list()
            arguments[module]["in_list"] = (
                workflow_runner.emapper_core_protein_files_reflist
            )
            run_eggnog(arguments[module])

        if task == "Calculate_SMBGCs_on_Group_representatives":
            module = "smbgc"
            arguments[module]["fasta_files"] = (
                workflow_runner.genomic_fasta_files_representatives
            )
            run_smbgc(arguments[module])

        if task == "Calculate_SMBGCs_on_entire_set":
            module = "smbgc"
            arguments[module]["fasta_files"] = workflow_runner.genomic_fasta_files
            run_smbgc(arguments[module])

        if task == "Calculate_CAZymes_on_Group_representatives":
            module = "cazy"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files_representatives
            )
            arguments[module]["input_type"] = input_type
            run_cazy(arguments[module])

        if task == "Calculate_CAZYmes_on_entire_set":
            module = "cazy"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files
            )
            arguments[module]["input_type"] = input_type
            run_cazy(arguments[module])

        if task == "Calculate_VFs_on_Group_representatives":
            module = "virulence"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files_representatives
            )
            arguments[module]["input_type"] = input_type
            run_virulence(arguments[module])

        if task == "Calculate_VFs_on_entire_set":
            module = "virulence"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files
            )
            arguments[module]["input_type"] = input_type
            run_virulence(arguments[module])

        if task == "Calculate_AMR_on_Group_representatives":
            module = "amr"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files_representatives
            )
            arguments[module]["input_type"] = input_type
            run_amr(arguments[module])

        if task == "Calculate_AMR_on_entire_set":
            module = "amr"
            arguments[module]["fasta_files"] = (
                workflow_runner.cds_or_protein_fasta_files
            )
            arguments[module]["input_type"] = input_type
            run_amr(arguments[module])

    return None


def main():
    parser = setup_parser()
    args = vars(parser.parse_args())

    if args["module"] == "species_demarcation":
        run_species_demarcation(args)

    elif args["module"] == "orthologues":
        fasta_in_dir = Path(args["in"])
        if not validate_directory(fasta_in_dir):
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_orthologues(args)

    elif args["module"] == "core":
        run_core(args)

    elif args["module"] == "phylogenomic":
        fasta_in_dir = Path(args["fasta_dir"])
        if not validate_directory(fasta_in_dir):
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_phylogenomic(args)

    elif args["module"] == "eggnog":
        fasta_in_dir = Path(args["fasta_dir"])
        if not validate_directory(fasta_in_dir):
            logging.error("eggNOG: fasta directory is invalid")
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_eggnog(args)

    elif args["module"] == "smbgc":
        fasta_in_dir = Path(args["fasta_dir"])
        if not validate_directory(fasta_in_dir):
            logging.error("Antismash: fasta directory is invalid")
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_smbgc(args)

    elif args["module"] == "virulence":
        fasta_in_dir = Path(args["fasta_dir"])
        if not validate_directory(fasta_in_dir):
            logging.error("Virulence: fasta directory is invalid")
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_virulence(args)

    elif args["module"] == "cazy":
        fasta_in_dir = Path(args["fasta_dir"])
        if not validate_directory(fasta_in_dir):
            logging.error("CAZY: fasta directory is invalid")
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_cazy(args)

    elif args["module"] == "amr":
        fasta_in_dir = Path(args["fasta_dir"])
        if not validate_directory(fasta_in_dir):
            logging.error("AMR: fasta directory is invalid")
            return
        fasta_files = list(fasta_in_dir.glob("*"))
        args["fasta_files"] = fasta_files
        run_amr(args)

    elif args["module"] == "download":
        run_download(args)

    elif args["module"] == "databases":
        run_databases(args)

    elif args["module"] == "workflow":
        run_workflow(args)
