"""
Author: Marios Nikolaidis
Git: https://github.com/m-nikolaidis
email: marionik23@gmail.com
"""
import argparse
from pathlib import Path
from datetime import datetime

from tqdm import tqdm

from pypgcf import config
from pypgcf import checks
from pypgcf.core import Core_identifier
from pypgcf.eggnog import eggNOGInstaller, eggNOGParser, eggNOGRunner
from pypgcf.orthologues import Orthologues_identifier
from pypgcf.phylogenomic import Phylogenomic
from pypgcf.smbgc import smBGCInstaller, smBGCLocalRunner, smBGCParser
from pypgcf.species_demarcation import SpeciesDemarcator
from pypgcf.download_genomes import GenomeDownloader


def main():
    parser = argparse.ArgumentParser(
        prog="pyPGCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="pyPGCF: PhyloGenomic, Core and Fingerprint analysis software",
        allow_abbrev=False,
    )
    subparsers = parser.add_subparsers(
        help="The various modules of the program", dest="module", required=True
    )

    # 1. species_demarcation module
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
        "--debug",
        help="Used for debugging purposes",
        action="store_true",
    )

    species_demarcation_fastani = species_demarcation.add_argument_group(
        "FastANI options"
    )
    species_demarcation_fastani.add_argument(
        "--fastani_cores",
        help="Number of cores for FastANI",
        default=config.species_demarcation_cores,
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--kmer",
        metavar="N",
        help="kmer size (<= 16)",
        default=config.species_demarcation_kmer,
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--fraglen",
        metavar="N",
        help="Fragment length",
        default=config.species_demarcation_fraglen,
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--minfraction",
        metavar="N",
        help="Minimum fraction",
        default=config.species_demarcation_minfraction,
        type=float,
    )

    species_demarcation_mcl = species_demarcation.add_argument_group("MCL options")
    species_demarcation_mcl.add_argument(
        "--inflation",
        help="Inflation parameter for MCL",
        default=config.species_demarcation_mcl_inflation,
    )
    species_demarcation_mcl.add_argument(
        "--mcl_cores",
        help="Number of cores for MCL",
        default=config.species_demarcation_cores,
        type=int,
    )

    # orthologues module
    orthologues = subparsers.add_parser(
        "orthologues",
        help="orthologues module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Calculate the orthologues of the dataset using one reference strain",
    )
    orthologues_basic = orthologues.add_argument_group("Basic options")
    orthologues_basic.add_argument(
        "-in", metavar="in", help="Fasta directory", required=True
    )
    orthologues_basic.add_argument(
        "-o", metavar="out", help="Output directory", required=True
    )
    ref_exclusive = orthologues_basic.add_mutually_exclusive_group(required=True)
    ref_exclusive.add_argument("-ref", help="Reference strain")
    ref_exclusive.add_argument("-ref_list", help="List of reference strains to use")
    orthologues_basic.add_argument(
        "--type", help="Type of input [prot (DIAMOND) or nucl (BLASTN)]", default="prot"
    )
    orthologues_blast = orthologues.add_argument_group("DIAMOND/BLASTN options")
    orthologues_blast.add_argument(
        "--cores", help="Number of cores", default=config.orthologues_cores
    )
    orthologues_blast.add_argument(
        "--evalue", help="E-value cut-off", default=config.orthologues_evalue
    )
    orthologues_blast.add_argument(
        "--dmnd_sensitivity",
        help="Sensitivity settings for DIAMOND i.e. very_sensitive",
        default=config.orthologues_dmnd_sensitivity,
    )
    orthologues_blast.add_argument(
        "--no_filter_orthologues", help="Do not filter orthologues", action="store_true"
    )

    # core module
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
    core.add_argument("-o", help="Output directory", required=True)
    core.add_argument(
        "--species",
        help="Input species assignment matrix (excel format; from 'species_demarcation' module)",
    )
    core.add_argument(
        "--core_perc",
        help="Percent presence of a protein/gene in a cluster to be considered core",
        default=config.core_core_perc,
    )

    # phylogenomic module
    phylogenomic = subparsers.add_parser(
        "phylogenomic",
        help="phylogenomic module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    phylogenomic.add_argument("-fasta_dir", help="Input fasta directory", required=True)
    phylogenomic.add_argument(
        "-in", help="Input orthology matrix (from 'orthologues' module)", required=True
    )
    phylogenomic.add_argument("-o", help="Output directory", required=True)
    phylogenomic.add_argument(
        "--cores", help="Number of cores", default=config.phylogenomic_cores
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
    phylogenomic_tree = phylogenomic.add_argument_group("IQTree2 options")
    phylogenomic_tree.add_argument(
        "--tree_model", help="Specific evolutionary model for tree calculation"
    )

    # eggnog module
    eggnog = subparsers.add_parser(
        "eggnog",
        help="eggnog module",
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
    eggnog.add_argument("-fasta_dir", metavar="fasta_dir", help="Input fasta directory")
    eggnog.add_argument("-o", metavar="o", help="Output directory")
    eggnog.add_argument("--debug", help="Print debug information", action="store_true")
    eggnog_mapper = eggnog.add_argument_group("eggNOG mapper options")
    eggnog_mapper.add_argument(
        "--cores",
        metavar="N",
        help="Number of cores",
        default=config.eggnog_cores,
        type=int,
    )
    eggnog_mapper.add_argument(
        "--pident", metavar="N", help="Percent identity", default=config.eggnog_pident
    )
    eggnog_mapper.add_argument(
        "--qcov", metavar="N", help="Query coverage", default=config.eggnog_qcov
    )
    eggnog_mapper.add_argument(
        "--scov", metavar="N", help="Suject coverage", default=config.eggnog_scov
    )
    eggnog_mapper.add_argument(
        "--nucl", help="Specify that the input type is CDS", action="store_true"
    )
    eggnog_installer = eggnog.add_argument_group("Installer options")
    eggnog_installer.add_argument(
        "--install", help="Install the eggNOG database", action="store_true"
    )

    # smbgc module
    smbgc = subparsers.add_parser(
        "smbgc",
        help="smbgc module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    smbgc.add_argument(
        "-fasta_dir",
        metavar="fasta_dir",
        help="Directory with genomic fasta files",
    )
    smbgc.add_argument("-o", metavar="output_dir", help="Output directory")
    smbgc.add_argument(
        "--cores",
        metavar="N",
        help="Number of cores",
        default=config.smbgc_cores,
        type=int,
    )
    smbgc.add_argument(
        "--strictness",
        metavar="strictness",
        help="antiSMASH strictness settings (loose, relaxed, strict)",
        default=config.smbgc_strictness,
        type=str,
    )
    smbgc.add_argument(
        "--genefinding_tool",
        metavar="tool",
        help="Gene finding tool used by antiSMASH",
        default=config.smbgc_genefinding_tool,
        type=str,
    )
    smbgc.add_argument(
        "--install",
        help="Install the required antiSMASH databases",
        action="store_true",
    )
    smbgc.add_argument("--debug", help="Print debug information", action="store_true")
    # smbgc.add_argument("--remote", help="Submit queries to antiSMASH web service", action="store_true")

    # download module
    download = subparsers.add_parser(
        "download",
        help="download module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Used to download genomes from NCBI datasets",
    )
    download.add_argument(
        "-o", metavar="out", help="Output directory", required=True
    )
    download.add_argument(
        "-t", metavar="taxon", help="Taxon to search", required=True
    )
    # ref_exclusive = orthologues_basic.add_mutually_exclusive_group(required=True) use this to provide specific gffs from file
    # ref_exclusive.add_argument("-ref", help="Reference strain")
    # ref_exclusive.add_argument("-ref_list", help="List of reference strains to use")
    download.add_argument(
        "--source", help="Assembly source [either: RefSeq or GenBank]", default="RefSeq"
    )
    download.add_argument(
        "--level", help="Assembly level of genomes [available options: contig, scaffold, chromosome, complete]", default="chromosome,complete"
    )
    download.add_argument("--keep_plasmids", help="Keep plasmids for analyses", action="store_true")
    download.add_argument("--keep_download", help="Keep downloaded archive", action="store_true")
    download.add_argument("--debug", help="Print debug information", action="store_true")


    # Parse arguments
    args = vars(parser.parse_args())
    if args["module"] == "species_demarcation":
        in_dir = Path(args["in"])
        if not checks.check_if_dir_exists(in_dir):
            print(f"{in_dir} does not exist")
            return
        if checks.check_if_dir_is_empty(in_dir):
            print(f"{in_dir} is empty")
            return
        out_dir = Path(args["o"])
        fastani_cores = args["fastani_cores"]
        debug = args["debug"]
        kmer = int(args["kmer"])
        fraglen = int(args["fraglen"])
        minfraction = float(args["minfraction"])
        inflation = float(args["inflation"])
        mcl_cores = args["mcl_cores"]
        demarcator = SpeciesDemarcator(
            in_dir,
            out_dir,
            fastani_cores,
            kmer,
            fraglen,
            minfraction,
            inflation,
            mcl_cores,
            debug,
        )
        demarcator.assign_species()

    if args["module"] == "orthologues":
        fasta_in_dir = Path(args["in"])
        if not checks.check_if_dir_exists(fasta_in_dir):
            print(f"{fasta_in_dir} does not exist")
            return
        if checks.check_if_dir_is_empty(fasta_in_dir):
            print(f"{fasta_in_dir} is empty")
            return
        out_dir = Path(args["o"])
        ref = args["ref"]
        ref_list = args["ref_list"]
        input_type = args["type"]
        cores = int(args["cores"])
        evalue = args["evalue"]
        dmdnd_sensitivity = args["dmnd_sensitivity"]
        no_filter_orthologues = args["no_filter_orthologues"]
        orthologues_identifier = Orthologues_identifier(
            fasta_in_dir,
            out_dir,
            ref,
            ref_list,
            input_type,
            cores,
            evalue,
            dmdnd_sensitivity,
            no_filter_orthologues,
        )
        orthologues_identifier.calculate_orthologues()

    if args["module"] == "core":
        out_dir = Path(args["o"])
        species_file = args["species"]
        if species_file != None:
            species_file = Path(species_file)
        og_matrix_in = args["in"]
        og_matrix_list_f = args["ref_list"]
        if og_matrix_list_f == None:
            og_matrix_list = [og_matrix_in]
        else:
            og_matrix_list = []
            with open(og_matrix_list_f, "r") as rf:
                for line in rf:
                    line = line.rstrip()
                    og_matrix_list.append(line)
        og_matrix_list = [Path(og_matrix_in) for og_matrix_in in og_matrix_list]
        for og_matrix_in in tqdm(og_matrix_list, ascii=True, desc="Calculating core"):
            core_perc = float(args["core_perc"])
            core_identifier = Core_identifier(
                og_matrix_in, out_dir, species_file, core_perc
            )
            core_identifier.calculate_core()

    if args["module"] == "phylogenomic":
        fasta_dir = Path(args["fasta_dir"])
        if not checks.check_if_dir_exists(fasta_dir):
            print(f"{fasta_dir} does not exist")
            return
        if checks.check_if_dir_is_empty(fasta_dir):
            print(f"{fasta_dir} is empty")
            return
        og_matrix_in = Path(args["in"])
        out_dir = Path(args["o"])
        cores = int(args["cores"])
        no_keep_fasta = args["no_keep_fasta"]
        tree_model = args["tree_model"]
        debug = args["debug"]
        phylogenomic = Phylogenomic(
            og_matrix_in, cores, fasta_dir, out_dir, no_keep_fasta, tree_model, debug
        )
        phylogenomic.run_phylogenomic()

    if args["module"] == "eggnog":
        install = args["install"]
        debug = args["debug"]
        if install:
            installer = eggNOGInstaller(debug)
            installer.download_databases()
            return
        fasta_dir = args["fasta_dir"]
        out_dir = args["o"]
        if fasta_dir == None:
            print(f"-fasta_dir is needed")
            return
        if out_dir == None:
            print(f"-o is needed")
            return
        fasta_dir = Path(fasta_dir)
        if not checks.check_if_dir_exists(fasta_dir):
            print(f"{fasta_dir} does not exist")
            return
        if checks.check_if_dir_is_empty(fasta_dir):
            print(f"{fasta_dir} is empty")
            return
        out_dir = Path(out_dir)
        cores = int(args["cores"])
        pident = args["pident"]
        qcov = args["qcov"]  # Query coverage
        scov = args["scov"]  # Subject coverage
        core_proteins_file = args["in"]
        core_protein_files_reflist = args["in_list"]
        nucl = args["nucl"]
        if core_proteins_file == None and core_protein_files_reflist == None:
            print("-in or -in_list are needed")
            return
        if core_proteins_file != None:
            core_proteins_file_list = [core_proteins_file]
        else:
            core_proteins_file_list = []
            with open(core_protein_files_reflist, "r") as rf:
                for line in rf:
                    core_proteins_file_list.append(line.rstrip())
        core_proteins_file_list = [
            Path(core_proteins_file) for core_proteins_file in core_proteins_file_list
        ]
        for core_proteins_file in core_proteins_file_list:
            core_proteins_file = Path(core_proteins_file)
            runner = eggNOGRunner(
                fasta_dir,
                core_proteins_file,
                out_dir,
                cores,
                pident,
                qcov,
                scov,
                nucl,
                debug=debug,
            )
            runner.execute_eggnog_mapper()
        parser = eggNOGParser(fasta_dir, core_proteins_file_list, out_dir, debug)
        parser.gather_eggnog_results()

    if args["module"] == "smbgc":
        install = args["install"]
        debug = args["debug"]
        if install:
            installer = smBGCInstaller(debug)
            installer.install_databases()
            return
        print(
            f"Starting antiSMASH run: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        fasta_dir = args.get("fasta_dir")
        out_dir = args.get("o")
        if fasta_dir is None or out_dir is None:
            print("Both -fasta_dir and -o parameters are required")
            return
        fasta_dir = Path(args.get("fasta_dir", "."))
        out_dir = Path(args.get("o", "."))
        cores = int(args["cores"])
        strictness = args["strictness"]
        genefinding_tool = args["genefinding_tool"]
        # Perform all the checks before running
        if not checks.check_if_dir_exists(fasta_dir):
            print(f"{fasta_dir} does not exist")
            return
        if checks.check_if_dir_is_empty(fasta_dir):
            print(f"{fasta_dir} is empty")
            return
        if not checks.is_valid_antismash_strict(strictness):
            return
        if not checks.is_valid_genefinding_tool(genefinding_tool):
            return
        local_runner = smBGCLocalRunner(
            fasta_dir, out_dir, cores, strictness, genefinding_tool, debug
        )
        local_runner.analyze_genomes()
        parser = smBGCParser(out_dir, cores)
        parser.gather_results()
        print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")


    if args["module"] == "download":
        out_dir = Path(args["out"])
        if not checks.check_if_dir_exists(out_dir):
            print(f"{out_dir} does not exist")
            return
        taxon = args["taxon"]
        debug = args["debug"]
        assembly_source = args["source"]
        assembly_level = args["level"]
        keep_plasmids = args["keep_plasmids"]
        keep_archive = args["keep_download"]
        downloader = GenomeDownloader(taxon=taxon, 
                                      out_dir=out_dir, 
                                      assembly_level=assembly_level, 
                                      assembly_source=assembly_source, 
                                      keep_plasmids=keep_plasmids
                                      debug=debug
                                      )
        print(f"Starting download: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        downloader.download_hydrated()
        downloader.extract_dataset_zip()
        downloader.create_output_directories()
        downloader.process_gbff_files()
        downloader.write_annotations()
        downloader.write_16S_fasta()
        if not keep_archive:
            downloder.remove_dataset_archive()
        print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
