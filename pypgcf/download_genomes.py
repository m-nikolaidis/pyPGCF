from os import system
from pathlib import Path
from zipfile import ZipFile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pandas import DataFrame

from tqdm import tqdm

from pypgcf import dispatchers, utils, checks

FILETYPE_TO_DIR = {
    "cds": "CDS_fasta_files",
    "protein": "Protein_fasta_files",
    "genome": "Genomic_fasta_files",
}

GBFF_INFO_COLUMNS = ["organism", "strain", "isolation_source", "host", "note", "country", "collection_date", "Perc_pseudogenes", "Genome_len", "N%", "GC%", "Num_plasmids"]

class GenomeDownloader:
    """
    Class tasked with downloading the gbff files from a given taxon
    This class will also extract the various gene, protein and genomic sequences from the gbff file
    """
    # TODO: Need to add option to download specific GCF IDs
    def __init__(
        self,
        taxon: str,
        out_dir: Path,
        assembly_level: str = "chromosome,complete",
        assembly_source: str = "RefSeq",
        keep_plasmids: bool = True,
        debug: bool = False,
    ):
        self.taxon = taxon
        self.assembly_level = assembly_level
        self.out_dir = out_dir
        self.assembly_source = assembly_source
        self.keep_plasmids = keep_plasmids
        self.debug = debug
        self.gbff_source_info = {}
        self.plasmids_per_genome = {} 
        self.rrna_16S_records = {}

        if not checks.is_valid_assembly_source(self.assembly_source):
            raise ValueError(f"{self.assembly_source} is not a valid value. Valid values: RefSeq or GenBank")

    def download_hydrated(self):
        cmd = " ".join(
            [
                "datasets download genome taxon",
                str(self.taxon),
                "--include gbff",
                f"--assembly-level {self.assembly_level}",
                f"--assembly-source {self.assembly_source} --filename {self.out_dir}/dataset.zip",
            ]
        )
        if self.debug:
            cmd += " --preview"
        res = dispatchers.execute_command(cmd)
        return res

    def download_from_file(self):
        ...

    def extract_dataset_zip(self) -> None:
        filename = self.out_dir / "dataset.zip"
        if not checks.check_if_file_exists(filename):
            raise FileNotFoundError(f"{filename} doesn't exist")
        with ZipFile(filename, 'r') as zf:
            zf.extractall(self.out_dir)
        return None

    def create_output_directories(self) -> None:
        for filetype, outdir_suffix in FILETYPE_TO_DIR.items():
            outdir = self.out_dir / outdir_suffix
            outdir.mkdir(exist_ok=True, parents=True)
        return None


    def save_gbff_dir_to_memory(self) -> None:
        gbff_files: dict = {}
        indir = self.out_dir / "ncbi_dataset" / "data"
        subdirs = list(indir.glob("GC*"))
        for subdir in subdirs:
            genome = subdir.stem
            gbff_files[genome] = subdir / "genomic.gbff"
        self.gbff_files = gbff_files


    def gbff_file_to_genomic_fasta(self, records: list, fout: Path):
        SeqIO.write(records, fout, "fasta")
    
    def gbff_file_to_protein_fasta(self, records: list, fout: Path):
        proteins = []
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                locus_tag = feature.qualifiers.get("locus_tag", ["X"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["X"])[0]
                translation = feature.qualifiers.get("translation", ["X"])[0]
                product = feature.qualifiers.get("product", ["X"])[0]
                if translation == "X" and "pseudo" in feature.qualifiers: # Skip pseudogenes
                    continue
                seq = Seq(translation)
                description = protein_id + " " + product
                item = SeqRecord(name=locus_tag, id=locus_tag, description=description, seq=seq)
                proteins.append(item)
        SeqIO.write(proteins, fout, "fasta")
    
    def gbff_file_to_cds_fasta(self, records: list, fout: Path):
        cds_records = []
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                locus_tag = feature.qualifiers.get("locus_tag", ["X"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["X"])[0]
                translation = feature.qualifiers.get("translation", ["X"])[0]
                seq = feature.extract(record.seq)
                product = feature.qualifiers.get("product", ["X"])[0]
                if translation == "X" and "pseudo" in feature.qualifiers: # Skip pseudogenes
                    continue
                description = protein_id + " " + product
                item = SeqRecord(name=locus_tag, id=locus_tag, description=description, seq=seq)
                cds_records.append(item)
        SeqIO.write(cds_records, fout, "fasta")

    def gather_gbff_source_info(self, records: list) -> dict:
        data = {}
        rec = records[0] # Do not need to use other records, all have the same source qualifiers
        for field_name, field_value in rec.features[0].qualifiers.items():
            data[field_name] = ";".join(field_value)
        return data

    def get_16S_from_gbff(self, records: list) -> SeqRecord:
        data = [None, "LocusTag", None, 0]
        for record in records:
            features = record.features
            for feature in features:
                if feature.type != "rRNA":
                    continue
                product = feature.qualifiers.get("product", [None])[0]
                if product != "16S ribosomal RNA":
                    continue
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                feature_length = len(feature.location)
                feature_coordinates = feature.location
                if feature_length > data[3]:
                    data = [record, locus_tag, feature_coordinates, feature_length]

        record, locus_tag, location, _ = data
        if location is None:
            return SeqRecord()
        if location.strand == -1:
            seq = location.extract(record.seq).reverse_complement()
        else:
            seq = location.extract(record.seq)
        rna_record = SeqRecord(name=locus_tag, id=locus_tag, seq=seq)
        return rna_record

    
    def calculate_perc_pseudogenes(self, records: list) -> float:
        num_pseudogenes = 0
        num_genes = 0
        perc_pseudogenes = 100.0
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                num_genes += 1
                locus_tag = feature.qualifiers.get("locus_tag", ["X"])[0]
                protein_id = feature.qualifiers.get("protein_id", ["X"])[0]
                translation = feature.qualifiers.get("translation", ["X"])[0]
                product = feature.qualifiers.get("product", ["X"])[0]
                if translation == "X" and "pseudo" in feature.qualifiers: # Skip pseudogenes
                    num_pseudogenes += 1
                    continue
                if num_genes != 0:
                    perc_pseudogenes = round(num_pseudogenes / num_genes * 100, 2)
        return perc_pseudogenes
        # TODO: There is too much repeating code.
        # Should read from disk only once to create the parserand
        # These functions should accept the SeqRecord as input and return lists
        # The IO should be put to separate functions (one for writing fasta, and one for reading)
        # Maybe reading the whole file into memory is better.
        # Would never encounter huge files that do not fit in memory
        # Would never open too many files and cause a memory throttle

    def calculate_GC_and_N(self, records: list) -> tuple[int, float, float]:
        total_unknown_nucl = 0
        total_gc = 0
        total_genome_len = 0
        total_unknown_perc = 0.0
        total_gc_perc = 0.0
        for record in records:
            seq = record.seq
            if seq is not None:
                seqlen: int = len(record.seq)
                g_count: int = record.seq.count("G")
                c_count: int = record.seq.count("C")
                a_count: int = record.seq.count("A")
                t_count: int = record.seq.count("T")
                known_nucl_count = g_count + c_count + a_count + t_count
                total_unknown_nucl += seqlen - known_nucl_count
                total_gc += g_count
                total_gc += c_count
                total_genome_len += seqlen
        total_unknown_perc = round(total_unknown_nucl / total_genome_len * 100, 2)
        total_gc_perc = round(total_gc / total_genome_len * 100, 3)
        return total_genome_len, total_unknown_perc, total_gc_perc
    
    def identify_plasmids_from_gbff(self, records: list) -> list:
        plasmids = []
        for record in records:
            has_plasmid_in_desc = False
            has_plasmid_in_source = False
            if "plasmid" in record.description:
                has_plasmid_in_desc = True
            source = record.features[0]
            plasmid_field = source.qualifiers.get("plasmid")
            if plasmid_field is not None:
                has_plasmid_in_source = True
            if has_plasmid_in_source or has_plasmid_in_desc:
                plasmids.append(record.name)
        return plasmids

    def process_gbff_files(self) -> None:
        for genome, file in tqdm(self.gbff_files.items(), desc="Extracting information"):
            parser = SeqIO.parse(file, "genbank")
            records = [r for r in parser]
            plasmids = self.identify_plasmids_from_gbff(records)
            if not self.keep_plasmids:
                records = [r for r in records if r.name not in plasmids]
            genome_fout = self.out_dir /  FILETYPE_TO_DIR.get("genome") / (genome + ".fna")
            protein_fout = self.out_dir /  FILETYPE_TO_DIR.get("protein") / (genome + ".faa")
            cds_fout = self.out_dir /  FILETYPE_TO_DIR.get("cds") / (genome + ".fna")
            self.gbff_file_to_genomic_fasta(records, genome_fout)
            self.gbff_file_to_protein_fasta(records, protein_fout)
            self.gbff_file_to_cds_fasta(records, cds_fout)
            rrna_16S_record = self.get_16S_from_gbff(records)
            tmp_info = self.gather_gbff_source_info(records)
            perc_pseudogenes = self.calculate_perc_pseudogenes(records)
            genome_len, unknown_nucl, gc_perc = self.calculate_GC_and_N(records)
            self.gbff_source_info[genome] = tmp_info
            self.gbff_source_info[genome]["Perc_pseudogenes"] = perc_pseudogenes
            self.gbff_source_info[genome]["Genome_len"] = genome_len
            self.gbff_source_info[genome]["N%"] = unknown_nucl
            self.gbff_source_info[genome]["GC%"] = gc_perc
            self.gbff_source_info[genome]["Num_plasmids"] = len(plasmids)
            self.rrna_16S_records[genome] = rrna_16S_record
            # self.plasmids_per_genome[genome] = plasmids
        return None

    def write_annotations(self) -> None:
        df = DataFrame.from_dict(self.gbff_source_info, orient="index")
        fout = self.out_dir / "Genome_information.xlsx"
        # TODO: Assign problematic assemblies
        # TODO: Keep only certain columns if they exist in the dataframe
        if not df.empty:
            for column in GBFF_INFO_COLUMNS:
                if column in df.columns:
                    continue
                df[column] = "X" # Assign the missing columns
            df = df[GBFF_INFO_COLUMNS]
        df.to_excel(fout, na_rep="X")
        return None

    def write_16S_fasta(self):
        records = []
        for genome, record in self.rrna_16S_records.items():
            record.description = record.name[:]
            record.id = genome
            records.append(record)
        fout = self.out_dir / "16S_sequences.fna"
        SeqIO.write(records, fout, "fasta")

    def remove_dataset_archive(self):
        to_remove = [
            self.out_dir / "dataset.zip",
            self.out_dir / "README.md",
            self.out_dir / "ncbi_dataset"
        ]
        for item in to_remove:
            utils.recursive_unlink(item)
