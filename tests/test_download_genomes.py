
from dataclasses import dataclass
from pathlib import Path
import unittest
import pytest
from pypgcf import download_genomes, utils
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

def calculate_num_records_in_fasta(file: Path) -> int:
    num = 0
    parser = parse(file, "fasta")
    for _ in parser:
        num += 1
    return num


@dataclass
class Data:
    taxon = 1234
    taxon_str = "Nitrospira"
    outdir = Path("../test_data/")
    assembly_source = "RefSeq"
    assembly_levels = "chromosome,complete"
    debug = True

# TODO: Implement test with wrong assembly source
# TODO: Implement test with wrong assembly levels

class TestModule(unittest.TestCase):

    def test_download_hydrated(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        res = downloader.download_hydrated()
        # downloaded_file = Path(data.outdir / "dataset.zip")
        self.assertIs(res, 0)
        # self.assertIs(downloaded_file.exists(), True)
        
    def test_create_output_directories(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        downloader.create_output_directories()
        expected_directories = [
            data.outdir / "Genomic_fasta_files",
            data.outdir / "Protein_fasta_files",
            data.outdir / "CDS_fasta_files",
        ]
        unexpected_directory = data.outdir / "GBFF_files"
        for expected_directory in expected_directories:
            self.assertIs(expected_directory.exists(), True)
            expected_directory.rmdir()
        self.assertIs(unexpected_directory.exists(), False)
        # TODO: When I add option to create either protein or CDS files I should change this test

    def test_extract_dataset_zip(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        expected_directory = Path(downloader.out_dir / "ncbi_dataset")
        readme_file = Path(downloader.out_dir / "README.md")
        self.assertIs(expected_directory.exists(), False)
        downloader.extract_dataset_zip()
        self.assertIs(expected_directory.exists(), True)
        utils.recursive_unlink(expected_directory)
        utils.recursive_unlink(readme_file)


    def test_save_gbff_to_memory(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        expected_directory = Path(downloader.out_dir / "ncbi_dataset")
        readme_file = Path(downloader.out_dir / "README.md")
        downloader.extract_dataset_zip()
        downloader.save_gbff_dir_to_memory()
        self.assertIs(type(downloader.gbff_files), dict)
        self.assertIs(len(downloader.gbff_files.keys()), 2)

        # Cleanup
        utils.recursive_unlink(expected_directory)
        utils.recursive_unlink(readme_file)

    def test_process_gbff_files_with_plasmids(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            keep_plasmids=True,
            debug=data.debug,
        )
        downloader.gbff_files = {
            "GCF_022870945": Path("../test_data/GCF_022870945.gbff"), 
            "GCF_902375805": Path("../test_data/GCF_902375805.gbff")
        }
        downloader.create_output_directories()
        downloader.process_gbff_files()

        expected_directories = [
            data.outdir / "Genomic_fasta_files",
            data.outdir / "Protein_fasta_files",
            data.outdir / "CDS_fasta_files",
        ]
        suffixes = [".fna", ".faa", ".fna"]
        num_records = {
            "Genomic_fasta_files": {
                "GCF_022870945": 3,
                "GCF_902375805": 10,

            },
            "Protein_fasta_files": {
                "GCF_022870945": 2467,
                "GCF_902375805": 3563,
            },
            "CDS_fasta_files": {
                "GCF_022870945": 2467,
                "GCF_902375805": 3563,
            },

        }
        for expected_directory, suffix in zip(expected_directories, suffixes):
            g1 = "GCF_022870945"
            g2 = "GCF_902375805"
            f1 = expected_directory / (g1 + suffix)
            f2 = expected_directory / (g2 + suffix)
            self.assertIs(f1.exists(), True)
            self.assertIs(f2.exists(), True)
            num_f1_records = calculate_num_records_in_fasta(f1)
            num_f2_records = calculate_num_records_in_fasta(f2)
            exp_num_f1_records = num_records.get(expected_directory.stem, {}).get(g1)
            exp_num_f2_records = num_records.get(expected_directory.stem, {}).get(g2)
            if exp_num_f1_records is None or exp_num_f2_records is None:
                raise ValueError("Not possible")
            self.assertEqual(num_f1_records, exp_num_f1_records)
            self.assertEqual(num_f2_records, exp_num_f2_records)
            utils.recursive_unlink(expected_directory)

    def test_process_gbff_files_without_plasmids(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            keep_plasmids=False,
            debug=data.debug,
        )
        downloader.gbff_files = {
            "GCF_022870945": Path("../test_data/GCF_022870945.gbff"), 
            "GCF_902375805": Path("../test_data/GCF_902375805.gbff")
        }
        downloader.create_output_directories()
        downloader.process_gbff_files()

        expected_directories = [
            data.outdir / "Genomic_fasta_files",
            data.outdir / "Protein_fasta_files",
            data.outdir / "CDS_fasta_files",
        ]
        suffixes = [".fna", ".faa", ".fna"]
        num_records = {
            "Genomic_fasta_files": {
                "GCF_022870945": 1,
                "GCF_902375805": 10,

            },
            "Protein_fasta_files": {
                "GCF_022870945": 2440,
                "GCF_902375805": 3563,
            },
            "CDS_fasta_files": {
                "GCF_022870945": 2440,
                "GCF_902375805": 3563,
            },

        }
        for expected_directory, suffix in zip(expected_directories, suffixes):
            g1 = "GCF_022870945"
            g2 = "GCF_902375805"
            f1 = expected_directory / (g1 + suffix)
            f2 = expected_directory / (g2 + suffix)
            self.assertIs(f1.exists(), True)
            self.assertIs(f2.exists(), True)
            num_f1_records = calculate_num_records_in_fasta(f1)
            num_f2_records = calculate_num_records_in_fasta(f2)
            exp_num_f1_records = num_records.get(expected_directory.stem, {}).get(g1)
            exp_num_f2_records = num_records.get(expected_directory.stem, {}).get(g2)
            if exp_num_f1_records is None or exp_num_f2_records is None:
                raise ValueError("Not possible")
            self.assertEqual(num_f1_records, exp_num_f1_records)
            self.assertEqual(num_f2_records, exp_num_f2_records)
            utils.recursive_unlink(expected_directory)

    def test_calculate_GC_and_N_base(self):
        records = [
            SeqRecord(name="1", seq=Seq("AGCT")),
            SeqRecord(name="2", seq=Seq("CGAT")),
            SeqRecord(name="3", seq=Seq("NNNN"))
        ]
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        result = downloader.calculate_GC_and_N(records)
        self.assertEqual(result, (12, 33.33, 33.33333))

    def test_calculate_GC_N_no_sequence(self):
        records = []
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        result = downloader.calculate_GC_and_N(records)
        self.assertEqual(result, (0, 0.0, 0.0))

    def test_calculate_GC_and_N_all_known_nucleotides(self):
        records = [
            SeqRecord(name="1", seq=Seq("AAA")),
            SeqRecord(name="2", seq=Seq("CCC")),
            SeqRecord(name="3", seq=Seq("GGG")),
            SeqRecord(name="4", seq=Seq("TTT"))
        ]
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        result = downloader.calculate_GC_and_N(records)
        self.assertEqual(result, (12, 0.0, 50.0))

    def test_calculate_GC_and_N_all_unknown_nucleotides(self):
        records = [
            SeqRecord(name="1", seq=Seq("NNN")),
            SeqRecord(name="1", seq=Seq("NNNNNN")),
        ]
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        result = downloader.calculate_GC_and_N(records)
        self.assertEqual(result, (9, 100.0, 0.0))

    def test_get_16S_from_gbff(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        downloader.gbff_files = {
            "GCF_022870945": Path("../test_data/GCF_022870945.gbff"), 
            "GCF_902375805": Path("../test_data/GCF_902375805.gbff")
        }
        observed_sequences = {}
        for genome, gbff_file in downloader.gbff_files.items():
            parser = parse(gbff_file, "genbank")
            records = [r for r in parser]
            rrna_16S_record = downloader.get_16S_from_gbff(records)
            observed_sequences[genome] = [rrna_16S_record.name, len(rrna_16S_record.seq)]

        # def get_16S_from_gbff(self, records: list) -> SeqRecord:
        expected_sequences = {
            "GCF_022870945": ["LVJ81_RS08520", 1542],
            "GCF_902375805": ["FXY68_RS00045", 1538]
        }

        for key in expected_sequences:
            exp_name = expected_sequences.get(key)[0]
            obs_name = observed_sequences.get(key)[0]
            exp_length = expected_sequences.get(key)[1]
            obs_length = observed_sequences.get(key)[1]
            self.assertEqual(exp_name, obs_name)
            self.assertEqual(exp_length, obs_length)

    # def test_write_annotations(self):
    #     data = Data()
    #     downloader = download_genomes.GenomeDownloader(
    #         taxon=data.taxon,
    #         out_dir=data.outdir,
    #         debug=data.debug,
    #     )
    #     downloader.gbff_files = {
    #         "GCF_022870945": Path("../test_data/GCF_022870945.gbff"), 
    #         "GCF_902375805": Path("../test_data/GCF_902375805.gbff")
    #     }
    #
    #     self.gbff_source_info = {'GCF_022870945': {'organism': 'Vitreoscilla stercoraria', 'mol_type': 'genomic DNA', 'strain': 'SAG 1488-6', 'isolation_source': 'Dung', 'host': 'Bos taurus', 'culture_collection': 'DSM:513', 'type_material': 'type strain of Vitreoscilla stercoraria', 'db_xref': 'taxon:61', 'collection_date': '1990', 'Perc_pseudogenes': 0, 'Genome_len': 1000000, 'N%': 0, 'GC%': 55}, 'GCF_902375805': {'organism': 'Vitreoscilla massiliensis', 'mol_type': 'genomic DNA', 'isolate': 'MGYG-HGUT-01520', 'isolation_source': 'human gut', 'db_xref': 'taxon:1689272', 'note': 'contig: NZ_LN898212.1', 'Perc_pseudogenes': 0, 'Genome_len': 1000000, 'N%': 0, 'GC%': 55}}

    def test_write_annotations_empty_input(self):
        data = Data()
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=data.outdir,
            debug=data.debug,
        )
        downloader.gbff_files = {
            "GCF_022870945": Path("../test_data/GCF_022870945.gbff"), 
            "GCF_902375805": Path("../test_data/GCF_902375805.gbff")
        }
        downloader.write_annotations()
        expected_file = data.outdir / "Genome_information.xlsx"
        self.assertIs(expected_file.exists(), True)
        df = pd.read_excel(expected_file, index_col=0)
        self.assertIs(df.empty, True)
        

    def test_write_annotations_invalid_input(self): # Is this possible?
        ...

    @pytest.mark.skip(reason="Gives too many requests error")
    def test_pipeline(self):
        data = Data()
        outdir = data.outdir / "tmp"
        outdir.mkdir(exist_ok=True)
        downloader = download_genomes.GenomeDownloader(
            taxon=data.taxon,
            out_dir=outdir,
            debug=False,
        )
        downloader.download_hydrated()
        downloader.extract_dataset_zip()
        downloader.create_output_directories()
        downloader.save_gbff_dir_to_memory()
        downloader.process_gbff_files()
        downloader.write_annotation()
        downloader.write_16S_fasta()

    # def test_gbff_file_to_genomic_fasta
    #     data = Data()
    #     downloader = download_genomes.GenomeDownloader(
    #         taxon=data.taxon,
    #         assembly_levels=data.assembly_levels,
    #         out_dir=data.outdir,
    #         assembly_source=data.assembly_source,
    #         debug=data.debug,
    #     )
    #     downloader

        # gbff_(file: Path, outdir: Path):
 

    # def test_download_hydrated(self):
    #     data = Data()
    #     downloader = download_genomes.GenomeDownloader(
    #         taxon=data.taxon,
    #         assembly_levels=data.assembly_levels,
    #         out_dir=data.outdir,
    #         filetypes=data.filetypes,
    #         assembly_source=data.assembly_source,
    #         debug=data.debug,
    #     )
    #     res = downloader.download_hydrated()
    #     # downloaded_file = Path(data.outdir / "dataset.zip")
    #     self.assertIs(res, 0)
    #     # self.assertIs(downloaded_file.exists(), True)
