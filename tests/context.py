import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pypgcf.cazy as cazy
import pypgcf.cli as cli
import pypgcf.config as config
import pypgcf.core as core
import pypgcf.databases as databases
import pypgcf.download_genomes as download_genomes
import pypgcf.eggnog as eggnog
import pypgcf.orthologues as orthologues
import pypgcf.phylogenomic as phylogenomic
import pypgcf.smbgc as smbgc
import pypgcf.species_demarcation as species_demarcation
import pypgcf.utils as utils
import pypgcf.virulence as virulence
import pypgcf.workflow as workflow
import pypgcf.amr as amr

resources_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../resources"))
