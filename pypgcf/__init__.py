"""
pyPGCF: PhyloGenomic, Core and Fingerprint analyses package
author: Marios Nikolaidis
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _package_version

try:
    __version__ = _package_version("pypgcf")
except PackageNotFoundError:  # running from a source tree that is not installed
    __version__ = "0+unknown"
