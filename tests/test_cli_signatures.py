"""
Every run_<module> function must call its analysis classes with arguments
that the real constructors accept. The other CLI tests mock these classes,
so a renamed keyword would otherwise only fail on a real run.
"""

import argparse
import inspect
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from pypgcf import cli

CLASSES = [
    "SpeciesDemarcator", "Orthologues_identifier", "Core_identifier",
    "Phylogenomic", "EggNOGRunner", "EggNOGParser", "smBGCLocalRunner",
    "smBGCParser", "VF_analyzer", "CAZY_analyzer", "AMR_analyzer",
    "GenomeDownloader", "SMBGC_installer", "EGGNOG_installer",
    "CAZY_installer", "AMR_installer", "VF_installer",
]
MODULES = [
    "species_demarcation", "orthologues", "core", "phylogenomic", "eggnog",
    "smbgc", "virulence", "cazy", "amr", "download", "databases",
]


def _parse_minimal_args(module: str) -> dict:
    """Parse the module's required options with placeholder values."""
    parser = cli.setup_parser()
    subparsers = next(
        a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
    )
    subparser = subparsers.choices[module]
    required = [a for a in subparser._actions if a.required and a.option_strings]
    for group in subparser._mutually_exclusive_groups:
        if group.required:
            required.append(group._group_actions[0])
    argv = [module]
    for action in required:
        value = action.choices[0] if action.choices else "x"
        if action.dest == "install":
            value = "all"  # construct every database installer
        argv += [action.option_strings[0], str(value)]
    args = vars(parser.parse_args(argv))
    args["fasta_files"] = [Path("a.fa")]
    return args


class TestRunFunctionSignatures(unittest.TestCase):
    def test_constructor_calls_match_real_signatures(self):
        real = {name: getattr(cli, name) for name in CLASSES}
        for module in MODULES:
            calls = []

            def recorder(name):
                def record(*args, **kwargs):
                    calls.append((name, args, kwargs))
                    return MagicMock()

                return MagicMock(side_effect=record)

            mocks = {name: recorder(name) for name in CLASSES}
            with patch.multiple(cli, **mocks), patch.object(
                cli, "validate_directory", return_value=True
            ), patch.object(cli, "read_excel"), patch.object(
                cli.utils, "join_pandas_dataframes"
            ):
                getattr(cli, f"run_{module}")(_parse_minimal_args(module))

            self.assertTrue(calls, f"run_{module} constructed no classes")
            for name, args, kwargs in calls:
                with self.subTest(module=module, cls=name):
                    inspect.signature(real[name]).bind(*args, **kwargs)


if __name__ == "__main__":
    unittest.main()
