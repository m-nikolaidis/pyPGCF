from pathlib import Path
from unittest.mock import patch

import pytest
from pandas import DataFrame

from pypgcf.amr import AMR_analyzer


@pytest.mark.parametrize("type_column", ["Type", "Element type"])
def test_parse_results_supports_new_and_legacy_headers(tmp_path, type_column):
    analyzer = AMR_analyzer(
        fasta_files_list=[],
        database_dir=None,
        out_dir=tmp_path,
        input_type="prot",
        available_cores=1,
        blast_cores=1,
        num_concurrent=1,
        debug=False,
    )
    results_dir = analyzer.out_dir / "amrfinder"
    results_dir.mkdir(parents=True)
    (results_dir / "genome.txt").write_text(
        f"Protein id\tElement symbol\t{type_column}\tClass\n"
        "protein1\tarsC\tSTRESS\tARSENIC\n"
        "protein2\tfosB\tAMR\tFOSFOMYCIN\n"
        "protein3\tvirulence\tVIRULENCE\tNA\n"
    )

    with patch.object(DataFrame, "to_excel", autospec=True) as write_excel:
        analyzer.parse_results()

    result, destination = write_excel.call_args.args
    assert list(result.index) == ["protein2"]
    assert result["Element symbol"].tolist() == ["fosB"]
    assert result["Genome"].tolist() == ["genome"]
    assert destination == analyzer.out_dir / "AMRfinder_results.xlsx"
