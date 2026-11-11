# Excel row-limit TSV fallback — design

## Problem

Several pyPGCF modules write tabular results to `.xlsx` via `DataFrame.to_excel`
or `ExcelWriter`. Excel workbooks have a hard row limit (1,048,576 rows) and
become unreliable to open/edit well before that. For analyses over large
genomes or genera, some of these tables (functional-annotation hit tables,
pangenome core-protein membership tables, per-genome summary tables) can
realistically approach or exceed that range. When a table is large enough,
pyPGCF should write a `.tsv` instead of an `.xlsx`, and anything downstream
that reads the file back needs to transparently handle either format.

## Inventory of affected call sites

**Leaf outputs** (nothing in the pipeline reads these back):

| File | Output | Row driver |
|---|---|---|
| `cazy.py` | `CAZY_families.xlsx` | genomes |
| `amr.py` | `AMRfinder_results.xlsx` | AMR gene hits across all genomes |
| `smbgc.py` | `antiSMASH_results.xlsx` | BGC regions found |
| `virulence.py` | `VF_results_50pident50qcov.xlsx` | VF gene hits across all genomes |
| `eggnog.py` (`write_counts_dictionary_to_excel`) | `COG_categories_of_proteomes.xlsx` (2 sheets: counts, percent) | proteomes analyzed |
| `eggnog.py` (`write_hypergeometric_dfs`) | `eggNOG_hypergeometric.xlsx` (1-2 sheets: core, fingerprint) | proteomes analyzed |

**Non-leaf outputs** (read back elsewhere in the pipeline — this is why the
fix is a shared interface change, not a per-file tweak):

| Writer | Output | Reader |
|---|---|---|
| `core.py` | `{ref}_core.xlsx` / `{ref}_species_core.xlsx` (rows = every protein in the pangenome) | `eggnog.py`'s `EggNOGParser.read_core_protein_tables`, via the eggnog `-in`/`-in_list` CLI args |
| `download_genomes.py` | `Genome_information.xlsx` (rows = genomes) | `cli.py`'s species-demarcation flow (join + rewrite) |
| `species_demarcation.py` | `FastANI_species_clusters.xlsx` (rows = genomes) | same join in `cli.py` |
| `workflow.py` (`write_species_file_for_core`) | one-column temp file (`Group`, indexed by genome) | `cli.py`'s `run_core`, via the `-species` CLI arg |

## Threshold

**50,000 rows.** Applied uniformly to every call site above via one shared
helper, including the currently low-risk ones (cazy, smbgc, eggnog counts),
for consistency and future-proofing even though their row counts are
currently bounded by genome/proteome count rather than protein count.

## Shared helpers (`pypgcf/utils.py`)

```python
def write_table(
    df: DataFrame,
    fout: Path,
    *,
    index: bool = True,
    max_rows: int = 50_000,
    **to_kwargs,
) -> Path:
    """Write df to fout (.xlsx) unless it exceeds max_rows, in which case
    write a sibling .tsv instead. Returns the path actually written."""
```

```python
def write_table_sheets(
    sheets: dict[str, DataFrame],
    fout: Path,
    *,
    index: bool = True,
    max_rows: int = 50_000,
    styler_fn: Union[Callable[[DataFrame], "Styler"], None] = None,
) -> list[Path]:
    """Write all sheets to one .xlsx workbook (applying styler_fn per sheet,
    if given) when every sheet is within max_rows. If ANY sheet exceeds
    max_rows, write every sheet instead as a separate
    <fout.stem>_<sheetname>.tsv file (plain, unstyled) so the workbook never
    ends up split between formats. Returns the list of paths written."""
```

```python
def read_table(path: Path, *, index_col: Union[int, None] = None) -> DataFrame:
    """Read path as .tsv (tab-separated) or .xlsx based on its actual
    extension. If path doesn't exist, probe path.with_suffix(".xlsx") then
    path.with_suffix(".tsv") before raising FileNotFoundError. Covers both
    a user-supplied exact filename and internal code that hardcodes a
    .xlsx path whose real file may have been written as .tsv."""
```

Rationale for `read_table`'s two modes: within a single process (e.g.
`workflow.py` writing a temp file and `cli.py` reading
`workflow_runner.core_species_file` right after), the exact path is known
from `write_table`'s return value, so no probing is needed. Across separate
CLI invocations or user-supplied paths (eggnog's `-in`/`-in_list`, core's
`-species`), no in-memory state survives, so the reader must detect the
format from disk.

## Known trade-off: lost cell coloring on the hypergeometric fallback

`eggnog.py`'s `write_hypergeometric_dfs` currently colors over/under
-represented COG categories via
`df.style.apply(highlight_pvalue_on_output)` before writing to Excel. TSV
cannot carry cell styling. If that table ever exceeds 50,000 rows (which
would require 50,000 reference proteomes in a single run — essentially
theoretical at today's usage), the fallback `.tsv` keeps the numeric
p-value/fold-change columns but silently drops the coloring. No compensating
"significant" text column is added — this edge case is not worth the added
complexity today, but is worth revisiting if usage patterns change.

## Per-file changes

- `cazy.py`, `amr.py`, `smbgc.py`, `virulence.py`, `core.py`,
  `download_genomes.py`, `species_demarcation.py`: replace `df.to_excel(...)`
  with `write_table(df, ...)`, forwarding existing kwargs (`index=False`,
  `na_rep="X"`, etc.).
- `eggnog.py`:
  - `write_counts_dictionary_to_excel` → `write_table_sheets({"counts": df,
    "percent": df_perc}, fout)`.
  - `write_hypergeometric_dfs` → `write_table_sheets(hypergeom_results, fout,
    styler_fn=lambda df: df.style.apply(self.highlight_pvalue_on_output,
    axis=1))`.
  - `read_core_protein_tables` → `read_table(table_f, index_col=0)` in place
    of `read_excel(table_f, index_col=0)`.
- `workflow.py`'s `write_species_file_for_core`: capture `write_table`'s
  return value and assign it to `self.core_species_file` (instead of
  assuming `.xlsx`).
- `cli.py`:
  - `run_core`'s species-file read: `read_table(species_file, index_col=0)`.
  - The species-demarcation join-and-rewrite of `Genome_information.xlsx` /
    `FastANI_species_clusters.xlsx`: `read_table` for both reads,
    `write_table` for the rewrite.
  - Help text for `-in`/`-in_list` (eggnog subparser) and `-species` (core
    subparser) updated to mention `.xlsx` or `.tsv`.

## Testing

- `test_utils.py`: new tests for `write_table` (under/over threshold
  branches, kwarg forwarding), `write_table_sheets` (all-under-threshold →
  one workbook; any-over-threshold → per-sheet TSVs; `styler_fn` applied only
  on the xlsx path), and `read_table` (reads by actual suffix; probes
  `.xlsx` then `.tsv` when the given path doesn't exist; raises
  `FileNotFoundError` when neither exists).
- Update existing tests that assert `to_excel`/`read_excel` was called
  directly, in: `test_cazy.py`, `test_amr.py`, `test_smbgc.py`,
  `test_virulence.py`, `test_core.py`, `test_download_genomes.py`,
  `test_species_demarcation.py`, `test_workflow.py`, `test_eggnog.py`,
  `test_cli_runs.py`.

## Out of scope

- `workflow.py`'s parameter workbook (`self.param_file`, read via
  `read_excel` in `load_parameters`) — this is a small user-authored config
  sheet, not a results table, and isn't a candidate for a row-count-based
  format switch.
- Re-encoding the hypergeometric significance coloring as a text column for
  the TSV fallback path (see trade-off above).
