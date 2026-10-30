# pyPGCF Phase 1 plan: installable, testable and safe to run

Status as of 2026-10-18. Branch: `claude/phylogenomics-pipeline-upgrade-3v99i4`.

## Baseline (before any changes)

Measured in a Python 3.11 environment with pytest 9.1 and the runtime packages installed by hand.

- `pyproject.toml` declares no runtime `dependencies`. A plain `pip install .` produces a package that cannot import (pandas, scipy, biopython, tqdm, requests and openpyxl are all missing).
- The test suite only runs with `PYTHONPATH=tests`, because test modules use `from context import ...`.
- `tests/test_databases.py` cannot be imported (`from databases import ...`).
- Result: 136 passed, 35 failed, 20 skipped (test_databases.py excluded).
- `pypgcf/pipeline.py` duplicated `pypgcf/workflow.py`, and nothing imported it.
- `utils.execute_command` builds shell strings and runs them with `os.system`.
- Database location is derived from `sysconfig` (`DB_BASE_DIR`).
- Version numbers disagree: `pyproject.toml` 2.0.1, `recipe.yaml` 2.0.1.2, README install line 2.0.0.

## Work items

### 1. Remove the duplicate orchestrator — DONE (pending commit)
- Deleted `pypgcf/pipeline.py`. `workflow.py` is the only orchestrator. The test suite passes without it.

### 2. Declare the real dependencies — DONE (pending push)
- `pyproject.toml` lists the Python runtime dependencies (pandas, openpyxl, tqdm, biopython, requests, scipy) and `requires-python = ">=3.11,<3.13"`.
- `pixi.toml` and `recipe.yaml` match that list, and also declare BLAST, DIAMOND and FastTree, which were missing.
- Version is set once in `pyproject.toml` (2.0.1). `pypgcf.__version__` reads it from package metadata. `recipe.yaml` had 2.0.1.2 and is now 2.0.1.
- `tests/test_packaging.py` fails if the version, Python range or runtime dependencies drift between the three files.
- Verified: a clean export of HEAD installs with `pip install .` on Python 3.11 and 3.12 and `import pypgcf.cli` works. The full suite passes.
- Still open: `pixi.lock` needs regenerating with `pixi lock` (pixi was not available here).

### 3. Triage and fix the failing tests — DONE (pending commit)
- `pytest` runs from the repo root with no `PYTHONPATH` (`pythonpath = ["tests"]` in `pyproject.toml`).
- Result: 182 passed, 32 skipped, 0 failed (baseline before work: 136 passed, 35 failed, 20 skipped, with `test_databases.py` not collected).
- Mock patch targets rewritten from `context.<module>` / `databases.<name>` to the real `pypgcf.<module>` paths.
- Stale tests updated to current signatures and method names (CAZY, Phylogenomic, SpeciesDemarcator, smBGC runner, utils DB commands, CLI run fixtures, mkdir expectations).
- Markers registered: `external`, `data`, `integration`, `slow`. The NCBI download test skips unless the `datasets` CLI is installed.
- Skips with explicit reasons: tests that need `data/` or `../data/` (not in the repository), and the VF size-mismatch test (see open decisions).

#### Real bugs found and fixed during triage
- `phylogenomic.py`: `concurrent_jobs` was assigned only inside the `else` branch of the IQ-TREE model check. With the default model (`None`), the attribute was never set and the alignment step would fail with `AttributeError`.
- `cazy.py`: `_clean_dmnd_dbcansub_output` joined a `set`, so family order changed between Python processes (hash randomisation). Output is now sorted.

#### Decisions (answered 2026-10-18)

1. **VF size check:** keep warning-only. Remote sizes are not reliable, so a mismatch does not block installs. The VF size-mismatch test stays skipped with its reason.
2. **Test data:** commit a small fixture subset under `tests/data` (a few genomes plus expected outputs) so the `data`-marked tests run in CI. Not started yet.
3. **Python range:** widen to `>=3.11,<3.13`, and test both versions in CI. Applies to item 2 and the conda recipe.
4. **Excel parameter files:** stay supported through Phase 1. YAML config comes in Phase 2.
