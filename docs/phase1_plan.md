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

### 1. Remove the duplicate orchestrator — DONE
- Deleted `pypgcf/pipeline.py`. `workflow.py` is the only orchestrator. The test suite passes without it.

### 2. Declare the real dependencies — DONE
- `pyproject.toml` lists the Python runtime dependencies (pandas, openpyxl, tqdm, biopython, requests, scipy) and `requires-python = ">=3.11,<3.13"`.
- `pixi.toml` and `recipe.yaml` match that list, and also declare BLAST and DIAMOND, which were missing. FastTree was dropped: phylogenomic now uses `iqtree --fast` (`pyPGCF phylogenomic --fast`), so no separate FastTree dependency is needed.
- Version is set once in `pyproject.toml` (2.0.1). `pypgcf.__version__` reads it from package metadata. `recipe.yaml` had 2.0.1.2 and is now 2.0.1.
- `tests/test_packaging.py` fails if the version, Python range or runtime dependencies drift between the three files.
- Verified: a clean export of HEAD installs with `pip install .` on Python 3.11 and 3.12 and `import pypgcf.cli` works. The full suite passes.
- Still open: `pixi.lock` needs regenerating with `pixi lock` (pixi was not available here).

### 3. Triage and fix the failing tests — DONE
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

### 4. Replace `os.system` with `subprocess`; add a tool preflight — DONE
- `utils.execute_command` uses `subprocess.run`. A string still runs through the shell; a list of arguments runs without one. It returns the exit code (not `os.system`'s wait status).
- Phylogenomic (MUSCLE, IQ-TREE) and orthologues (DIAMOND, BLAST) build argument lists.
- `config.MODULE_TOOLS` lists the binaries per module and `config.TOOL_PACKAGES` the conda package for each. `utils.check_tools` reports every missing binary with the `mamba install` command.
- `pyPGCF <module>` checks its tools before any work starts. A workflow checks the tools for all selected tasks before gathering input.
- Bugs fixed on the way: CDS orthologues built `blastn makedb` and `blastn blastp` commands (wrong input-type checks), and orthologues `--input_type` accepted any string, so `CDS` silently ran DIAMOND on nucleotides.
- Still on shell strings: species demarcation (MCL pipeline with redirection), eggNOG, smBGC, AMR, CAZy, virulence, download and database installers. They work unchanged; migrate them when those modules are next touched.

### 5. Split the dependency suite into core and full environments — DONE
- pixi environments: `default` (core: FastANI, MCL, BLAST, DIAMOND, MUSCLE, IQ-TREE, NCBI datasets), `full` (adds eggNOG-mapper, antiSMASH, dbCAN, AMRFinderPlus) and `test` (core plus pytest). One solve group.
- The NCBI datasets CLI is in core rather than full: it is one small binary and downloading genomes is the first pipeline step.
- `pixi.lock` is generated with pixi 0.62.2 to keep lock format v6; pixi 0.82 upgrades it to v7.
- Verified: the core environment installs and runs the suite with the real tools, locally and in CI.
- Found on the way: `pixi run pyPGCF` never ran the CLI (no `__main__` guard). `python -m pypgcf` now works and the task uses it.
- The conda recipe is still one full package (see open questions).

### 6. Explicit database paths — DONE
- `config.database_dir(module, cli_value)` is the only place a database path is resolved: the CLI option, then `PYPGCF_DB_DIR`, then the earlier default. It logs the result (visible with `--debug`).
- Kept the existing option names (`--db` per module, `-db` for workflow, `--database_dir` for installs) instead of adding `--db_dir`, so existing commands keep working.
- Defaults are unchanged, including eggNOG's legacy `site-packages/site-packages/data` path, so existing installs keep working.
- Bugs fixed on the way: CAZy and AMR read from `CAZY/CAZY` and `AMR/AMR` while installers wrote `CAZY` and `AMR`; `pyPGCF databases` without `--database_dir` crashed; `pyPGCF amr`, `pyPGCF eggnog` (after the mapper step) and `pyPGCF cazy` crashed with wrong constructor arguments or missing options; `pyPGCF amr` without `-out` crashed instead of printing usage.
- `tests/test_cli_signatures.py` binds every run function's constructor calls against the real signatures, so this class of crash fails the tests.

### 7. CI — DONE
- `.github/workflows/ci.yml` runs on every push and pull request:
  - `tests`: Python 3.11 and 3.12, flake8 for syntax errors and undefined names, pytest, `pyPGCF -h`.
  - `external-tools`: the pixi `test` environment, a check that it provides every core module tool, and the suite with the real tools.
  - `package`: builds `recipe/recipe.yaml` with rattler-build (installs it into a fresh environment and runs the recipe tests) and runs `pixi build`.
- All jobs green. Full flake8 and black checks are not enabled: the existing code has 569 flake8 findings and 18 files black would reformat (see open questions).

### 8. Bioconda packaging check — DONE
- Recipe moved to `recipe/recipe.yaml`. At the repository root, rattler-build copied the whole repository into the package (8.8 MB for 40 KB of code); now 49 KB.
- Added an `about` section (GPL-3.0-or-later, homepage, summary), `noarch: python` and build number 0. Tests run `pyPGCF -h`, `pyPGCF core --help` and `python -m pypgcf -h`.
- Verified with rattler-build locally and in CI: the package builds, installs into a fresh environment with all run dependencies, and passes its tests.
- Not done: submitting to Bioconda needs a tagged release with a source URL and checksum (see open questions).

## Sequence

- Week 1: items 1, 2, 3.
- Week 2: items 4, 5, 6.
- Week 3: items 7, 8.

## Open questions

1. Bioconda recipe: one full package (all modules), or separate core and full packages?
2. Bioconda submission needs a tagged GitHub release (for example v2.0.1) to point the recipe at. Create one?
3. Code style: apply black to the whole codebase? It touches 18 files, so it cannot be split into commits under 100 lines; one mechanical commit per file is possible.
