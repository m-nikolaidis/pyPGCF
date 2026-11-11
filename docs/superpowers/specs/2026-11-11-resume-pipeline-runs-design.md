# Resumable module runs — design

## Problem

Several pyPGCF modules dispatch one external-tool invocation per input
fasta/genome file across a process pool (`multiprocess_dispatch` in
`pypgcf/utils.py`), often hundreds or thousands of units in one run. If the
run fails partway through — the motivating example: analyzing 5,000
proteomes with `virulence` and running out of disk halfway — there is
currently no way to continue. The user has to diagnose the problem, then
rerun the entire module from scratch, including every unit that already
completed successfully.

## Scope

In scope: `virulence.py`, `amr.py`, `cazy.py`, `smbgc.py`, `eggnog.py`
(`EggNOGRunner`), `orthologues.py`. These all share the same shape: one
external-tool command per input file, with a dedicated (if not always
single-file) output per unit.

Explicitly out of scope for this spec:
- `download_genomes.py` and `phylogenomic.py` — same general problem, but
  deferred to a later spec.
- Workflow-level stage resume (e.g. `pypgcf workflow` skipping an entire
  already-finished module on restart, not just skipping within one module).
  This spec only covers skipping already-completed *units* within a single
  module invocation.

## Why a log file, not output-file-existence checks

Checking "does the expected output already exist" seems simpler, but two of
the six modules (`cazy.py` via `run_dbcan`, `smbgc.py` via `antismash`)
write a whole output *directory* per unit, not one file — and even for the
single-file cases (amr, virulence, eggnog, orthologues), a file left behind
by a crash mid-write looks the same, on a cursory existence check, as a
genuinely complete one. A separate append-only progress log, written to
only after a unit's command exits 0, avoids both problems: it's agnostic to
what shape a unit's output takes, and it can't record a crash as a success.

## New module: `pypgcf/progress.py`

```python
class ProgressTracker:
    def __init__(self, log_file: Path, *, force: bool = False):
        """If force, truncate/remove log_file before reading it, so every
        unit is treated as not-done. Otherwise read existing entries into
        self._completed (a set of unit ids)."""

    def is_done(self, unit_id: str) -> bool:
        """True if unit_id is already recorded as completed. Only ever
        called from the main process, before dispatch starts — no
        concurrent access, no locking needed here."""

    def mark_done(self, unit_id: str) -> None:
        """Append unit_id to log_file. Called from worker processes after
        a unit's command exits 0 — see locking, below."""
```

```python
def run_tracked_command(
    unit_id_and_cmd: tuple[str, Union[str, list]],
    tracker: ProgressTracker,
    debug: bool = False,
) -> int:
    """Run cmd via execute_command; on exit code 0, call
    tracker.mark_done(unit_id). Module-level (not a closure) so it can be
    passed into multiprocess_dispatch via functools.partial."""
```

### Concurrency / locking

`multiprocess_dispatch` already forces a `fork` mp_context (explicitly on
macOS; `fork` is Linux's default), so each worker is a forked copy of the
process that called it — including its own reference to the same
`ProgressTracker` instance, but no shared memory and no shared file handle.
`mark_done` must therefore be safe to call concurrently from independent
processes writing to the same file:

- Each call opens `log_file` in append mode, takes an exclusive
  `fcntl.flock`, writes one line, flushes, and closes — no long-lived open
  handle is kept across calls or across processes.
- `fcntl` is POSIX-only, which is fine since pyPGCF's `pixi.toml` already
  only targets `linux-64`/`osx-64`/`osx-arm64`.
- This is a deliberate choice over relying solely on `O_APPEND`'s
  usual single-short-write atomicity: that guarantee is weaker on network
  filesystems (NFS), which is a realistic environment for a run this size
  on shared/HPC storage.
- `is_done` reads are not locked: they only happen once, in the main
  process, before any worker has started writing.
- Duplicate lines (e.g. if a log somehow gets a unit recorded twice) are
  harmless — `is_done`'s backing set is built by deduplicating on read.

## Per-module integration

Each module's run method currently looks like:

```python
cmds = [create_cmd(f) for f in self.fasta_files]
_ = multiprocess_dispatch("system", cmds, self.concurrent_jobs, ...)
```

This becomes (shape, not exact code):

```python
tracker = ProgressTracker(self.out_dir / ".pypgcf_progress.log", force=self.force)
units = [f for f in self.fasta_files if not tracker.is_done(f.stem)]
skipped = len(self.fasta_files) - len(units)
if skipped:
    logging.info(f"{skipped}/{len(self.fasta_files)} already completed, skipping")
jobs = [(f.stem, create_cmd(f)) for f in units]
_ = multiprocess_dispatch(
    partial(run_tracked_command, tracker=tracker, debug=self.debug),
    jobs,
    self.concurrent_jobs,
    ...,
)
```

### Per-module unit identity

| Module | Unit | Notes |
|---|---|---|
| `amr.py` | one fasta file (stem) | single output file |
| `virulence.py` | one fasta file (stem) | single output file |
| `cazy.py` | one fasta file (stem) | output is a whole directory — fine, the log doesn't care |
| `smbgc.py` | one genome file (stem) | output is a whole directory |
| `eggnog.py` (`EggNOGRunner`) | one fasta file (stem) | see caveat below |
| `orthologues.py` | one forward-or-reverse BLAST command, per (reference, fasta file) pair | unit id e.g. `f"{ref}:{fasta_stem}:fwd"` / `:rev`; forward and reverse are tracked independently |

### `eggnog.py` caveat: `clean_unessecary_output` must only touch units run this invocation

`EggNOGRunner.clean_unessecary_output()` currently loops over **all**
`self.fasta_files` unconditionally, renaming each one's
`<stem>_eggNOG_results.csv.emapper.annotations` to its final name. If a
fasta file was skipped this invocation because it was already completed
(and already renamed) in a prior run, re-running this rename on it raises
`FileNotFoundError` — the intermediate file is long gone. Fix: scope
`clean_unessecary_output` to only the fasta files that were actually
(re)run in *this* invocation (the `units` list above), not the full
`self.fasta_files`.

(Separately noticed, not part of this spec: `clean_unessecary_output` also
builds a `files_to_remove` list for the `.hits`/`.seed_orthologs` sibling
files but never actually unlinks them. Pre-existing, unrelated latent bug —
not fixing it here.)

## CLI changes

Add a `--force` flag to each of the six subparsers (`amr`, `virulence`,
`cazy`, `smbgc`, `eggnog`, `orthologues`) in `cli.py`, following the
existing per-subparser style (each subparser already adds its own
`--debug` individually — `--force` follows the same pattern, no new shared
-arg abstraction). Default `False`. When `True`, the `ProgressTracker` is
constructed with `force=True`, which truncates/removes the existing log
before building the unit list — so every unit is (re)run, and the log
starts fresh (avoiding unbounded duplicate-line growth across repeated
forced reruns).

## Testing

- `test_progress.py` (new): `ProgressTracker.is_done`/`mark_done`/`force`
  behavior; `run_tracked_command` appends only on exit code 0; concurrent
  `mark_done` calls from multiple processes don't corrupt the log (a
  `ProcessPoolExecutor`-based test hammering `mark_done` from several
  workers and asserting every line survives intact).
- Update `test_amr.py`, `test_virulence.py`, `test_cazy.py`,
  `test_smbgc.py`, `test_eggnog.py`, `test_orthologues.py` to cover: units
  already in the log are skipped (not dispatched), `--force` dispatches
  everything regardless of the log's contents, and (eggnog only)
  `clean_unessecary_output` is only called with newly-run files.
- Update `test_cli_subparsers.py`/`test_cli_runs.py` for the new `--force`
  flag on each of the six subcommands.

## Out of scope (explicitly)

- `download_genomes.py`, `phylogenomic.py` — candidates for a follow-up
  spec using the same `progress.py` mechanism.
- Workflow-level (multi-module) stage resume.
- Any attempt to validate that a previously-completed unit's output is
  still *correct* for the current run's parameters (e.g. a different
  `--evalue` on a resumed run) — the log is trusted as-is. If a user changes
  parameters between a failed run and a resumed one, already-completed
  units are not re-run even though their results reflect the old
  parameters. Not addressed here; worth a user-facing warning in a future
  pass if it proves confusing in practice.
