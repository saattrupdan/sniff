# PTR-MS analysis

Sniff is an open-source replacement for proprietary PTR-MS Viewer. It ships a Python
CLI that reads IONICON IoniTOF `.h5` files, detects and quantifies ion peaks, proposes
time segments, and serves a browser-based expert review. This file contains durable
agent-facing contracts; `README.md` and the CLI remain the user documentation.

## Stack

- Python 3.9 or newer, packaged with setuptools through `pyproject.toml`.
- Runtime dependencies: NumPy and h5py.
- The `sniff` console entry point resolves to `sniff.analyze:main`.
- The review UI is generated and served by Python; there is no separate frontend build.

## Layout

| Path | Purpose |
| --- | --- |
| `src/sniff/analyze.py` | CLI parsing, command handlers, CSV output, and orchestration. |
| `src/sniff/ptrms.py` | HDF5 loading, extraction, segmentation, and quantification. |
| `src/sniff/formula_id.py` | Formula enumeration and candidate scoring. |
| `src/sniff/viz.py` | Self-contained browser review UI and localhost server. |
| `src/sniff/gen_rate_constants.py` | Rebuilds the bundled rate-constant JSON. |
| `src/sniff/reference/` | Scientific references and package data shipped with the CLI. |
| `packaging/` | PyInstaller spec and frozen-entry point; the `package` workflow builds a folder bundle per OS. |

## Agent-facing contracts

- **CLI boundary:** use `sniff <command> --help` and the CLI's JSON output for
  inspection and analysis. Do not perform ad-hoc HDF5 analysis in an agent or script;
  keep scientific decisions in the package's CLI and library.
- **Review ownership:** the browser review belongs to the human reviewer. Prepare the
  config before opening it, then hand over the URL or app. Do not use browser
  automation to inspect or edit the review or click Done.
- **Peak scope:** decide comprehensive versus targeted output before selecting peaks.
  Comprehensive output retains every credible real channel, including fragments,
  isotopes, reagent ions and water clusters; targeted output is only for an explicit
  named panel. Never silently reduce a general export to familiar VOCs.
- **Blank files:** if the CLI reports no significant reagent ion and no peaks, treat
  the run as blank/no-beam/aborted. Report it and stop; never lower thresholds or
  invent analytes from noise.
- **Background diagnostic:** after ranges are curated, inspect the `background`
  diagnostic. Channels stronger in backgrounds (S/B below 1) or with an upward
  background trend are background or contamination: relabel or drop them rather than
  shipping them as analytes.
- **Range labels:** preserve chronological, deterministic labels: `sample_01`,
  `sample_02`, and `background_01`, `background_02`, numbered independently. Do not
  ask users to name automatically detected plateaus.
- **Concentration communication:** distinguish file-derived scale from a project or
  standards calibration. State when K is uncalibrated and treat humidity-sensitive
  compounds as indicative unless their calibration supports more. Never imply that a
  plausible number is accurate without evidence.
- **Identification limits:** m/z and formula candidates are proposals, not proof of
  chemical identity. Preserve honest unknowns, report ambiguity and overlap, and do
  not present a library match or candidate score as a calibrated probability.
- **Mass-axis authority:** a valid `CALdata/Mapping` with three or more references is
  authoritative. Do not apply a second water-cluster/iodobenzene affine correction;
  those measured peaks are temporal-stability diagnostics on this path. Two-row Mapping
  and Spectrum fallback files retain the mandatory two-reference affine correction.
- **Tolerance layers:** preserve the distinction between broad 200 ppm formula/name
  proposals for expert investigation and the run-validated 5–10 ppm assignment radius.
  Wider proposals stay visible but must never become automatic identities.
- **Fragmentation evidence:** use condition-matched PTR Library pathways and measured
  parent/fragment co-variation only to rerank existing candidates or suggest fragment
  roles. It must not create formulas, override mass eligibility, prove an isomer, or be
  described as MS/MS evidence.

## Repository workflow

- Work directly on `main` in this repository; do not create a feature branch.
- After every completed change, rebuild and install `/Applications/Sniff.app` so Dan
  can test it by opening Sniff. The orchestrator may perform this step outside an
  isolated builder worktree.

## Running it

For development, use the checkout's project environment:
```bash
uv sync
uv run sniff --help
uv run sniff rates water
```

Use real IoniTOF data only when exercising file-dependent commands. `.h5` files can be
about 1 GB, and a full analysis commonly takes about a minute.

## Validation

Run the automated tests, lint check, build, and CLI smoke checks after a change. The
browser regression is especially relevant when changing identification display; it
requires the `agent-browser` CLI (`npm i -g agent-browser` and `agent-browser install`)
in addition to the package's normal Python dependencies:

```bash
uv run pytest
uv run ruff check --select F,I src tests scripts
uv run python scripts/smoke_viz.py
uv run sniff --help
uv run sniff inspect --help
uv run sniff peaks --help
uv run sniff segments --help
uv run sniff analyze --help
uv run sniff viz --help
uv run sniff app --help
uv run sniff calibrate --help
uv run sniff compare --help
uv run sniff rates h2o    # the bundled library has no water entry; h2o returns matches
uv build
```

For scientific or HDF5-processing changes, also run the affected command on a suitable
local fixture and inspect its JSON diagnostics or CSV output. Do not commit measurement
files, generated review HTML, configs, or result CSVs.

Packaging is checked separately because it is slow and pulls its own toolchain:
Remove stale `dist/` and `build/` output, then run
`uv run --with pyinstaller pyinstaller --noconfirm packaging/sniff-app.spec` and
`uv run python scripts/smoke_frozen.py "dist/Sniff.app/Contents/MacOS/sniff"` on
macOS, or `dist/sniff/sniff.exe` on Windows, which starts the frozen bundle and asserts
it serves the review page. The executable is at the bundle root; PyInstaller
uses the app's `Contents/Resources/` on macOS and `_internal/` on Windows for its
package data and libraries. Run it when `packaging/`, dependencies, or the app server
change; the Windows half of it can only be verified on a Windows runner.

## Conventions

- Keep compatibility with Python 3.9; do not introduce Python 3.10+ syntax without first
  raising `requires-python` deliberately.
- Write documentation and new prose in British English, wrapped at 88 characters.
- Preserve JSON on stdout for discovery commands and send progress logs to stderr.
- Keep the CLI deterministic. Chemistry assignment and segment curation remain explicit
  agent decisions; do not add a one-shot automatic workflow.
- Update `README.md` when flags, output fields, workflow, or scientific interpretation
  change. The `analyze` object is resolved with CLI override > curated
  config > legacy default; keep the Methods provenance and authoritative Done rerun
  wording aligned with that implementation.
- Use Conventional Commits, following the parent dotfiles repository.

## Gotchas

- `src/sniff/` is the installable package. Keep package-internal imports
  relative and use `importlib.resources` for bundled data; do not reintroduce flat
  top-level modules.
- Run checkout commands with `uv run sniff`; package modules are imported as
  `sniff.*`.
- `src/sniff/reference/rate_constants.json` is generated from
  `src/sniff/reference/ptrlibrary.csv` by
  `uv run python -m sniff.gen_rate_constants`. Change the source or generator,
  regenerate the JSON, and review both files together rather than hand-editing entries.
- `src/sniff/reference/compound_catalogue.sqlite3` is exported from the complete
  WebBook crawl state, not hand-edited. Run `crawl_webbook_species.py reparse` and
  require `ready_to_classify: true`, then run `build_compound_catalogue.py export-full`.
  Never commit the full crawl database or cached HTML pages.
- Reference Markdown, CSV, and JSON files under `src/sniff/reference/` are
  package data. Keep `pyproject.toml` in sync when adding a new bundled file type.
- `viz` reviews an already curated config; it must not silently perform peak or segment
  detection. The delivered CSV is always produced by the analysis path.
- Preserve 1-based, inclusive cycle ranges and deterministic chronological labels:
  `sample_01`, `sample_02`, and `background_01`, `background_02`, numbered separately.
- Do not replace missing calibration or transmission data with plausible-looking values.
  Surface degraded accuracy through the existing diagnostics and NaN behaviour.
- Do not weaken noise, overlap, apex, humidity, blank-file, or sample/background checks
  merely to produce more populated output. These are scientific safeguards.
- Sniff replaces proprietary PTR-MS Viewer. Documentation must not instruct users to
  generate, validate, or repair results in that tool; an existing reference CSV may only
  be used for calibration or comparison.
