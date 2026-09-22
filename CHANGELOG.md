# Changelog

All notable changes to `sniff` are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- The bundled offline catalogue now includes a bounded, response-hashed PubChem PUG
  REST snapshot for formula-supported name and structure proposals. PubChem evidence is
  source-separated, locally packaged, never queried at runtime, and cannot create or
  promote a formula assignment.
- New schema-3 reviews use joint `formula-envelope-v2` fitting across connected assigned
  parent and M+1/M+2 channels. Rank, condition, weighted residual and uncertainty gates
  withhold inseparable or inconsistent isotope corrections; auxiliary channels remain
  excluded from analyte rows and legacy reviews keep `formula-v1` arithmetic.
- Peak discovery now reports versioned candidate coverage and separate alternative-ion
  proposals for condition-compatible hydration, water loss, charge transfer, hydride
  abstraction and evidence-backed multiply charged envelopes. Isotope and known-fragment
  channels inherit parent compound proposals without becoming automatic neutral-analyte
  assignments. Versioned selected roles persist through review, while canonical
  unresolved-peak components and one-per-component Raw weighting prevent overlap from
  inflating candidate coverage.
- New analyses learn an empirical peak shape from clean isolated channels and use
  bounded, non-negative run/interval fits for overlapping groups. Rank, conditioning,
  correlation, residual and fallback diagnostics prevent unresolved components from
  receiving plausible-looking independent concentrations.
- Accepted formulas now derive exact natural M+1/M+2 auxiliary channels, show
  expected/observed isotope evidence, and support guarded isotope-spillover correction
  without adding duplicate analyte rows. Monoisotopic-abundance scaling remains withheld
  unless the calibration basis explicitly supports it.
- Hand-drawn review peaks now snap to a measured apex and request formula candidates
  from the local server. **Use table on another file** conservatively transfers
  identities to an explicitly selected run while retaining every credible new detection.
- New app review panels now give every peak with an available unique candidate an
  editable best guess: a preferred or canonical PTR Library name when available,
  otherwise its best candidate formula. Experiment-specific compounds influence initial
  ranking and supply the preferred provisional isomer name; global matching prevents an
  automatic compound/formula family from being assigned twice and uses the next-best
  available candidate instead.
- Sniff now ships a generated, read-only SQLite catalogue derived from the complete NIST
  Chemistry WebBook sitemap. Its filtered 20,106 formula families and 93,946 ordinary
  neutral species provide offline exact-mass, formula, name, CAS and InChI lookup while
  keeping external names visibly separate from PTR-specific evidence. The maintainer
  crawler caches every canonical page for offline classification, resumes from per-page
  checkpoints, refreshes robots policy daily, observes NIST's host-wide five-second
  delay and stores no HTML or spectra in the package.
- Every detected peak now has either a strict formula candidate or an explicit
  scientific interpretation. Known reagent ions, evidence-backed possible isotope
  channels and detector artefacts are excluded from automatic analyte assignment;
  unmatched channels are labelled as unresolved ions rather than being given fabricated
  compounds.
- Peak discovery now detects high-prominence detector echoes that evade local noise
  rules. Classification requires a delay recurring behind several taller peaks in
  calibrated time-bin space plus parent/satellite level and cycle-change co-variation;
  supporting-only evidence remains visible and no formula-bearing peak is suppressed.
  Curated reviews also expose complete-range sample/background evidence as a separate
  non-compound role without changing chemical assignments.
- PTR Library source rows now retain product-ion masses, branching percentages, E/N,
  instrument and citation provenance. Under compatible measured reaction conditions,
  parent/fragment level and change co-variation can rerank existing formula proposals
  and link possible fragment channels without creating candidates or overriding the
  exact-mass assignment gate.

### Changed

- Formula enumeration and natural-isotope evidence now include silicon and iodine,
  covering common siloxane backgrounds and iodinated reference ions without widening
  either mass tolerance. Legacy Gaussian overlap projection now uses a guarded,
  scale-relative regularisation and withholds non-finite components instead of emitting
  overflow-driven traces.
- Valid three-or-more-point `CALdata/Mapping` fits are now authoritative mass axes.
  Sniff no longer applies a second water-cluster/iodobenzene translation and scale that
  could move those references hundreds of ppm after a single-digit-ppm Mapping fit.
  Version-1 review coordinates migrate through their recorded historical transform so
  they retain the same physical timebins.
- Formula discovery now separates broad 200 ppm expert-review proposals from a
  run-validated 5–10 ppm assignment radius derived from independent Mapping residuals.
  Wider NIST/catalogue matches remain visible with their ppm errors but are never
  assigned automatically; poor calibration leaves every result as a proposal rather than
  widening the scientific acceptance boundary.
- Fresh app reviews now fill remaining identity blanks after trace extraction refines
  the measured peak apex. Reagent markers first use the strict 12 mDa window; three
  consistent marker anchors may then recover another known reagent channel with the
  same small centroid displacement without moving the mass axis or relaxing formula
  tolerances. The former vague m/z 30.994 O₂⁺/NO⁺ region label is now the specific NO⁺
  (¹⁵N) isotope, preventing nearby formaldehyde and reagent-region satellites from
  receiving misleading reagent labels.
- New app-generated configs use analysis schema 3 with `empirical-v1` fitting and joint
  `formula-envelope-v2` isotope handling. Schema-2 configs retain `formula-v1`, and
  legacy configs retain their previous Gaussian and isotope-off arithmetic.

## [0.8.0] - 2026-09-12

### Added

- Sniff now checks GitHub Releases on startup and prompts when a compatible update is
  available. **Install update** downloads, verifies and opens the platform installer in
  one click; **Install later** defers the prompt until the next launch.
- Linux x86_64 releases now include an Ubuntu/Debian `.deb` installer and a portable
  `.tar.gz` bundle, both built and smoke-tested on Ubuntu 22.04. Linux deliberately
  opens Sniff's local interface in the default browser to avoid distro-specific GUI
  toolkit dependencies; the `.deb` adds the application-menu entry, icon, native file
  chooser dependency and `/usr/bin/sniff` command.

### Fixed

- Browser-based file selection now invokes KDE's `kdialog` with its supported open-file
  option, so Linux systems without `zenity` can still choose an H5 run.

## [0.7.1] - 2026-09-12

### Fixed

- The compound autocomplete now scrolls within a shorter suggestion list instead of
  making the entire pre-analysis modal scroll.

## [0.7.0] - 2026-09-12

### Added

- After an H5 file is selected, the app now offers an optional compound-of-interest
  modal. It validates typed or pasted names against the bundled PTR Library, supports
  autocomplete and comma/newline-separated input, and records recognised compounds as a
  transparent contextual prior for candidate ranking. The modal can be skipped, and the
  prior never forces an identification or suppresses other detected peaks.

### Fixed

- The frozen-app smoke now isolates launch modes. On hosted macOS, where every second
  invocation of a windowed executable can stall before Python starts, CI tests the
  installed PKG once and verifies its review and start screens in that process. Windows
  CI and local macOS checks retain the independent no-argument and window probes.

### Changed

- The app landing page now shows the installed Sniff version discreetly in its footer.
- Reopening an unchanged H5 in the app now reuses the fully validated mass-axis
  calibration stored in its config. A versioned file-identity fingerprint prevents a
  changed, replaced or copied H5 from inheriting stale calibration evidence; uncertain
  and legacy cases still recalibrate. The opening stage now distinguishes loading H5
  data for a saved review from computing data for a new one.
- Removed the browser review checklist. New configs no longer generate checklist data,
  and older top-level or nested checklist fields are discarded while unrelated config
  fields are preserved.

## [0.6.1] - 2026-09-10

### Fixed

- **The Windows desktop app no longer leaves a command prompt open.** The Start Menu
  launches a windowed executable, while a separate `sniff-cli.exe` preserves terminal
  commands and JSON output; console-free diagnostics continue in Sniff's log file.
- **Restarting Sniff now returns to the opening screen.** Selecting a previously
  reviewed H5 file still reuses its same-basename JSON config instead of detecting peaks
  and intervals again.
- **Windows paths no longer fill the review header.** The top bar shows only the H5
  filename for both Windows and POSIX paths.
- **The Intervals card no longer shows automatic gap-joining diagnostics.** The
  underlying interval provenance remains in the config and CLI diagnostics without
  occupying the review interface.

## [0.6.0] - 2026-09-10

### Fixed

- **Review edits now survive immediate closes and app restarts.** Autosave flushes the
  latest config as a review page or desktop window closes, rejects delayed older writes,
  reports failed saves honestly, and the next app launch resumes the active run.
- **The packaged app icon now fits the macOS Cmd+Tab switcher.** Its coloured tile is
  drawn inside the platform safe area instead of filling the entire icon canvas and
  appearing oversized beside other apps.
- **Supplied mass axes now require complete calibration evidence.** Production
  boundaries reject contradictory coefficients, models, anchors, peak-quality fields,
  timebins and persistence records instead of trusting a forged applied diagnostic.
- **Mass-axis persistence ignores unrelated corrupt bins.** Raw-cycle anchor checks now
  use only finite samples in each anchor window, while rejecting anchors without two
  usable cycles and preserving the existing block-persistence threshold.

### Changed

- **The guided tour now opens automatically only for the first analysis.** Sniff stores
  that app-wide onboarding state under `~/.sniff`, so later analyses and app restarts do
  not repeat it; the help button can still replay it manually.
- **The analysis loading sheet gives its stage text more room.** Extra space before the
  progress bar keeps the status from feeling cramped.
- **The Sniff mark is now a focused mass-spectrum peak.** The web mark, start-screen
  artwork, Dock icon and Windows icon share the same clean peak without a nose or
  mascot. The macOS app also applies its bundled icon directly at launch so Cmd+Tab
  cannot retain stale artwork.
- **The native desktop app now opens maximised.** It uses the full usable screen as a
  normal window while retaining its explicit width and height arguments for
  compatibility.
- **Mass-axis calibration is now one validated, cancellable operation per open.** Its
  progress is included in the opening bar, raw-cycle persistence is scanned once, and
  the validated axis is reused by detection, review preparation and extraction.
- **Calibration and review round-trips now fail closed and preserve provenance.** Raw
  cycle persistence is mandatory for both internal anchors, caller-supplied axes must be
  applied internal calibrations, and browser saves retain authored peak/range fields
  while updating only the fields edited in the review.
- **Opening a run now has a little more life and a smoother finish.** The progress sheet
  carries restrained, accessible peak, ion and breath motion without covering the file,
  stage, percentage, ETA or Cancel. Completion settles at 100%, hands off into the
  review with a guarded zoom/fade transition, and respects reduced-motion preferences.
- **Agent-facing workflow contracts now live in `AGENTS.md`.** The retired workflow
  document no longer duplicated the CLI boundary, review ownership, peak scope, blank
  handling, diagnostic checks, range labels, concentration caveats or identification
  limits.

- **Review saves now require the corrected-axis marker.** Legacy configs remain readable
  only long enough to migrate after successful calibration; new saves and exports cannot
  silently reintroduce file-axis masses.
- **Mass calibration now fails closed on two required internal standards.** The
  operational water calibrant (37.033) and protonated iodobenzene (204.951) must both
  pass prominence, S/N, ambiguity, proximity and raw-cycle persistence checks. Missing
  or implausible anchors produce structured errors rather than silently using the HDF5
  axis; corrected forward/inverse mapping remains shared by detection, extraction,
  identification, analysis and browser interaction. Unmarked old configs are migrated
  once from file-axis masses to corrected-axis masses and marked with their domain and
  version.
- **The product is now Sniff.** The installable package and source directory are
  `sniff`, and the command is `sniff`. Desktop state lives under `~/.sniff`; existing
  `~/.ptr-ms/recent.json` and beside-file `.ptr.json` data remain readable and are never
  deleted. Installer payloads and generated artefacts use the Sniff names while the
  stable macOS bundle identifier, MSI upgrade code and component identities remain
  unchanged.
- **The Sniff opening screen is now browse-only.** It uses the native file dialog rather
  than path entry, while retaining Browse and the open-review affordance. Recent-file
  history remains in the API and storage layer rather than being displayed in the UI; an
  unavailable native dialog is reported without exposing a path-field fallback.
- **Frozen one-dir bundles keep the launcher separate from package contents.** Windows
  stores the Python payload under `_internal/`, while the macOS app uses its standard
  `Contents/Resources/` area; both retain the visible `sniff` executable path.

## [0.5.0] - 2026-09-08

### Added

- **Opening a file runs behind a progress you can watch and cancel.** The start screen
  now puts a half-minute open into a full-screen sheet — the file's name, the stage, a
  determinate bar, a rough ETA and **Cancel** — instead of a line of text and a spinner
  while the reviewer wonders whether the app is wedged. The bar is driven by the work:
  reading the run is 28.6 s of the ~33 s an open takes on the 2 GB / 20,725-cycle
  fixture (14.6 s reading every cycle once and 13.9 s re-reading the intervals to
  re-centre peaks on them), so it gets 89 % of the bar and reports cycles read rather
  than a smoothed guess, and the quick phases in front of it are what the other 11 % is
  made of. Cancelling stops the analysis at the next block, closes the file and returns
  you to a usable start screen: the session goes back to empty with no error string and
  nothing half-written, and the same file can be opened again immediately. When an open
  finishes, the page goes to the review by itself — so **Open the review** is no longer
  a button you wait beside, and Back does not land on a sheet for a file that is already
  open. `GET /api/state` gained `progress` and `cancellable`, and `POST /cancel` stops
  an open in flight (and is a no-op when there is none). Analysis is untouched: with no
  callbacks attached, `extract_traces` produces numerically identical traces.
- **The app can open in its own desktop window.** `sniff app --window` runs the review
  in a single window with no address bar, `pywebview` being an extra
  (`uv sync --extra desktop`) rather than a dependency, and a packaged bundle uses the
  window by default because a double-clicked app has no terminal to read an address out
  of. Closing the window stops the server and **Browse this computer…** uses the
  window's own dialog when there is one; a machine without the extra, or without a
  display, says so once and serves a browser tab exactly as before. The `package`
  workflow installs `.[desktop]` and the PyInstaller spec bundles `webview` when it is
  present, so the `.pkg` and `.msi` ship the window; the frozen smoke runs `--window` on
  a runner with no display to prove the fallback rather than assume the window.
- **The start screen was rebuilt, and gained a real file dialog.** Recent runs are
  listed with their size, when you last opened them and whether a config exists yet;
  **Browse this computer…** opens the desktop's own file picker, because the server runs
  where the files are, which is the one thing a web page normally cannot do. The page is
  a single-column layout with proper type, focus rings, a dark scheme and a progress bar
  while a run loads.
- **Double-clicking the bundled app opens it.** Finder starts the bundle with no
  arguments, and the plain `sniff` command line answers that with usage text and exit
  code 2 — invisibly, in a windowed bundle. A runtime hook turns a bare launch inside a
  bundle into `sniff app`, and when no browser opens the address is written to
  `~/.sniff/log.txt` rather than vanishing.
- **Packaging: `packaging/sniff-app.spec`, `packaging/make_msi.py`,
  `scripts/smoke_frozen.py`, and a `package` workflow.** PyInstaller builds a bundle a
  reviewer can run with no Python installed (one-dir by choice — one-file unpacks into
  `%TEMP%` on every start and is what antivirus tools object to), wrapped as a `.app`
  that a `.pkg` installs on macOS and an `.msi` on Windows. The MSI's component list is
  generated from the built folder with path-derived ids and GUIDs, so an upgrade
  replaces and removes exactly the right files, and it targets the MIT-licensed WiX v3
  rather than v6+, which will not run until the Open Source Maintenance Fee EULA is
  accepted. PyInstaller cannot cross-compile, so the workflow builds on native macOS and
  Windows runners and uploads one installer per platform, with a `v*` tag publishing
  them as a GitHub Release. The smoke script starts the finished bundle against a tiny
  synthetic file and fails unless it serves the review page, since `sniff --help` would
  pass on a bundle that can do nothing else.
- **"Stop the app" on the start screen, and a log file for bundle runs.** A
  double-clicked app has no terminal to press Ctrl-C in, so the page can shut the server
  down itself (`POST /shutdown`, and `SIGTERM` now quits the same way); with no console,
  its URL and errors are appended to `~/.sniff/log.txt`.
- **App mode: `sniff app`.** A persistent local review app for people who want the tool
  rather than the chat. It opens on a start screen of recent files, takes a file from
  there, and stays up between files. Each file's config lives beside it under the same
  stem (`sniff.h5` → `sniff.json`, with an existing `<stem>-analysis-config.json`
  honoured), so reopening a reviewed file returns exactly what was saved; a file that
  has never been reviewed gets the deterministic peak and interval pipeline written to
  that path, with a checklist that says plainly that nothing has been curated yet and
  which calls are still a human's. The primary button is **Export** instead of Done: it
  writes `<stem>.csv` beside the file and leaves the app open for more work. With
  `--agent URL` (or `SNIFF_AGENT_URL`) a newly generated config is offered to an agent
  for curation first, and the deterministic config is kept — visibly — whenever that
  endpoint is missing, slow or unhelpful. The server binds to 127.0.0.1 and the only
  outbound request is to the endpoint the user named.
- `analyze.auto_peaks` / `analyze.auto_ranges`: the detection pipeline is now callable
  without argparse, which is what lets the app build a config on the user's behalf.
- The Peaks sidebar tick is now sample-specific and follows the **average over** choice:
  with one sample selected it is that sample's own tick, with the whole run selected it
  is the aggregate — ticked for every sample interval, empty for none, a dash for some.
  Clicking the aggregate box only ever flicks it between ticked and empty, so a dash
  becomes a full tick rather than a dead end. The Details view shows one box per sample
  interval in every row, and the heading names the scope. A compound in only some
  samples records `samples` in the config; a compound in every sample needs no new
  field, and the summary output is unchanged either way.

### Changed

- **The start screen no longer has a Stop the app button.** Closing the window is the
  way out of a windowed app, and the same quit remains available in a browser tab —
  where there is no window to close — as one quiet footer link, with `POST /shutdown`
  and Ctrl-C unchanged.

- **The desktop app is called Sniff.** `PTR-MS Review` described an instrument and an
  activity, which is right for a command line and bland for an icon in a Dock. The name
  and the mark now live in `sniff/brand.py` instead of being spelled four ways, and the
  mark itself is drawn by `packaging/make_icons.py` — a teal tile, the mass-spectrum
  trace, one warm nose above the tallest peak — which also rasterises it into the
  `.icns` and `.ico` the installers attach, so the Dock icon is no longer the generic
  page Finder invented. On macOS the bundle identifier changed with the name, so the
  `.pkg` installs alongside an old copy rather than over it; `packaging/README.md` gives
  the one command that removes it. The Windows upgrade code did not change, so that one
  does upgrade in place. The command remains `sniff` and the distribution remains
  `sniff`.

- **The Peaks sidebar says less.** Its heading is just **Peaks** — no dot and no
  `· sample_03` suffix; the compound name and the row of numbered boxes above the list
  are gone, as are the numbered boxes in the Details rows and the **sample/background**
  switch. What remains is the interval selector, which already decided what the tick
  boxes and the mass spectrum speak about, so the **average over** dropdown has left the
  Mass spectrum header rather than leaving a second control saying the same thing.
  Nothing became unreachable: an interval's class is set in the **Intervals** card, and
  which samples a compound belongs to is read by selecting those intervals. Changing the
  class there also stops losing the per-sample ticks on the way back — the row handler
  moved the class before asking `setSampleClass` what it was leaving, so reclassifying
  an interval to background and back used to leave the compound out of it.
- **Panel close buttons sit on the right.** The ✕ in the Configuration, Method and
  Checklist panels was pushed by a spacer that only had a rule inside a `.card`, and
  these panels are not cards.

- **The macOS artifact is a `.pkg`, and `packaging/README.md` now explains both
  installers.** A disk image asked the reviewer to drag a bundle into Applications; a
  product archive installs it, records what it wrote, and installs from a command line,
  which is how CI now proves the artifact rather than the build folder.
  `packaging/make_pkg.py` writes the `productbuild` distribution — title, version, a
  macOS 11.0 floor and the architecture, with no paths and no timestamps, so the same
  checkout gives the same XML twice — and the `package` workflow pairs it with
  `pkgbuild --component … --install-location /Applications`, then installs the result
  with `sudo installer -pkg` and smokes `/Applications/Sniff.app/Contents/MacOS/sniff`,
  exactly as the Windows job smokes both `dist/sniff` and `C:\Program Files\Sniff`. The
  guide has both command pairs, the reason WiX v3.14 is pinned (v6 and later are gated
  behind the Open Source Maintenance Fee), how to inspect a `.pkg` (`lsbom`, `pkgutil`)
  and remove one (there is no uninstaller), what a double-clicked bundle gets, and what
  is still missing: no Developer ID signing or notarization, and no Authenticode, so an
  unsigned `.pkg` still trips Gatekeeper on first run and SmartScreen warns on Windows.
- **Adjacent plateaus are now joined on evidence, not on a cycle count.** A gap between
  two same-class plateaus merges only when it covers under ~60 s of acquisition (never
  fewer than 30 cycles) _and_ never left the phase its neighbours are in: its highest
  cycle stays below the higher neighbour's level times 2, read against the run's own
  background, and for a sample its lowest cycle also stays above the lower neighbour's
  level divided by 2. A background has no lower test, since a background cannot fall out
  of itself — a dropout toward zero is still the same blank, and splitting it would cost
  the longer reference interval a blank exists to provide. A wobble inside one sample
  therefore merges at 1 s/cycle and at 5 s/cycle alike, while a sample that fell back to
  the background — or a gap that spiked out of the phase — stays a break however short
  it is. An opposite-class plateau between two segments is still a hard boundary, and
  each gap is judged against the plateau it abuts rather than the running average of a
  partly merged interval, so the verdict does not depend on which plateau came first. On
  a 20,725-cycle IoniTOF run this gave 12 sample and 11 background intervals where a
  reviewer curated 12 and 10 by hand, and where the length rule at 30 cycles gave 17 and
  7 — it joined two samples whenever the gap between them happened to be short.
  `--merge-high-gap N` survives as a cap override and `0` still means never join high
  plateaus; `sniff segments --merge-high-gap` and `sniff analyze --merge-high-gap` now
  default to the automatic test instead of off.
- **Every merge explains itself.** `merged_gaps` now carries, per gap, its length, the
  gap's minimum and maximum level and a reason (`level held`, `fell to baseline`,
  `adjacent`, or `length only` on the legacy path), and the review app says the same in
  one line on the Intervals card — `joined 2 wobbles, level held (≤ 28 cycles)` — via
  the config's `merge_note`. The reviewer in `sniff app` never sees a command line, so a
  silent join of their intervals would have been unfalsifiable. A merged interval's
  level is the mean of its plateaus weighted by plateau cycles, so the cycles between
  them cannot drag the reported level toward the baseline.
- The interval in scope and its sample/background class now sit in the Peaks sidebar,
  one click away on either tab instead of only in the Intervals card. Switching class
  renames the interval, because the analysis reads the class from the interval name;
  switching it back restores the name and the recorded per-sample selections.
- The composite VOC curve behind the Signal over time trace is labelled in the plot, and
  its legend entry is a switch kept in the config as `viz.show_disc`.
- The Intervals card shrinks as far as the splitter is dragged: the plot is no longer
  capped at 560 px, which left the card never smaller than half a tall window.
- The Intervals card updates while an interval edge is dragged and keeps its rows in
  chronological order, so the table no longer lags behind the plot.
- The Mass spectrum "Average over" list follows interval renames, recolouring and
  resizes instead of showing stale names, and re-averages the spectrum when the interval
  it points at changes shape.
- The Raw / Corrected / Conc / µg selector now also drives the mass spectrum and the
  sidebar values. In Conc and µg the sidebar figure is the mean of the compound's own
  converted trace over the cycles being shown — the number the CSV reports as `Average`
  for that interval. The shared spectrum axis carries only the conversion every compound
  shares, and says so; the per-compound humidity correction stays per compound.
- A compound name never contradicts its identification: an auto-generated
  `unknown m/z …` label on a peak with an assigned formula is replaced by that formula,
  a hand-drawn peak is named from the library only within 10 mDa of a library mass, and
  a label naming a different formula than the assigned one is flagged.
- Interval labels must be unique, because they key the per-sample selection.

### Fixed

- **Exporting no longer ends on a dead end.** The results dialog used to turn its own
  button into a disabled label reading "Opened ✓ — you can close this tab" — which does
  nothing in an app that has no tab to close — and a failed export left an error card
  with no way out at all. Revealing the CSV now closes the dialog and returns you to the
  review, **Keep reviewing** closes it without revealing, and both routes leave the
  **Export** button ready to run again. The one-shot `sniff viz` flow is unchanged.

- **The installed app opens its own window again.** A packaged bundle was opening a
  browser tab and saying nothing useful about it. `desktop.py` reaches pywebview through
  `importlib.import_module()`, which PyInstaller never follows, so although pywebview
  itself was collected, `bottle` and `proxy_tools` — which `webview/__init__` imports at
  module scope — were not. The first `import webview` failed with
  `No module named 'bottle'`, and that was reported as _"the desktop extra is not
  installed"_ while the extra sat inside the bundle. The spec now bundles pywebview's
  declared dependencies and the GUI toolkit its installed backend imports; the message
  blames the packaging rather than the user; and `/api/state` reports `surface` as
  `window` or `browser`, because a `console=False` bundle logs to no terminal and the
  user had no way to tell. The frozen smoke test reads that field now — it used to grep
  the log for the word "window", which a silent bundle satisfied by saying nothing, and
  which the failure line also contained.

- **An autosave no longer invents compound names.** A peak with neither a label nor a
  formula needs something to draw on the spectrum, so the review page displays a
  mass-derived stand-in (`m17.032`) — and wrote it back into the config as though
  someone had assigned it. A freshly detected file therefore gained 132 names that
  looked curated but were not, and its CSV printed the mass twice. The stand-in is now
  display-only, and the saved config keeps the name empty until someone chooses one.
- **A file opened from the recents list no longer appears twice on the start screen**:
  once in its own panel with **Open the review**, once as a plain row. The recents list
  now says which entry is the open one, and leaves it out.
- **The app no longer reports itself busy the moment after it reports itself ready.**
  Readiness was announced before the in-flight flag was cleared, so a client that acted
  on "ready" could have a close or an export refused.
- **An error you caused stays on screen.** The status poll cleared the message within a
  couple of seconds, so a failed open left no trace of what went wrong.
- ↑/↓ move the peak selection down and up the order the sidebar is showing. They had
  always walked the stored m/z order, so with the list sorted by abundance or label they
  jumped to compounds that were nowhere near the highlighted row.
- The Raw / Corrected / Conc / µg buttons now sit in the same place on both tabs: the
  x-axis selector moved to the right-hand slot that **average over** occupies on the
  Mass spectrum tab, and a long interval list can no longer squeeze the tab buttons
  sideways as a side effect.
- Reclassifying an interval from sample to background now reaches the saved config. It
  used to change only the in-memory colour: the class is read back from the interval
  name at load, so an unrenamed interval came back as a sample. A compound's `samples`
  list also survives the trip — before, an interval classed away and back silently
  dropped compounds that were in only some samples.

## [0.4.0] - 2026-08-24

### Added

- Added a check/uncheck-all toggle to the Peaks sidebar.
- Added alphabetical label ordering alongside m/z and abundance ordering.
- Added a persistent draggable splitter between the plot and context card.

### Changed

- `viz` now waits indefinitely by default; `--timeout` remains available as an explicit
  opt-in limit.
- Improved responsive sizing and alignment of the Peaks sidebar and context cards.
- Removed redundant card guidance and improved the guided-tour button's light-mode
  contrast.

## [0.3.0] - 2026-08-23

### Added

- Added a peaks-sidebar ordering selector for m/z order or descending abundance, where
  abundance is the mean per-cycle integrated Raw signal. The compact list now shows only
  the active sort field; details shows both m/z and abundance. Sidebar values follow the
  Mass spectrum tab's selected average-over interval.

### Changed

- Absolute-time displays now apply the file's `UTC_Offset` so plot, crosshair, and
  interval times use the lab PC's local wall-clock time.

## [0.2.1] - 2026-08-22

### Fixed

- Restored the detailed scientific sections in the browser review Methods panel while
  retaining its live effective-settings and provenance summary.

## [0.2.0] - 2026-08-22

### Added

- Added the configurable `viz` x-axis unit selector. It accepts exactly `cycle`,
  `relative`, and `absolute`, with CLI/config precedence and cycle-based range and CSV
  persistence. Relative time uses valid `PCTime` values before a duration-based
  fallback; absolute time requires valid timestamps within ISO UTC years 0000–9999.

## [0.1.2] - 2026-08-22

### Fixed

- Fixed browser-review preparation failing when `SPECdata/AverageSpec` contains
  non-finite bins; `viz` now treats those rare corrupt bins as zero, matching peak
  detection and trace extraction.

## [0.1.1] - 2026-08-21

### Fixed

- Fixed the former too-many-values-to-unpack error when `CALdata/Mapping` contains more
  than two anchors; three or more `(m/z, timebin)` anchors are now accepted only after a
  well-conditioned least-squares fit and a finite reconstructed-mass residual check of
  at most 100 ppm. This is a deliberately generous corruption/model-consistency ceiling,
  not an accuracy claim. Mapping still falls back to usable per-cycle `CALdata/Spectrum`
  coefficients when absent or unusable.

## [0.1.0] - 2026-08-21

### Added

- Initial standalone Python package for processing IONICON IoniTOF PTR-MS and PTR-TOF
  `.h5` files.
- `sniff` command-line interface with commands for inspection, peak detection,
  segmentation, analysis, visual review, calibration, comparison, and rate constants.
- Browser-based review workflow for curating detected peaks and time segments.
- Formula identification, transmission correction, concentration conversion, and bundled
  PTR-MS reference data.
- Package metadata, resource loading, automated tests, and CLI smoke checks for
  installation in isolated environments.
