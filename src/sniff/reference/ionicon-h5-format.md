# IONICON IoniTOF HDF5 format & the PTR-MS Viewer algorithm

Reverse-engineered from an IoniTOF PTR-MS breath run and its PTR-MS Viewer CSV
export. This documents the file layout, the calibration data, and the exact
Raw → Corrected → Concentration chain the pipeline reproduces.

## File layout

Root attributes hold instrument/run metadata: `Single Spec Duration (ms)`
(cycle length), `Timebin width (ps)`, `Pulsing Period (ns)`, `FileCreatedTime*`,
`UTC_Offset` (lab-PC offset in seconds), `InstrumentType = IoniTof`, drift settings,
etc.

| Path | Shape | Meaning |
|---|---|---|
| `SPECdata/Intensities` | (n_cyc, n_bins) | **Raw mass spectra**, one row per cycle, columns = TOF timebins (cps). The bulk of the file (~1 GB, gzip, chunked one-row-per-chunk). |
| `SPECdata/AverageSpec` | (n_bins,) | Run-average spectrum — used for peak detection / apex finding. |
| `SPECdata/Times` | (n_cyc, 4) | col 0 = 1-based cycle index; col 2 = acquisition time (IONICON epoch). |
| `SPECdata/PCTime` | (n_cyc, 1) | PC Unix timestamp per cycle; add root `UTC_Offset` for lab-PC wall-clock time. |
| `CALdata/Mapping` | (N, 2), N ≥ 2 | `(m/z, timebin)` anchors for mass calibration. Exactly two valid anchors are solved directly; three or more valid, well-conditioned anchors are fit by least squares. |
| `TRACEdata/TraceRaw` | (n_cyc, n_pk) | Acquisition-time pre-computed peak traces (raw cps). |
| `TRACEdata/TraceCorrected` | (n_cyc, n_pk) | Pre-computed transmission-corrected traces. |
| `TRACEdata/TraceConcentration` | (n_cyc, n_pk) | Pre-computed concentration traces (ppb). |
| `TRACEdata/TraceInfo` | (8, n_pk) | Per-trace metadata: row1 = label, row2 = centre m/z, rows3–4 = m/z window. |
| `PTR-Transmission/Masses_Factors` | (5,2,21) | `[0,0,:]` = m/z nodes, `[0,1,:]` = relative transmission factors. |
| `PTR-PrimaryIons/*` | | Primary-ion definitions (H3O⁺ monitored via the m21 H₃¹⁸O⁺ isotope ×500, etc.). |
| `AddTraces/PTR-Reaction/Data` | (n_cyc, 6) | Per-cycle drift params: `Udrift`, `p_drift`, `T-Drift_Act`, `E/N`, primary-ion index. Column names in the sibling `Info` dataset. |
| `AddTraces/DataCollection/Data` | (n_cyc, 5) | Per-cycle `ACQ_SRV_MassCal_a/b` and spec timing. |
| `AddTraces/PTR-Instrument/Data` | (n_cyc, 75) | Full instrument telemetry (voltages, temperatures, flows, turbos). |

The `viz` x-axis selector accepts exactly `cycle`, `relative`, and `absolute`.
Relative-time display uses elapsed acquisition time from `SPECdata/PCTime` when valid,
falling back to the spectrum duration when it is finite and positive, or one second
otherwise. Absolute-time display validates every `SPECdata/PCTime` value and requires a
four-digit ISO year (0000–9999) after applying root `UTC_Offset` when available. If the
dataset is missing, malformed, or outside that range, absolute display is unavailable;
relative display remains available when its domain can be represented as finite increasing
values.

## Key finding: Viewer re-processes from raw spectra

The PTR-MS Viewer CSV values are **not** copies of the pre-computed `TRACEdata`.
Viewer re-integrates the raw `SPECdata/Intensities` with its own peak list and
calibration. Evidence: the CSV's exact max values never appear in any `TRACEdata`
array, and the CSV mass labels (theoretical masses) don't match the file's
built-in trace centres. So faithful reproduction must start from the raw spectra.

## Mass calibration

TOF relation is `timebin = a·√(m/z) + b`. `CALdata/Mapping` is an `(N, 2)`
array of `(m/z, timebin)` anchors with at least two rows. With exactly two valid
anchors (finite values, positive distinct masses, and a finite positive result),
solve `a, b` directly:

```
a = (tb2 − tb1) / (√m2 − √m1)
b = tb1 − a·√m1
```

With three or more anchors, all rows must be finite and positive. After sorting by
mass, both masses and timebins must be strictly increasing. The two-column design
must have rank 2 and condition number no greater than `ε⁻¹/²`, which limits
float64 round-off amplification to roughly the square root of machine epsilon.
Only then fit `a, b` by least squares, minimising the residuals of
`timebin = a·√m + b` across all anchors. For three or more anchors, reconstruct each
anchor mass with the inverse below and accept the fit only when every reconstructed
mass has a finite absolute relative error of at most 100 ppm. This is a deliberately
generous corruption/model-consistency ceiling, not an accuracy claim. In either
case, invert the fit as:

```
m/z = ((timebin − b) / a)²
```

A valid Mapping with three or more rows is the authoritative axis. Its least-squares
`a,b` fit and per-row residuals are retained as provenance, and no later mass-domain
translation or scaling is applied. The run-average water-cluster and iodobenzene peaks
may still be measured across cycle blocks, but their movement is a temporal-stability
diagnostic only and cannot replace or modify the Mapping calibration.

If `CALdata/Mapping` is absent or unusable, Sniff falls back to usable per-cycle
`CALdata/Spectrum` coefficients. Exactly two Mapping rows follow the same fallback
calibration path. In those cases only, the sanitised run-average spectrum must contain
the operational water calibrant at 37.033 and protonated iodobenzene at 204.951. Their
sub-bin centres must be prominent, high-S/N, unambiguous, conservatively positioned and
persistent in at least five of eight deterministic cycle blocks. The accepted fallback
mapping is:

```
m_corrected = scale·m_file + offset
m_file = (m_corrected − offset) / scale
```

A failed mandatory fallback anchor is reported with its status, reason, prominence/S/N
and block persistence. Config schema version 2 migrates version-1 coordinates through
their recorded old inverse transform, preserving the selected physical timebins while
removing the obsolete second correction. Per-cycle `MassCal_a/b` also exist in
`AddTraces/DataCollection`, but barely differ from the global fit in the examined files.

Formula identification uses two deliberately different ppm limits. A broad 200 ppm
window generates formula and catalogue proposals for expert investigation. With at
least three valid Mapping rows, their reconstructed baseline-mass residuals provide an
independent check on model fit: the 95th-percentile absolute residual plus a 2 ppm
small-sample margin defines the assignment radius, clamped to 5–10 ppm. Only candidates
inside that run-specific radius may become automatic defaults. Wider candidates retain
their signed ppm errors and remain explicit reviewer hypotheses. Fewer than three
Mapping rows, or residuals above 10 ppm, disable automatic formula assignment without
hiding the broad proposal list. Block-to-block internal-reference movement remains a
temporal-stability diagnostic rather than an accuracy estimate.

Candidate discovery first tests direct `[M+H]⁺` formulas. Separate versioned ion
hypotheses may then propose hydrated/dehydrated products under measured H₃O⁺ conditions,
charge transfer or hydride abstraction when matching reagent evidence exists, and `z=2`
or `z=3` only after two fractional-spacing isotope satellites co-vary with the parent.
These alternatives retain the same ppm limits but are never assignment-eligible.
Known isotope and fragment channels can carry their parent candidate list while remaining
non-analyte channels. A versioned selected role preserves that evidence in the review
config but is invalidated when its peak is moved. Candidate-coverage categories partition
canonical components, merging neighbours that the physical-resolution diagnostic says
are unresolved. Signal-weighted coverage uses the strongest mean Raw trace per component
and is withheld if any component has no separable trace. Detector artefacts already
rejected by peak-shape diagnostics remain outside the denominator.

## The four quantities

### 1. Raw [cps]
Sum of `Intensities` over the peak's m/z window. Window is apex-centred with a
resolution-based half-width `hw = m / (2·R)`, R ≈ 1200 (≈ ±1 FWHM; the physical
resolution is R ≈ 2400, so this captures ~95 % of a Gaussian peak). Reproduces
isolated CSV peaks to <2 %.

### 2. Corrected
`Corrected = Raw / Transmission(m/z)`, transmission linearly interpolated (and
end-clamped) from the `PTR-Transmission` curve. Verified against the file's own
pre-computed traces: `TraceCorrected/TraceRaw` exactly equals `1/T(m)` with this
curve. Viewer uses a slightly different transmission curve (its own project
setting), giving a ~3–5 % systematic offset.

### 3. Concentration [ppb]
Standard primary-ion-normalised model: **`Conc = Corrected × K / I_primary(t)`**,
where `I_primary(t)` is the per-cycle reagent-ion signal (H₃O⁺ monitored via its
H₃¹⁸O⁺ isotope at m/z ≈ 21, since m/z 19 saturates) and `K` is a single
calibration constant. Dividing by `I_primary(t)` tracks reagent-ion drift over the
run. Verified: on the reference exports, `ref_sensitivity(t) × I_primary(t)` is
constant per experiment (K ≈ 18.3 and 16.6 on the two files) even though the
sensitivity itself drifts ~15 % over each run.

**K is the one irreducible calibration constant** and is *not* uniquely fixed by
the raw file — a specific PTR-MS Viewer project uses its own sensitivity setting
(K differed by ~15 % between the two experiments). Defaults:
- Default `K` is derived from the file's own pre-computed concentration
  (`K = median_t[(TraceConc/TraceCorrected)(t) × I_primary(t)]`), reproducing the
  *acquisition* calibration.
- `calibrate` fits `K` against a reference CSV
  (`K = median[ref_conc × I_primary / Corrected]`), matching a Viewer project to
  ≈ 3 %.
Raw and Corrected do not depend on K. Formula-derived isotope handling preserves both
quantities: spillover-adjusted signal affects concentration only. Dividing by the
formula's monoisotopic fraction is permitted only when calibration metadata says its
signal basis is `total`; a standards-derived or unknown K may already contain that
factor and is not silently adjusted.

### 4. Concentration [µg/m³]
`Conc_µg = Conc_ppb × M_neutral / Vₘ`, where `M_neutral = m_ion − m_proton`
(1.007276) and `Vₘ` is the molar volume at the drift/inlet temperature:
`Vₘ = 22.414 · (T_drift[K] / 273.15)` ≈ 28.9 L/mol at 80 °C. Verified: the
CSV's µg/ppb ratio equals `M_neutral/28.90` across all masses to <0.1 %.

## Overlapping peaks (isobaric interference)

The central difficulty of PTR-TOF. A window wide enough for accurate area on isolated
peaks reaches into neighbours < ~0.05 m/z away. Version-2 analyses group peaks within
`cluster_gap` (0.2 m/z), learn an empirical line shape from clean isolated peaks in the
same run, and fit bounded shared centre and width changes on the run and interval
average spectra. A dependency-free non-negative least-squares solve then obtains
amplitudes per cycle and rescales them to the window-sum Raw definition. Numerical
rank, condition, component correlation and residual diagnostics decide whether the
components are
independently identifiable. Unreliable values become unavailable rather than being
clipped into plausible concentrations. The previous fixed-centre Gaussian model remains
available as `gaussian-v1` and is the reported fallback when an empirical profile cannot
be established.

Assigned formulas also derive exact natural M+1/M+2 auxiliary channels. Schema-3
analyses build connected components wherever parent and isotope observations share a
channel, then fit all parent amplitudes jointly in transmission-corrected signal space.
The weighted non-negative solution propagates the design covariance and reports rank,
condition, residual, fitted-cycle count and median relative uncertainty. Missing,
rank-deficient, excessive-residual, negative or high-uncertainty solutions are withheld.
This changes concentration input only: Raw and Corrected remain the observed extracted
signals. Auxiliary channels are evidence and correction inputs, not independent analyte
rows. Older schema-2 reviews retain their sequential `formula-v1` correction.

## What is NOT in the raw file (must be supplied)

- **The target peak list** — which compounds to quantify.
- **The time ranges** — which cycle windows are which sample. (Here they were
  recovered from the CSV's own `Cycle` variable rows; in general they are the
  operator's annotations.)
- **Viewer's exact transmission curve and sensitivity constant** — its private
  calibration. The file's own values are used as a physically valid default.

## Open-source ecosystem

| Tool | Lang | Notes |
|---|---|---|
| **PyTRMS** (ionicon-analytik) | Python | Reads IONICON `.h5`, traces → pandas. Vendor library. |
| **ptairMS** | R/Bioconductor | Purpose-built for exhaled-breath PTR-TOF biomarker discovery; raw `.h5` → peak tables. |
| **PTRwid** | IGOR Pro | Untargeted peak detection, internal m/z calibration, deconvolution (Tofwerk instruments). |
| **PeakCalc** | — | Normalised peak areas for PTR-TOF. |
| Ionicon Data Analyzer / Tofware | commercial | Full HR peak fitting. |
