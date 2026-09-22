# Scientific methodology audit

This audit records what Sniff can infer from an IONICON IoniTOF file, what remains a
heuristic, and what evidence is required before a result is promoted. It applies to the
methodology shipped with mass-axis algorithm 3 and the
`unresolved-component-reasons-v1` diagnostic.

## Evidence hierarchy

Sniff deliberately separates four levels:

1. **Measured signal:** calibrated m/z, fitted Raw trace, transmission-corrected trace,
   overlap diagnostics and sample/background behaviour.
2. **Formula proposal:** exact-mass arithmetic, elemental plausibility and isotope
   evidence. A broad proposal is not an assignment.
3. **Reaction-path proposal:** an explicit protonated, hydrated, dehydrated,
   charge-transfer, hydride-abstraction, fragment or charge-state hypothesis supported
   by the recorded reagent context and, where required, cycle co-variation.
4. **Compound proposal:** source-labelled names for a supported formula. Exact mass
   cannot distinguish structural isomers, and a catalogue name is not identification.

No database result, normalised score share or visually plausible concentration promotes
a result between these levels.

## Audit ledger

| Area | Literature position | Sniff implementation and verdict |
| --- | --- | --- |
| Mass calibration | PTR-ToF formula work commonly uses internal exact-mass references and reports mass accuracy in ppm. The `ptairMS` benchmark uses m/z 21.022, 203.943 and 330.84 as calibration ions. | `CALdata/Mapping` remains the file calibration. Its formula-assignment accuracy is now tested by leave-one-reference-out prediction, not residuals from the same fitted points. A degraded Mapping may receive a bounded affine correction only when exact water-cluster and iodobenzene molecular-ion references are both resolved and persistent. Automatic assignment remains disabled when held-out error exceeds 10 ppm. **Corrected.** |
| Water and iodobenzene references | H3O+(H2O) has exact m/z about 37.028405. Iodobenzene's calibration molecular ion is about 203.942993; 204.950 is the protonated neutral and 204.946 is near the M+1 region, not the usual molecular-ion calibrant. | Earlier Sniff used operational 37.033 and protonated iodobenzene 204.951. On `ptr.h5` the latter selected the M+1 satellite near 204.995. Algorithm 3 uses the exact ions and records the old coordinates only for migration. **Corrected.** |
| Peak detection | Resolution is mass-, shape-, abundance- and signal-to-noise-dependent; no fixed separation guarantees two resolvable ions. | Local maxima and physical-FWHM suppression are deterministic heuristics. Schema-4 overlap fits share one empirical shape across the run, select non-negative temporal regularisation by held-out spectral prediction and must improve on reduced models. Rank-deficient, poorly predictive or unsupported components are withheld. **Conservative heuristic.** |
| Formula search | Accurate mass and isotope patterns support molecular formulae, not structures. Kind-Fiehn-style element rules are useful filters but are not universal chemistry. | Sniff uses exact local masses, non-negative integer formulae, integer DBE and loose physical bounds. Typical VOC ratios are soft priors. The 200 ppm list is review-only; automatic assignment requires run validation at 5-10 ppm. **Supported with explicit limits.** |
| Isotope evidence | Isotope spacing and abundance can reject or support formulae, but neighbouring ions and transmission alter raw ratios. | Joint envelope fitting is non-negative and withholds rank-deficient, high-residual, high-uncertainty or boundary solutions. A satellite can now override a coincidental direct formula only with spacing, corrected abundance and primary-normalised level/change co-variation. Preview ratios remain lower-bound evidence where channels overlap. **Corrected and conservative.** |
| PTR reaction chemistry | H3O+ commonly forms MH+, but hydration, dehydration and fragmentation depend on E/N, pressure, humidity and molecule. NO+ and O2+ use different pathways. | The recorded stable reagent and E/N context gates reaction hypotheses. Alternative pathways remain separate and never assignment-eligible. They are now retained beside direct formulae; a non-H3O active reagent disables automatic `[M+H]+` assignment. **Corrected.** |
| Fragment evidence | Product-ion distributions depend strongly on operating conditions and are not MS/MS spectra. Common ions such as m/z 69 have documented interferences from larger aldehydes, alkenes and cycloalkanes. | Sniff uses condition-matched PTR Library pathways only to rerank or link existing formulae. Parent/fragment traces must co-vary after transmission correction and primary-ion normalisation. Links cannot create a formula, establish an isomer or override mass gates. **Supported.** |
| Detector echoes | Detector ringing is instrumental and must not be inferred merely from low m/z or formula absence. | Classification requires recurring calibrated delays, at least three independent response parents, stable response gain, a candidate amplitude inside the run-wide prediction interval, blocked held-out prediction and primary-normalised level/change correlation. The candidate parent family is excluded from its own model. Supporting-only patterns remain visible, and formula or saved-review evidence vetoes suppression. **Conservative response model.** |
| Quantification | The first-order equation assumes `k[M]dt << 1`, negligible reagent depletion and a known sensitivity. Rate constant alone does not account for transmission or the fraction of signal remaining as MH+. | Raw and transmission-corrected signals are measured. File-derived K reproduces a file scale; per-compound k is only approximate scaling. Reagent, artefact, background, fragment, isotope, inseparable-overlap and alternative-ion roles retain signal but receive NaN analyte concentrations. Standards calibration is required for absolute claims. **Corrected wording and role handling.** |
| Humidity | Humidity response is compound- and instrument-specific and can redistribute protonated, hydrated and fragment ions. | The power-law option is disabled by default. `p=1` is an upper-bound model, not a universal correction; absolute use requires a measured exponent and reference condition. **Acceptable optional model.** |
| Background | Published PTR workflows commonly normalise by reagent ion and subtract measured blanks when quantitative net concentrations are required. | Classification now compares transmission-corrected, primary-normalised traces. Sniff reports gross concentrations and background evidence separately; it does not silently subtract a possibly non-contemporaneous blank. Explicit provenance-bearing blank subtraction remains future work. **Conservative, documented limitation.** |
| Mass concentration | Converting ppb to micrograms per cubic metre requires a stated molar-volume temperature and pressure basis. | Sniff records whether molar volume came from drift temperature, a configured value or fallback. `Conc [ug]` follows the legacy Viewer-compatible column name; Methods provenance must be used for its basis. **Compatibility limitation.** |
| Segmentation | Plateaus are operational sample/background intervals, not chemical identities. Carry-over and monotonic backgrounds can violate simple high/low assumptions. | Gradient/level segmentation is deterministic and ranges remain human-reviewed. Background trend is diagnostic; labels stay chronological. **Acceptable heuristic.** |

## Quantification conditions

The usual simplified proton-transfer relationship is proportional to
`I(MH+) / I(H3O+) / (k dt)`. Li et al. show that it assumes first-order chemistry and
negligible primary-ion depletion, and that sensitivity also depends on ion transmission
and the fraction of product signal retained as MH+. Consequently:

- file-derived K is an instrument/file scale, not a project standards calibration;
- rate-constant scaling is approximate unless branching and sensitivity are measured;
- fragment, hydrate and alternative-ion channels need pathway-specific calibration;
- humidity-sensitive compounds remain indicative without humidity standards;
- unknown-formula channels cannot support a defensible mass concentration;
- detector saturation or material primary-ion depletion requires a standards-based or
  nonlinear method outside Sniff's present model.

## `ptr.h5` calibration consequence

The file's three Mapping points have small in-fit residuals but a 95th-percentile
leave-one-out prediction error of about 126 ppm. The run therefore cannot support
automatic 5-10 ppm assignments. Both internal exact-mass references are persistent,
so the bounded affine correction is applied to broad proposal coordinates. This adds
review candidates without changing the 200 ppm proposal gate or enabling automatic
identity assignment.

The old 204.951 reference selected the iodobenzene M+1 satellite on this run. Correcting
the reference to the 203.942993 molecular ion is the principal reason many formerly
unresolved homologous channels now enter the broad proposal window. These remain
formula and compound proposals only.

## Sources

- de Gouw, J. and Warneke, C. (2007), *Mass Spectrometry Reviews* 26, 223-257,
  <https://doi.org/10.1002/mas.20119>.
- Yuan, B. et al. (2017), *Chemical Reviews* 117, 13187-13229,
  <https://doi.org/10.1021/acs.chemrev.7b00325>.
- Pagonis, D., Sekimoto, K. and de Gouw, J. (2019), PTR reaction library,
  <https://doi.org/10.1007/s13361-019-02209-3>.
- Cappellin, L. et al. (2012), quantitative PTR-ToF-MS,
  <https://doi.org/10.1021/es203985t>.
- Li, F. et al. (2024), protonated, adduct and fragment response,
  <https://doi.org/10.5194/amt-17-2415-2024>.
- Vermeuel, M. P. et al. (2024), PTR-ToF interferences,
  <https://doi.org/10.5194/amt-17-801-2024>.
- Faraone, N. et al. (2019), formula/isotope and collision-rate analysis,
  <https://doi.org/10.5194/amt-12-5947-2019>.
- Roquencourt, C. et al. (2022), `ptairMS` benchmark and exact calibration ions,
  <https://doi.org/10.1093/bioinformatics/btac031>.
- Kind, T. and Fiehn, O. (2007), molecular formula constraints,
  <https://doi.org/10.1186/1471-2105-8-105>.
- Vlasenko, A. et al. (2010), formaldehyde humidity correction,
  <https://doi.org/10.5194/amt-3-1055-2010>.
