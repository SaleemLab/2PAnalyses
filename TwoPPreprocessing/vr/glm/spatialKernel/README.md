# Spatial kernels: position code or gain on vision?

## The question

In the original V1 GLM fits (Julien's VisNavGUI "GLM1"), a large group of cells
("spatial cells") have a position kernel. In that model, adding position improves
the cross-validated fit over a model with only visual and running-speed predictors.

That model is **additive**: Gaussian likelihood with an identity link, so position
adds a fixed amount of activity at each place. Two quite different things could
produce such a kernel:

1. **A real position signal:** the cell fires more at some places regardless of
   what it sees.
2. **A gain on the visual response:** the visual drive is multiplied by something
   that changes along the track, for example with running speed or with
   anticipation of a landmark. An additive model can only describe a gain as a
   position kernel, because it has no way to multiply.

The code in this folder asks two related questions:

- **Does the position effect survive a multiplicative model?** The refits use a
  Poisson log-link model, where every term multiplies the others, and add an
  explicit vision × speed interaction. If position still improves the held-out
  fit, there is a place-dependent change in activity that vision and speed
  (and their interaction) don't explain.
- **Is it additive or a gain?** The swap/omit analysis (step 3) tests this
  directly by comparing additive and gain predictions.

**Note on what the Poisson refits can and cannot show.** In a log-link model,
position also enters multiplicatively:

rate = exp(vision + speed + …) × exp(position)

So in every Poisson model here, including VSP, the position kernel is a gain
(in log units) that scales the whole predicted rate, the visual drive included.
These refits therefore test whether a place-dependent gain remains once
vision × speed is in the model. They do not, on their own, separate "a
position signal added to the visual response" from "a gain on the visual
response". That needs either step 3 or a model where position is added to the
rate under the same Poisson likelihood (not yet implemented). The only fits
here where position is additive are GLM1 and the refit-1 Gaussian fits, which
use an identity link.

## Data

- Two-photon recordings from Sonali's VR experiment: 12 sessions from 4 mice,
  listed in `sk_config.m`.
- The 3184 GLM1 spatial cells, taken from `V1_GLMData_table.csv` (`spatialcells == 1`).
- The activity is **deconvolved fluorescence**, not spike counts: it is
  non-integer, has an arbitrary scale and is overdispersed. This affects which
  likelihood-ratio test is valid (see below).
- The track is 200 cm long, with position in 2 cm bins.
- There are six trial types:
  - base;
  - two texture swaps;
  - three omits, where one element is left out.
- The visual layout is the same in every session:
  - two wall textures pass the animal at 37–43, 77–83, 117–123 and 157–163 cm;
  - six background elements repeat every 24 cm.

  `sk_landmarks.m` prints these positions from the design matrix.

Raw data and all outputs are outside the repo. The paths are set in
`sk_config.m`; outputs go to `D:\DATA\Sonali\GLM\TrialLevel` and figures to its
`figures\` subfolder.

## The models

Each model is built from blocks of predictors, using the same predictors as GLM1.

| Name | Predictors | Meaning |
|---|---|---|
| VS | onset/offset, speed, texture, background, track end | vision + speed |
| VSP | VS + position | the original GLM1 "full" model |
| VSI | VS + vision × speed | vision whose drive scales with speed |
| VSIP | VSI + position | position on top of vision × speed |

Two comparisons come out of these models:

- **Position without vision × speed:** VSP vs VS (`LLHrel`).
- **Position with vision × speed:** VSIP vs VSI (`LLHrelI`). This is the main test.

In the Poisson fits both comparisons treat position as a multiplicative gain;
they differ only in whether vision × speed is in the base model.

All refits keep as much of GLM1 as possible:
- the same sample subset;
- the **same 10 cross-validation folds**, recovered exactly by `sk_recover_folds`;
- each cell's own visual delay, read by `sk_read_best_delay`.

So any difference comes from the model, not from the data split.

## How significance is tested

Each refit reports two tests that position improves the held-out fit.

**Paired t-test over the 10 folds** (`p_*_pos`, `p_*_posI`)
- Tests whether position improves the held-out log-likelihood consistently across folds.
- A cell counts as significant if p < 0.05 and LLHrel > 1.

**Likelihood-ratio test, as in Julien's `RunGLMs.m`** (`LR_*`, `pLR_*`, `pLRq_*`)
- The statistic is 2 × the summed held-out log-likelihood gain, compared with a
  χ² distribution with one degree of freedom per position bin.
- Because the data are deconvolved activity, not Poisson counts, the raw
  statistic is too lenient.
- The code therefore also reports `pLRq_*`: the statistic divided by the
  measured dispersion φ (quasi-likelihood), where φ is the Pearson dispersion of
  the all-data fit.
- φ is typically 2–3, so the corrected test is the one to use.

## Files

### Setup and data export

| File | What it does |
|---|---|
| `sk_config.m` | Paths, the session list, condition names, the 2 cm bins and the GLM1 block order, shared by every other file. |
| `sk_export_session.m` | Reads one session's large VisNavGUI file (via `h5read`, without loading the 10 GB structs) and writes two compact files. `<session>_binned.mat` holds trial × bin activity, GLM1 predictions and occupancy. `<session>_samples.mat` holds the sample-level design matrix and activity needed for refitting. |
| `sk_load_session.m` | Loads the binned export as smoothed rates. It can swap in simulated data (`simAdd`: additive position; `simGain`: gain on vision) to check that an analysis can tell the two apart. |
| `sk_read_best_delay.m` | Reads each cell's visual delay chosen by GLM1. |
| `sk_recover_folds.m` | Recovers GLM1's trial-to-fold assignment, which isn't stored. Only one assignment of whole trials to folds reproduces the per-fold activity totals saved for every cell, so it is solved as a small integer program and then checked against all cells. |
| `sk_landmarks.m` | Prints where each texture and background element sits along the track. |

### Analyses on the binned data (steps 3–4)

These are model-free checks on the GLM1 fits, done before any refitting. Each
one is also run on simulated additive and gain data, so it is clear what each
result looks like.

| File | What it does |
|---|---|
| `sk_gain.m` | Estimates the gain profile g(x) = data / vision+speed prediction from base trials, regularised so near-silent bins don't dominate. |
| `sk_step3_swapomit.m` | Step 3. In swap and omit trials the visual drive changes but position doesn't. It predicts those trials in two ways: (a) adding the base-trial position kernel; (b) multiplying the new visual prediction by the base-trial gain. It then measures which prediction is closer (β: 0 = additive, 1 = gain). |
| `sk_step4_speed.m` | Step 4. Does the position residual scale with running speed? An additive kernel shouldn't; a gain on vision should. |
| `sk_run_steps34.m` | Runs steps 3 and 4 on real and simulated data and writes `steps34_summary.txt`. |

### Refit 1: Gaussian and Poisson, ridge (steps 1–2)

| File | What it does |
|---|---|
| `sk_refit_session.m` | Refits VS, VSP, VSI and VSIP for one session with both Gaussian (additive) and Poisson (multiplicative) likelihoods, using a ridge penalty. It saves the cross-validated scores, per-fold log-likelihoods and the position kernels. |
| `sk_run_refit.m` | Runs it over sessions (resumable) and combines them into `refit_spatialcells.csv` and `refit_kernels.mat`. |
| `sk_summarise_refit.m` | Summarises refit 1: whether position still helps once the model is Poisson and/or has vision × speed. |
| `sk_kernel_properties.m` | Checks whether earlier observations about the GLM1 kernels still hold in the refitted kernels: coverage, multiple peaks, co-location with visual peaks, ~25 cm spacing and the peak at onset. |

Refit 1 was run on 892 cells from three sessions.

### Refit 2: Poisson with ridge and elastic net (current)

Refit 2 fixes three problems found in refit 1:

1. **Omit trials.** The omit block is no longer part of the vision baseline. In
   GLM1 it is coded by distance, so it acts as a position code within omit
   trials and would compete with the position kernel. It is now fitted after
   position, on top of the fitted model.
2. **Onset and offset.** These are now 1.5 s temporal kernels built from a
   smooth basis. They replace GLM1's ~0.4 s lag block, which left the early
   part of each trial to the position kernel.
3. **Penalties.** Both ridge and **elastic net** are fitted (α = 0.9, with λ
   chosen by cross-validation separately for each model). Elastic net is closer
   to GLM1's lasso and gives sparse kernels. It uses its own glmnet-style
   solver: active-set coordinate descent with early stopping along the λ path.

| File | What it does |
|---|---|
| `sk_refit2_session.m` | The refit-2 fit for one session, with both tests (t-test and dispersion-corrected likelihood ratio). By default it fits only VSI and VSIP, because the models without vision × speed were already done in refit 1 and fitting them would double the run time. Set `opts.models = 1:4` to fit all four. |
| `sk_run_refit2.m` | Runs refit 2 over sessions in parallel and combines them into `refit2_spatialcells.csv` and `refit2_kernels.mat`. It is resumable: sessions with an existing CSV are skipped. It defaults to 4 workers, because the fits are limited by memory bandwidth and 8 or 16 workers were no faster on this machine. |

### Comparisons and figures

| File | What it does |
|---|---|
| `sk_plot_session_kernels.m` | Example kernels for one session: summary panels, the top cells that pass both tests, and the top cells that pass only the corrected likelihood-ratio test. |
| `sk_plot_session_kernel_images.m` | Sorted kernel images for one session (cells as rows): GLM1, refit 1 and refit 2. |
| `sk_plot_refit2_progress.m` | Counts and kernel images for every session refit 2 has finished, for checking on the run as it goes. |
| `sk_position_vxs.m` | Position without vs with vision × speed in the base model (Poisson VSP vs VS against VSIP vs VSI, refit 1): which cells keep a significant position term in each, and how much position improves the fit. Position is a multiplicative gain in both. |
| `sk_kernel_shapes.m` | Whether the Poisson position kernels change shape when vision × speed is added (VSP vs VSIP kernels, GLM1 for reference): depth, peak height and width, number of peaks, main peak position, and where peaks fall relative to the texture and 24 cm background landmarks. |

### Kernel-shape analysis (any kernel set)

These describe what the position kernels look like and how they relate to the rest
of each cell's activity, in a way that doesn't depend on which fit made them. When
improved kernels arrive (e.g. Julien's model), add a case to `sk_kernelset` and rerun;
every observation below can then be checked on the new kernels.

| File | What it does |
|---|---|
| `sk_cell_components.m` | Builds, once, a cache of per-cell quantities along the track that don't depend on the kernel fit, for all 6,473 cells: the measured map (and odd/even halves), GLM1's full, vision + speed and position predictions, GLM1's background, landmark and end-of-corridor contributions, the mean running speed per bin, and speed tuning. The refits only save position kernels, so GLM1 provides the vision and speed reference for every kernel set. |
| `sk_kernelset.m` | Loads a kernel set in one standard form: kernels, SEs where available, significance, cell information and the matching components. It also puts every kernel on a common scale, `G`: the fractional change in activity relative to the track mean (log kernels: exp(k − mean) − 1; additive kernels: (k − mean) / mean predicted activity). GLM1 kernels are stored per second and are converted to the per-sample scale of the maps. |
| `sk_kernel_features.m` | One row of features per cell: depth, main peak and trough, width and asymmetry of the main peak, number of peaks and troughs, roughness, onset and end effects, power at the 24 and 40 cm periods, peaks just before or after textures, and correlation/lag with the measured map, GLM1's vision and position predictions, and the speed profile. |
| `sk_run_kernel_features.m` | Runs the features on several kernel sets and writes `kernelFeatures/<set>.csv` and a side-by-side `kernel_features_summary.txt`. |
| `sk_landmark_aligned.m` | Cuts each kernel (z-scored) in ±20 cm windows around the four textures and compares the mean with a null in which each kernel is shifted circularly over 10–200 cm. It shows GLM1's vision + speed prediction, the measured map and running speed in the same windows, folds the kernels modulo 24 cm away from the textures (background period), and compares peak-to-texture distances with the null. Writes `landmark_aligned.png` and a summary. |
| `sk_plot_kernelset.m` | Population images for one kernel set: kernel, measured map, GLM1 vision + speed and GLM1 position prediction, z-scored, in the same row order (by peak, session or cluster). |
| `sk_plot_cell_cards.m` | One card per cell: measured map (with odd/even halves), GLM1 predictions, the kernel as a fractional change (± SE where available) and the mean speed. By default it shows the 16 significant cells with the deepest kernels beyond 10 cm. |
| `sk_kernel_clusters.m` | Recurring shapes: PCA and k-means (k = 6) on z-scored kernels beyond 10 cm. It reports each cluster's mean kernel with the matching measured map, vision prediction and speed, how clusters spread over sessions, and stability (centroids matched between random halves of cells and between mouse groups). |
| `sk_simulate_null.m` | Simulated activity with no position term: GLM1's vision + speed prediction, optionally shifted in time (the simulated cell responds later or earlier than the model's visual timing), with overdispersed Poisson noise at each cell's own dispersion. |
| `sk_run_sim_refit.m` | Fits the simulated cells of one session with the refit-1 pipeline at several timing shifts; writes `simrefit/sim0_*`, `simlag6_*`, `simlead6_*` (100 ms later / earlier). Load them as `sk_kernelset('sim0_vsip')` etc. |

Kernel-set names accept modifiers: `':all'` treats every cell as selected and `':s33'`
keeps one session, e.g. `'refit1_vsip:s33:all'`. `sk_run_kernel_features` and
`sk_landmark_aligned` take an output name, so comparisons don't overwrite the main outputs.

## Typical order

```matlab
for i = 1:12, sk_export_session(i); end   % once; writes the binned and samples files
sk_run_steps34                            % model-free checks (steps 3-4)
sk_run_refit; sk_summarise_refit          % refit 1 (Gaussian + Poisson, ridge)
sk_run_refit2                             % refit 2 (Poisson, ridge + elastic net); ~3-4 days for all cells
sk_plot_refit2_progress                   % while refit 2 runs
sk_position_vxs; sk_kernel_shapes         % position without vs with vision x speed
sk_cell_components; sk_run_kernel_features % kernel features for any kernel set
sk_landmark_aligned                       % kernels around the textures, with shuffle null
sk_kernel_clusters('glm1')                % recurring shapes (also draws the sorted images)
sk_plot_cell_cards('refit1_vsip', 16)     % example cells
sk_run_sim_refit(9, [0 6 -6])             % simulation control, one session (~0.5 h per shift, 8 workers)
sets = {'refit1_vsip:s33:all', 'sim0_vsip:all', 'simlag6_vsip:all', 'simlead6_vsip:all'};
sk_landmark_aligned(sets, [], 'landmark_aligned_sim'); sk_run_kernel_features(sets, 'kernel_features_sim')
```

## Results so far (26 Nov 2026)

These are provisional. Refit 2 is still running.

- **Steps 3–4.** Swap and omit trials sit between the two predictions (population
  β = 0.41, where 0 = additive and 1 = gain). On simulated data the same measure
  gives 0.01 for additive and 0.66 for gain. So the GLM1 position signal is
  partly a gain on vision and partly not. This is so far the only direct test of
  additive against gain.
- **Position without vs with vision × speed (refit 1, Poisson, 892 cells).**
  Position is a multiplicative gain in both models here.
  - Position remains significant in 74% of cells without vision × speed and 59%
    with it (59% and 48% by the corrected likelihood-ratio test).
  - Almost all cells significant with vision × speed are also significant without it.
  - Position adds about half as much once vision × speed is in the model (median
    0.0038 vs 0.0019 bits per unit activity).
- **Kernel shape (Poisson, without vs with vision × speed).**
  - The kernels are very similar (median correlation 0.92); with vision × speed
    they are about 15% shallower.
  - The peaks are just as sharp. More main peaks fall in the first 10 cm, as an
    onset effect.
  - In both, and in GLM1, pooled peak positions cluster around 75, 120 and 160 cm,
    near the textures. Counted peak by peak, though, peaks are about as common just
    after a texture as just before it (see the kernel-shape observations below).
  - Peaks also repeat with the 24 cm background period; this is weaker once
    vision × speed is in the model.
- **Refit 2 so far (2 of 12 sessions).**
  - The fraction of cells where position stays significant with vision × speed varies a lot
    between sessions: 35% in M25126_20260313 and 66% in M26004_20260318
    (elastic net, corrected likelihood-ratio test).

## Kernel-shape observations to re-test on new kernels

Medians over cells where position is significant in each fit
(`kernel_features_summary.txt`). GLM1: 3,184 cells; refit 1 Poisson without / with
vision × speed: 656 / 527 of 892; refit 2 elastic net with vision × speed: 158 of 463
(2 sessions). Each row should be rechecked when a new kernel set is added.

| Observation | GLM1 | Refit 1, no VxS | Refit 1, VxS | Refit 2 EN, VxS |
|---|---|---|---|---|
| Kernels are nearly uncorrelated with GLM1's vision + speed prediction (median r) | 0.09 | 0.07 | 0.09 | 0.06 |
| Kernels follow the measured map moderately (median r) | 0.68 | 0.64 | 0.59 | 0.60 |
| Kernels are anticorrelated with the mean running-speed profile (median r) | −0.07 | −0.12 | −0.19 | −0.28 |
| ... but hardly at all beyond 10 cm, i.e. it is the onset (median r) | −0.04 | −0.03 | −0.03 | −0.08 |
| Texture-aligned mean kernel (textures 2–4) peaks ~6 cm before the texture, dips at it, and peaks again ~4 cm after; all above the shuffle null | yes | yes | yes | yes |
| Kernel peaks 6–8 cm before a texture are over-represented (peak histogram vs null) | yes | yes | yes | yes |
| Kernels vary with a 24 cm period away from the textures (phase fold vs null) | p < 0.002 | p < 0.002 | p < 0.002 | p < 0.002 |
| Kernel clusters (beyond 10 cm) reproduce in random halves of cells (matched centroid r) | 0.90–0.99 | 0.59–0.89 | 0.54–0.88 | 0.08–0.73 |
| ... and between mouse groups (matched centroid r) | −0.10–0.88 | −0.14–0.54 | −0.18–0.50 | −0.10–0.42 |
| Peaks within 8 cm of a texture fall before : after it | 1207 : 1356 | 445 : 510 | 362 : 372 | 71 : 64 |
| Onset effect (0–10 cm vs 10–190 cm, fractional change) | 0.03 | 0.05 | 0.09 | 0.31 |
| Peaks per kernel | 2 | 4 | 4 | 2 |
| Main peak width at half height (cm) | 8 | 4 | 4 | 4 |

Notes:
- Around the textures the kernels have a consistent biphasic shape in every fit:
  a peak about 6 cm before the texture, a dip just where GLM1's vision prediction
  starts to rise, and a smaller peak about 4 cm after. Counting all peaks within
  ±8 cm gives similar numbers before and after, but the strongest over-representation
  is 6–8 cm before. Because the dip sits exactly at the visual rise, this could be
  the position term correcting the timing or width of the visual response rather than
  anticipation. The simulation control below argues against a pure artefact: the
  pipeline alone gives a weak before > after asymmetry but not the strong peaks.
- The onset transient (first 10 cm) is the largest single feature: about half the
  significant Poisson kernels have their main peak there. The animal is slowest at
  onset, which explains the negative correlation with the speed profile.
- GLM1 kernels fall into six stable shape groups beyond 10 cm (`clusters_glm1.png`):
  high at the start of the track, repeated peaks ~6 cm before each texture, a ramp
  along the track, the track end, and peaks just after the 2nd or 4th texture that
  coincide with the visual response there. How common each group is differs a lot
  between sessions and mice (e.g. the 167 cm group is mostly M26004), so the shapes
  are partly session-specific. The refit kernels are noisier and their clusters
  don't reproduce across mice (only 2 mice, 3 sessions).
- In the example cards, many large kernel peaks sit where the measured map has a peak
  that GLM1's vision + speed prediction under-predicts.
- Kernels are larger relative to the mean activity in GLM1 (median depth 1.04)
  than in the Poisson refits (0.47–0.56), and much smoother (roughness 0.2 vs 2.5–2.8).
  The 4 cm Poisson peak widths are at the bin resolution, so peak width is not
  yet a reliable feature there.

### Simulation control: cells with no position term (27 Nov 2026)

The 126 GLM1 spatial cells of M25126_20260313 were simulated from GLM1's vision + speed
prediction alone (with each cell's own overdispersion) and refitted with refit 1.
Three timings: the model's own, the cell responding 100 ms later, and 100 ms earlier.
Outputs: `kernelFeatures/kernel_features_sim.txt`, `landmark_aligned_sim.png` and summary.
One session, one noise seed.

| | Real | No position | 100 ms later | 100 ms earlier |
|---|---|---|---|---|
| Position significant, Poisson without VxS | 52/126 | 11 | 12 | 31 |
| Position significant, Poisson with VxS | 45/126 | 11 | 9 | 22 |
| Position significant, Gaussian | 33/126 | 3 | 2 | 23 |
| Kernel depth, all cells (median fractional change) | 0.34 | 0.18 | 0.19 | 0.21 |
| Correlation with the measured map (median r) | 0.56 | 0.06 | 0.07 | 0.14 |
| Same cell, simulated vs real kernel beyond 10 cm (median r) | | 0.01 | 0.03 | 0.05 |
| Textures 2–4: mean 0–8 cm before minus after (z), p vs null | 0.09, p = 0.03 | 0.11, p = 0.01 | 0.05, p = 0.12 | 0.12, p < 0.001 |
| Largest texture-aligned peak (z; null upper 95% ~0.15) | 0.37–0.75, 6–10 cm before | 0.14–0.20 | 0.14–0.29 | 0.12–0.26 |
| 24 cm phase fold: range, p vs null | 0.36, p < 0.002 | 0.20, p = 0.06 | 0.22, p = 0.01 | 0.20, p = 0.02 |

- With correct timing the pipeline gives about 9% false positives (Poisson), against
  36–41% significant in the real data. A cell responding earlier than the model's
  visual timing inflates this to 17–25%, and to 18% for the Gaussian fit; a later
  response does not. Some real "position" cells may therefore be visual-timing misfits.
- The before > after asymmetry around textures appears even with no position term, so
  on its own it is not evidence of a position signal. The strong peaks 6–10 cm before
  each texture (2–5 times the null bound) are not reproduced by any simulation.
- A weaker 24 cm structure also appears from the pipeline alone, at a different phase
  (0–8 cm vs 22 cm in the real kernels), so the 24 cm result needs this comparison
  rather than the shuffle null.
- Simulated kernels are about 60% as deep as the real ones and unrelated to each cell's
  real kernel, which follows the measured map. The real kernels are mostly not
  produced by the fitting procedure.
- To do: more sessions and seeds; repeat on Julien's improved model.

## Relation to Julien's code

Julien's VisNavGUI code (used for the original GLM1 fits) is kept locally in
`JulienCode/` at the repo root. It is in `.gitignore` because it is 1.2 GB.
**None of his files have been modified.** They were read to match the
original fits:

- **`GLMsAnalysis_Mai.m` / `RunGLMs.m`:** the GLM1 settings.
  - Gaussian lasso via glmnet, α = 1, 10 folds.
  - Position added by forward selection with a likelihood-ratio test.
  - Coefficients smoothed by one bin after fitting.
  - The visual delay chosen from 0–500 ms.
- **`createPosfile_MaiData_VisAngle4.m`:** how the visual predictors are built.
  - Visual angle = atan(distance / 6 cm), in 5° bins.
  - The omit block is binned by distance, not angle. That is why refit 2 fits it
    after position.
- **`VisSpkLinParams.mat`:** the parameters of the GLM1 run.

Only one piece of his method was taken over, and it was re-implemented here, not
copied:

- **The likelihood-ratio test from `RunGLMs.m`:** 2 × the summed held-out
  log-likelihood gain, compared with χ² with one degree of freedom per position
  bin. It is in `sk_refit2_session.m`, which adds the dispersion correction that
  deconvolved activity needs.

Two things to be aware of:

- **The code isn't the exact version used for these data.** The copy we have
  gives 22 angle bins, 13 background segments and 2 omit conditions, but the
  data have 19, 6 and 3.
- **The receptive-field mask and position-null options weren't used.** They are
  in his code, but the stored kernels show they weren't applied, so the refits
  don't use them either.
