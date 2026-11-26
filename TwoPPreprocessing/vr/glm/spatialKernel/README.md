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

The code in this folder tests whether the GLM1 position kernels survive when the
model is allowed to be **multiplicative**. That means a Poisson log-link model,
where every term scales the others, plus an explicit vision × speed interaction.
If position still improves the held-out fit in that model, it is more than a
gain on vision.

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

- **Additive position gain:** VSP vs VS (`LLHrel`).
- **Multiplicative position gain:** VSIP vs VSI (`LLHrelI`). This is the main test.

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
| `sk_refit2_session.m` | The refit-2 fit for one session, with both tests (t-test and dispersion-corrected likelihood ratio). By default it fits only VSI and VSIP, because the additive models were already done in refit 1 and fitting them would double the run time. Set `opts.models = 1:4` to fit all four. |
| `sk_run_refit2.m` | Runs refit 2 over sessions in parallel and combines them into `refit2_spatialcells.csv` and `refit2_kernels.mat`. It is resumable: sessions with an existing CSV are skipped. It defaults to 4 workers, because the fits are limited by memory bandwidth and 8 or 16 workers were no faster on this machine. |

### Comparisons and figures

| File | What it does |
|---|---|
| `sk_plot_session_kernels.m` | Example kernels for one session: summary panels, the top cells that pass both tests, and the top cells that pass only the corrected likelihood-ratio test. |
| `sk_plot_session_kernel_images.m` | Sorted kernel images for one session (cells as rows): GLM1, refit 1 and refit 2. |
| `sk_plot_refit2_progress.m` | Counts and kernel images for every session refit 2 has finished, for checking on the run as it goes. |
| `sk_add_vs_mult.m` | Additive vs multiplicative: which cells keep a significant position term in each, and how much position improves the fit (refit 1, all four models). |
| `sk_kernel_shapes.m` | Whether the kernels change shape between additive and multiplicative models: depth, peak height and width, number of peaks, main peak position, and where peaks fall relative to the texture and 24 cm background landmarks. |

## Typical order

```matlab
for i = 1:12, sk_export_session(i); end   % once; writes the binned and samples files
sk_run_steps34                            % model-free checks (steps 3-4)
sk_run_refit; sk_summarise_refit          % refit 1 (Gaussian + Poisson, ridge)
sk_run_refit2                             % refit 2 (Poisson, ridge + elastic net); ~3-4 days for all cells
sk_plot_refit2_progress                   % while refit 2 runs
sk_add_vs_mult; sk_kernel_shapes          % additive vs multiplicative
```

## Results so far (26 Nov 2026)

These are provisional. Refit 2 is still running.

- **Steps 3–4.** Swap and omit trials sit between the two predictions (population
  β = 0.41, where 0 = additive and 1 = gain). On simulated data the same measure
  gives 0.01 for additive and 0.66 for gain. So the GLM1 position signal is
  partly a gain on vision and partly not.
- **Additive vs multiplicative (refit 1, 892 cells).**
  - Position remains significant in 74% of cells in the additive model and 59%
    in the multiplicative model (59% and 48% by the corrected likelihood-ratio test).
  - Almost all multiplicative cells are also additive cells.
  - Position adds about half as much once vision × speed is in the model (median
    0.0038 vs 0.0019 bits per unit activity).
- **Kernel shape.**
  - Multiplicative kernels are very similar to additive ones (median correlation
    0.92) and about 15% shallower.
  - The peaks are just as sharp. More main peaks fall in the first 10 cm, as an
    onset effect.
  - In both models, peaks cluster just before the textures (around 75, 120 and
    160 cm) rather than on them.
  - Peaks also repeat with the 24 cm background period; this is weaker in the
    multiplicative model.
- **Refit 2 so far (2 of 12 sessions).**
  - The fraction of cells with significant multiplicative position varies a lot
    between sessions: 35% in M25126_20260313 and 66% in M26004_20260318
    (elastic net, corrected likelihood-ratio test).

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
