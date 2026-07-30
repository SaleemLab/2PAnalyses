function plotFinalFigures_Sonali(EXP, figname, singlecell_IDs)
%%
% PLOTFINALFIGURES_SONALI  Plot summary/diagnostic figures from a
% concatenated multi-session EXP structure (GLM + Maps).
%
% EXAMPLE USAGE
%   EXP_all = load('path_to_batch_file\batch_file.mat');
%   plotFinalFigures_Sonali(EXP_all, 'Ordered-Kernels-Space-Landmarks-BG-FullMap-all')
%   plotFinalFigures_Sonali(EXP_all, 'Resp-snake-base-omit2-omit3')
%   plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', {'M26004_20260318_cell#17', 'M26004_20260318_cell#43'})
%   plotFinalFigures_Sonali(EXP_all, 'RF-LLHrel-summary')
%   plotFinalFigures_Sonali(EXP_all, 'Texture-layout-50')
%   plotFinalFigures_Sonali(EXP_all, 'LLHi-w/oSpace')

%
% INPUTS
%   EXP            Structure resulting from concatenating per-session GLM
%                   and Maps results across recording sessions (see "THE
%                   EXP STRUCTURE" section below). For 'GLM-kernels' only,
%                   EXP must come from a SINGLE session, not a
%                   concatenation (see that figname's description).
%   figname        String selecting which figure to produce. One of the
%                   options listed in "FIGNAME OPTIONS" below. Several
%                   names are matched with strcmp (exact) and others with
%                   contains (substring/prefix), see each entry.
%   singlecell_IDs  (optional, default = []) Cell selection, required only
%                   by 'GLM-kernels' and 'Resp-singleCell'. Either:
%                     - a numeric vector of row indices into EXP.Spk
%                       (i.e. indices into kernel arrays matching single cell indices), or
%                     - a cell array of cell-name strings, e.g.
%                       {'M26004_20260318_cell#281'}, which gets resolved
%                       to indices via EXP.Spk.CellListString.
%
% OUTPUT
%   None returned; the function opens MATLAB figure windows.
%
% ===========================================================================
% FIGNAME OPTIONS
% ===========================================================================
%
% 'LLHi-w/oSpace'   (exact match)
%   Paired dot plot (2 categories x 2 conditions, with connecting lines)
%   comparing each cell's LLHi (log-likelihood increase over the constant-
%   mean model, in bits/spike) for the vision-only GLM ("No space",
%   columns 1 & 3) vs. the full GLM including the spatial/position
%   predictor ("Space", columns 2 & 4). Cells are split into:
%     columns 1-2 (green): "good" cells with a significant spatial kernel (p<=0.05)
%     columns 3-4 (grey):  "good" cells without a significant spatial fit (p>0.05)
%   Only goodcells (see "CELL SELECTION MASKS" below) are included.
%
% 'GLM-kernels'   (exact match)
%   Plots the native GLM kernels exactly as fit for a SINGLE recording 
%   session - i.e. this option must be called with a per-session EXP, 
%   not a concatenated/batch EXP
%   Requires singlecell_IDs (numeric or cell-name list);
%   one row of 4 subplots is drawn per requested cell:
%     subplot 1 - "Visual kernels": the visual-space
%       kernels (x-axis = visual space, deg, 0-110), overlaid for:
%       Landmark L1, Landmark L2, Background textures kf=1..6,
%       End of corridor, and the 3 "skipped landmark" (omission) kernels 
%       for L2/L3/L4. Shaded bands = SE (via seplot).
%     subplot 2 - "VR onset/offset": Tuning(iOnOff) kernels for VR onset
%       (feature 1) and VR offset (feature 2), x-axis = time since
%       event (s), 0-0.25 s.
%     subplot 3 - "Run speed kernel": Tuning(iSpd) kernel, x-axis = run
%       speed (cm/s), 2-50.
%     subplot 4 - "Spatial position kernel": Tuning(iPos) kernel,
%     x-axis = corridor position (cm), 0-200.
%   Each subplot's y-axis is the GLM kernel value (log link-function
%   scale, hence the y=0 reference dashed line rather than y=1).
%
% 'Texture-layout-<x>'
%   Does NOT use EXP at all - purely a geometry/illustration figure of the
%   task's fixed visual layout, built from the hard-coded layout constants
%   defined at the top of the function (L1pos_pct, L2pos_pct, BGseg_pct,
%   period_pct, corrW_pct, width_pct, corridor_cm). <x> is an optional
%   trailing integer giving the viewer's position along the corridor, as a
%   percentage of corridor length (e.g. 'Texture-layout-50'); if omitted
%   or not parseable, defaults to x = 20%.
%     subplot 1 - "Background textures and landmarks in corridor's space": 
%       a 1-D strip (x = corridor position in cm, 0-200) showing the repeating
%       background-texture segments (color gradient, period = period_pct
%       of corridor length) as colored patches, with Landmark L1 (col_L1)
%       and Landmark L2 (col_L2) patches overlaid at their fixed
%       positions/widths.
%     subplot 2 - "BG textures & landmarks in visual azimuth": the same
%       layout re-projected into the visual angle (deg, 0-110) seen by a
%       viewer standing at corridor position x=<x>%, using the geometric
%       mapping azmap(pct) = 90 - atand((pct-x)/halfW_pct), where
%       halfW_pct = corrW_pct/2 is half the corridor width (in % of its
%       length). This shows how the fixed-in-space texture/landmark
%       layout maps onto the time-varying visual receptive-field axis
%       used by the GLM's visual kernels.
%
% 'Resp-snake'   (matched via contains(figname,'Resp-snake'))
%   Heatmap ("snake plot") of trial-averaged responses across corridor
%   position, for each of the 6 trial conditions in condNames = {'base',
%   'swap23','swap34','omit2','omit3','omit4'}. If figname contains one or
%   more condition-name substrings (e.g. 'Resp-snake-base-omit2') only
%   those conditions are shown; otherwise all 6 are shown.
%   Top 6x(ncond+1) grid of imagesc panels (cells x position, color = 0-1
%   min-max normalized response per cell/row), 2 row-blocks (top row: all
%   goodcells; middle row: only the spatially-selective subset of goodcells,
%   i.e. signicells & goodLLHcells), each row-block has 3 column blocks:
%     "DATA": measured mean response
%     "VS":   GLM-predicted response, Vision+Speed model, no spatial predictor
%     "VSP":  GLM-predicted response, full Vision+Speed+Position model
%   Cells are ordered once (by peak position of the base-
%   condition DATA response) and that same row order is reused across all
%   conditions/rows for visual comparison.
%   Bottom row: GLM spatial kernels (EXP.GLMs{1}.Tuning(iPos)) for the
%   spatially-selective cells, in the same row order, plus per-condition
%   residual heatmaps (DATA-VS and DATA-VSP, each min-max normalized
%   using the DATA row's own range) to visualize where the model under/
%   over-predicts relative to the data for each condition.
%
% 'Resp-singleCell'   (exact match)
%   Requires singlecell_IDs. One row of subplots per requested cell, one
%   column per trial condition in condNames (6 columns: base, swap23,
%   swap34, omit2, omit3, omit4). Each subplot overlays:
%     "Data"   - measured mean response +/- SE (EXP.Maps, black)
%     "GLM"    - full-model (VSP) GLM-predicted mean response (EXP.Maps, col_Spatial)
%   x-axis = corridor position (cm), 0-200; same y-limit shared across all
%   subplots for a given cell (set from the max across all conditions).
%
% 'RF-LLHrel-summary'   (exact match)
%   A 3x6-tile summary figure (tiledlayout) of relationships between
%   spatial coding strength (LLHrel) and several other cell properties,
%   restricted throughout to goodcells, with spatial_mask = goodcells &
%   signicells & goodLLHcells used to contrast the spatially
%   selective subset:
%     1. Histogram of peak position (cm) of the spatial GLM kernel
%        (EXP.GLMs{1}.Tuning(iPos)) vs. peak position of the raw mean
%        response (EXP.Maps, condition 'base'), spatial cells only.
%     2. Scatter of spatial-kernel peak position vs. overall-response peak
%        position (identity line for reference), spatial cells only.
%     3. Violin plot of LLHrel vs. overall-response peak position, and a
%        second violin of LLHrel vs. spatial-kernel peak position, both
%        restricted to spatial cells (peak position binned in 20 cm bins
%        from 0-200 cm).
%     4. Overlaid histograms of LLHrel distribution, non-spatial (gray)
%        vs. spatial (col_Spatial) goodcells.
%     5. Bar plot of the fraction of goodcells that are spatially
%        selective, as a function of RF position (visual azimuth, deg;
%        EXP.GLMs{1}.Perf.bestRFpos).
%     6. Violin plot of LLHrel vs. RF position, all goodcells.
%     7. Overlaid histograms of RF position distribution, non-spatial vs.
%        spatial goodcells.
%     8. Bar plot of the fraction of goodcells that are spatially
%        selective, as a function of visual response latency
%        (EXP.GLMs{1}.Perf.bestVisdelay, seconds).
%     9. Violin plot of LLHrel vs. visual latency, all goodcells.
%    10. Overlaid histograms of visual-latency distribution, non-spatial
%        vs. spatial goodcells.
%
% 'Cluster-SpaceKernels'   (exact match)
%   Hierarchical clustering (see ClusterSpaceKernels local function below)
%   of the GLM spatial kernels (EXP.GLMs{1}.Tuning(iPos)) restricted to
%   spatially-selective cells (goodcells & signicells & goodLLHcells).
%   Distance metric/linkage method are dist_metric='cityblock' /
%   dist_method='average' (set near the top of the function); number of
%   clusters is n_Clusters (=5 by default). Produces a grid figure with
%   one row per cluster x 6 columns:
%     col 1: dendrogram (spans all rows; leaf order = optimalleaforder;
%            labels = cluster ID; colored by colorCutoff)
%     col 2: cluster mean +/- SE of (optionally smoothed) min-max-scaled
%            spatial kernels
%     col 3: imagesc of the (optionally smoothed) min-max-scaled spatial
%            kernels for that cluster's cells
%     col 4: imagesc of the raw (unsmoothed) min-max-scaled spatial
%            kernels
%     col 5: imagesc of the raw (min-max-scaled) overall mean response
%            (EXP.Maps mapDATA_idx, condition 'base') for the same cells
%     col 6: histogram of recording-session ID (EXP.Spk.series, offset by
%            10*EXP.Spk.animal if EXP.Spk.animal exists) for that
%            cluster's cells, to check clusters aren't driven by a single
%            session/animal.
%   Smoothing (clust_smthwin) and PCA pre-reduction (clust_nPCs) for the
%   clustering step are both set to 0 (disabled) by default near the top
%   of the function.
%
% 'Ordered-Kernels-<...>'   (matched via contains(figname,'Ordered-Kernels'))
%   Configurable figname: a hyphen-separated list of "tags" after
%   'Ordered-Kernels-' selects (a) which kernel/map type(s) to display as
%   heatmaps, (b) the reference kernel used to order/group cells, and (c)
%   optional display modifiers. figname is split on '-' and everything
%   after the first 2 tokens ('Ordered','Kernels') is parsed as a list of
%   tags, S. The FIRST tag in S is always used as the reference (ref) for
%   ordering/grouping; ALL recognized tags present anywhere in S are drawn
%   as separate heatmap columns. Recognized tags:
%     'Space'        - GLM spatial kernel, EXP.GLMs{1}.Tuning(iPos)
%     'Landmarks'    - GLM landmark kernel, EXP.GLMs{1}.Tuning(itex)
%     'BG'           - GLM background-texture kernel, EXP.GLMs{1}.Tuning(iBG)
%     'EOC'          - GLM end-of-corridor kernel, EXP.GLMs{1}.Tuning(iEOC)
%     'Vis'          - sum of the Landmarks + BG kernels
%     'FullMap'      - full-model (VSP) predicted response, base
%                      condition, EXP.Maps mapVSP_idx(1)
%     'Omit'         - the 3 omission-component GLM kernels (L2/L3/L4
%                      skipped-landmark responses), EXP.GLMs{1}.Tuning(iOmit)
%     'OmitMap'      - measured mean response (EXP.Maps mapDATA_idx group)
%                      for the omit2/omit3/omit4 conditions
%     'DataOmitMapDiff'   - measured-response difference, each omission
%                      condition minus base (EXP.Maps)
%     'VisOmitMapDiff'    - same difference but for the VS (Vision+Speed,
%                      no Position) model prediction (EXP.Maps)
%     'noOmitFullOmitMapdiff' - same difference but for the full-model
%                      prediction fit WITHOUT the omission predictor (EXP.Maps)
%     'swapMap'      - measured mean response for the swap23/swap34
%                      conditions (EXP.Maps)
%     'DataSwapMapDiff' - measured-response difference, each swap
%                      condition minus base (EXP.Maps)
%   Other (non-kernel-selecting) tags modify behavior rather than adding a
%   column:
%     'all'          - include ALL goodcells (default, when neither 'all'
%                      nor 'omission' is given, restricts to spatially
%                      selective cells: goodcells & signicells &
%                      goodLLHcells)
%     'omission'     - restrict to goodOmitcells instead (goodcells with a
%                      significant, LLH-improving omission component; see
%                      "CELL SELECTION MASKS"). Also switches the
%                      color-coded "spatial weight" sidebar to use
%                      LLHrel_omit instead of LLHrel.
%     'hist'         - add a per-group session-ID histogram next to each
%                      group's mean-kernel subplot
%     'peak'         - order/group cells by peak position of the
%                      reference kernel instead of the default similarity-
%                      based seriation (orderKernelsBySimilarity, which
%                      minimizes the sum of distances between successive
%                      cells in the ordering)
%   Example: 'Ordered-Kernels-Space-Landmarks-BG-hist' draws Space (used
%   as the ordering reference), Landmarks and BG kernel heatmaps,
%   similarity-ordered, restricted to spatially-selective cells, with
%   per-group session histograms.
%   Layout: leftmost column = a per-cell "spatial weight" strip (LLHrel or
%   LLHrel_omit, red ticks = spatially-selective cells, black = not, gray
%   line = local running median) plus a small colored scatter strip
%   showing each cell's recording session; then one heatmap column per
%   requested kernel/map tag (cells in rows, reordered; position in cm on
%   x-axis), with the reference column additionally annotated with
%   colored group brackets/labels; then, per contiguous "group" of Kn(1)-
%   Kn(2) similar adjacent cells (found by sliding-window SSE
%   minimization on the z-scored reference kernel, up to max_examples
%   groups), a small "mean +/- SE kernel" subplot for every requested
%   kernel/map tag overlaid (reference kernel highlighted in its own
%   color; non-reference kernels shown for shape comparison) and,
%   if 'hist' was requested, an adjacent histogram of session ID for that
%   group's cells.
%
% ===========================================================================
% CELL SELECTION MASKS (computed once near the top of the function, used
% throughout most figname options)
% ===========================================================================
%   goodcells       Cells passing basic quality control:
%                     - EXP.Maps{1}.Tuning.nlogL_pval <= p_th (the overall
%                       mean-response GLM fit, base condition, is
%                       significantly better than a constant-rate model)
%                     - EXP.GLMs{1}.Tuning(itex).pval <= p_th (the
%                       vision model fit is significant)
%                     - the landmark kernel (Tuning(itex)) is not flat
%                       (max ~= min across bins, i.e. the GLM didn't
%                       degenerate to a constant for this predictor)
%                   p_th = 0.05 (set near the top of the function).
%   signicells      Cells with a statistically significant spatial-
%                   position kernel: EXP.GLMs{1}.Tuning(iPos).pval<=p_th.
%   goodLLHcells    Cells for which adding the spatial predictor improves
%                   the model enough to matter: LLHrel >= 1+LLHrel_th
%                   (LLHrel_th = 0.015, set near the top of the function).
%   spatialcells    goodcells & signicells & goodLLHcells - the canonical
%                   "spatially selective" cell set (referred to as
%                   spatial_mask in several figname blocks).
%   goodOmitcells   Cells with a significant AND LLH-improving omission
%                   component: a likelihood-ratio test (lratiotest)
%                   comparing the best model (with omission terms) against
%                   the same model without them (32 degrees of freedom =
%                   16 omission parameters x 2, i.e. comparing the "no
%                   omit" reduced model to the full model), combined with
%                   LLHrel_omit > 1+LLHrel_th.
%
% ===========================================================================
% LLH / LLHi / LLHrel NOMENCLATURE
% ===========================================================================
%   nlogL    Negative log-likelihood of the fitted GLM (EXP.GLMs{1}.Perf.
%            nlogL), one column per model variant (see "GLM MODEL COLUMN
%            INDICES" below), one row per cell.
%   LLH      Log-likelihood (= -nlogL for a given model column).
%   LLHi     "LLH increase": the improvement in log-likelihood (in bits/
%            spike) of a given model over the constant-mean (null) model,
%            i.e. LLHi = -(nlogL(model) - nlogL(null)) / nSpikes / log(2).
%            nlogL(:,1) is the null/constant model; nSpikes is the spike
%            count used for the fit (EXP.GLMs{1}.Perf.nSpikes(:,1)).
%   LLHrel   "LLH relative": the ratio of LLHi for the full model
%            (including the spatial/Position predictor) over LLHi for the
%            vision-only model (without it): LLHrel = LLHi_pos / LLHi_vis.
%            LLHrel > 1 means adding the spatial predictor further
%            improves the fit beyond what vision alone explains.
%   Variants of LLHi/LLHrel (LLHi_noOmit, LLHrel_omit, LLHi_noL1L2,
%   LLHrel_L1L2, LLHi_noBG, LLHrel_BG, LLHi_noSpd, LLHrel_Spd,
%   LLHi_noOnOff, LLHrel_OnOff) follow the same logic but knock out one
%   predictor group at a time (omission terms, landmarks L1/L2, BG
%   textures, running speed, VR onset/offset respectively) from whichever
%   model (vision-only or full) is "best" for that cell (i.e. the
%   spatial/non-spatial model selected per-cell via spatialcells), to
%   quantify how much each predictor group contributes.
%
% ===========================================================================
% GLM MODEL COLUMN INDICES (columns of EXP.GLMs{1}.Perf.nlogL)
% ===========================================================================
%   The GLM is fit as a sequence of nested models, each adding/removing
%   predictor groups, with each model's nlogL stored in one column. The
%   function hard-codes the column indices for the specific model-fitting
%   sequence used in this pipeline:
%     column 1            : null / constant-mean model
%     column iVismodel=121 : full vision-only model (no spatial predictor)
%     column iVisnoOmit=100   : vision-only model, omission terms removed
%     column iVisnoL1L2=106   : vision-only model, L1/L2 landmark terms removed
%     column iVisnoBG=103     : vision-only model, BG-texture terms removed
%     column iVisnoSpd=110    : vision-only model, running-speed term removed
%     column iVisnoOnOff=115  : vision-only model, VR onset/offset term removed
%     column end (last)   : full model including the spatial/Position
%                            predictor
%     column iPosnoOmit=122   : full model, omission terms removed
%     column iPosnoL1L2=125   : full model, L1/L2 landmark terms removed
%     column iPosnoBG=124     : full model, BG-texture terms removed
%     column iPosnoSpd=126    : full model, running-speed term removed
%     column iPosnoOnOff=127  : full model, VR onset/offset term removed
%   These indices are specific to this GLM-fitting pipeline's predictor
%   ordering and will need updating if that pipeline's model sequence
%   changes.
%
% ===========================================================================
% THE EXP STRUCTURE
% ===========================================================================
% EXP is built by concatenating per-session results across recording
% sessions (except where noted). Only the fields actually populated/used
% by this script are documented below.
%
%    EXP.Spk                 1x1 struct. Per-cell metadata/identity,
%                            concatenated across sessions (one row/entry
%                            per cell across the whole concatenated
%                            dataset).
%     EXP.Spk.CellListString   Cell array of strings, one per cell, e.g.
%                              'M26004_20260318_cell#281'. Used to resolve
%                              cell-name inputs in singlecell_IDs to row
%                              indices.
%     EXP.Spk.series           Recording session ID (date) per cell.
%     EXP.Spk.animal           Animal ID per cell; combined
%                              with series as sessionid = series +
%                              10*animal to give a unique session
%                              identifier across animals (used in
%                              'Cluster-SpaceKernels' and
%                              'Ordered-Kernels-*').
%
%      EXP.GLMs              1x1 cell array (a single cell, EXP.GLMs{1},
%                            even though sessions are concatenated - the
%                            GLM results themselves are concatenated along
%                            the cell/row dimension inside that one
%                            struct, not as separate cell-array entries).
%      EXP.GLMs{1}.Tuning      Struct array, one element per GLM predictor
%                              group, in this fixed order (indices set
%                              near the top of the function):
%                                1 (iOnOff) : VR onset/offset
%                                2 (iSpd)   : running speed
%                                3 (itex)   : landmarks (L1, L2)
%                                4 (iBG)    : background textures
%                                5 (iEOC)   : end of corridor
%                                6 (iOmit)  : omission components (skipped
%                                             landmarks L2/L3/L4)
%                                7 (iPos)   : spatial/position kernel
%                                             (iPos = numel(Tuning), i.e.
%                                             always the LAST entry)
%                              Each Tuning(k) element has fields:
%       .meanrespModel  [cells x features x bins] (and a 4th
%                        dim for conditions in case of visual kernels) - 
%                       the fitted kernel value per cell,
%                       per feature column within that predictor group
%                       (e.g. feature 1/2 = L1/L2 for itex; feature
%                       1/2 = onset/offset for iOnOff; feature kf=1..6 for
%                       the 6 BG texture segments at iBG; feature
%                       1/2/3 = omit-L2/omit-L3/omit-L4 for iOmit), per
%                       spatial/temporal bin of that kernel's own x-axis
%                       (e.g. visual space 0-110 deg for itex/iBG/iEOC/
%                       iOmit; time 0-0.25s for iOnOff; speed 2-50 cm/s
%                       for iSpd; corridor position 0-200cm for iPos).
%       .SErespModel    Same shape as meanrespModel; standard error of the
%                       kernel estimate.
%       .pval           [cells x 1] significance of that predictor group's
%                       contribution to the fit.
%     EXP.GLMs{1}.Perf         Struct of overall per-cell GLM fit
%                              statistics:
%       .nlogL          [cells x model-columns] negative log-likelihood,
%                       one column per nested model variant (see "GLM
%                       MODEL COLUMN INDICES" above).
%       .nSpikes        [cells x >=1] spike count(s) used for the fit;
%                       column 1 used throughout this script.
%       .bestRFpos      [cells x 1] best-fit visual receptive-field
%                       position (azimuth, deg).
%       .bestVisdelay   [cells x 1] best-fit visual response latency (s).
%
%   EXP.Maps                 1x42 cell array. Each cell EXP.Maps{i} is a
%                            struct holding a [cells x bins] (or
%                            [cells x features x bins]) trial-averaged
%                            response map, computed either from the raw
%                            data or from a GLM model's predicted output,
%                            for ONE trial condition. The 42 cells are 7
%                            consecutive groups of 6, one entry per
%                            condition in condNames = {'base','swap23',
%                            'swap34','omit2','omit3','omit4'} (in that
%                            order) within each group:
%                              cells  1- 6 (mapDATA_idx)      : measured data (actual recorded mean response)
%                              cells  7-12 (mapVSP_idx)       : full GLM model prediction (Vision + Speed + Position)
%                              cells 13-18 (mapVS_idx)        : GLM model prediction without the spatial/Position predictor (Vision + Speed only)
%                              cells 19-24                    : landmark component of the model output
%                              cells 25-30                    : BG-texture component of the model output
%                              cells 31-36                    : omission component of the model output
%                              cells 37-42 (mapVSPnoOmit_idx) : full GLM model prediction with the omission predictor excluded
%                            
%                            Each EXP.Maps{i} struct has field:
%       .Tuning         Struct with fields:
%         .meanrespModel  [cells x bins] (or [cells x features x bins]) -
%                         trial-averaged response (data or model
%                         prediction) across corridor position (0-200cm)
%                         for that condition.
%         .SErespModel    Same shape; standard error.
%         .nlogL_pval     [cells x 1] significance of the
%                         overall mean-response fit vs. a constant-rate
%                         model; one of the goodcells QC criteria.
%
% IMPORTANT DISTINCTION - EXP.GLMs{1}.Tuning kernels vs. EXP.Maps:
%   The GLM Tuning kernels (EXP.GLMs{1}.Tuning) are the model's fitted
%   per-predictor response functions (e.g. "how does firing rate depend
%   on corridor position, all else being equal"), each expressed along
%   that predictor's own native axis (position, visual space, time,
%   speed). EXP.Maps are reconstructed trial-averaged response curves
%   along corridor position only, either from the real data or by
%   plugging a given subset of GLM kernels back through the model
%   (including, for the visual/BG/landmark kernels, the time-varying
%   visual latency that maps corridor position to visual space over the
%   course of a trial) - so EXP.Maps incorporate the visual latency
%   integration that the native Tuning kernels do not.
%
% ===========================================================================
% EXTERNAL / HELPER FUNCTIONS REQUIRED ON PATH
% ===========================================================================
%   seplot          Plots a shaded standard-error band around a mean trace
%                   (seplot(ax, x, mean, se, alpha)).
%   RedWhiteBlue    Custom diverging (red-white-blue) colormap function.
%   lratiotest      Likelihood-ratio test (MATLAB Econometrics Toolbox).
%   SpecialSmooth   (optional) smoothing function; used only if present on
%                   path and if clust_smthwin > 0.
%   
% Local helper functions defined within this file (see end of file):
%   plot_cat                     - paired dot/line plot for one category
%   plotViolin                   - per-bin kernel-density violin plot
%   ClusterSpaceKernels          - hierarchical clustering of kernels
%   orderKernelsBySimilarity     - similarity-based seriation of kernels

%% HANDLING SINGLECELL_IDS ARGUMENT IF PROVIDED
if nargin < 3
    singlecell_IDs = [];% no cell selection requested (fine for figname
    % options that don't need one, e.g. 'LLHi-w/oSpace', 'Resp-snake')
end

% Allow singlecell_IDs to be given either as numeric row indices (used
% as-is below) or as a cell array of cell-name strings (e.g.
% {'M26004_20260318_cell#281'}); in the latter case, resolve each name to
% its row index by matching against EXP.Spk.CellListString.
if iscell(singlecell_IDs)
    mask = cellfun(@(a) any(strcmp(a, singlecell_IDs)), EXP.Spk.CellListString);
    singlecell_IDs = find(mask);
end

%% THIS BLOCK ARE PARAMETERS YOU MAY CHANGE IF NEEDED
% THE OTHER PARAMETERS ARE TIGHTLY COUPLED TO THE ANALYSIS PIPELINE OR THE
% VR GEOMETRY
% ---------------------------------------------------------------------
% Tunable thresholds / parameters used throughout the function
% ---------------------------------------------------------------------
LLHrel_th = 0.015; % minimum LLHrel margin above 1 required to call a cell
                    % "spatial" (part of goodLLHcells below); guards
                    % against cells that pass significance only by a
                    % razor-thin margin
p_th = 0.05;        % significance threshold (p-value) used for all GLM
                    % predictor significance tests below (goodcells,
                    % signicells, goodOmitcells)
                    
% --- 'Cluster-SpaceKernels' clustering options ---
n_Clusters = 5;% 10            % number of clusters requested (k-means)
clust_nPCs = 0;% 20            % # PCs to reduce kernels to before
                                % clustering; 0 = skip PCA, cluster on the
                                % (z-scored) kernels directly
clust_smthwin = 0;%0;          % smoothing window (% of kernel length)
                                % applied before clustering; 0 = no smoothing

% --- 'Ordered-Kernels-*' grouping options ---
Kn = [4 4];%[1 1];% size range [min max] of each contiguous group of
                   % "most similar adjacent cells" averaged after ordering
max_examples = 9;%100;  % cap on the number of groups formed (so the figure
                         % doesn't grow unbounded for large cell counts)
                         

% Distance metric/linkage method shared by both the 'Cluster-SpaceKernels'
% hierarchical clustering and the 'Ordered-Kernels-*' similarity ordering
dist_metric = 'cityblock';%'cosine';%
dist_method = 'average';

% Scatter-point cosmetics (size, transparency) reused across scatter plots
ptSize = 40;
ptAlpha = 0.5;

%% VR GEMOETRY FOR TEXTURE FIGURE
% ---------------------------------------------------------------------
% Fixed task/corridor geometry (used only by the 'Texture-layout-<x>'
% figure, which does not read EXP at all - these constants describe the
% experimental rig, not anything inferred from the data)
% ---------------------------------------------------------------------
corridor_cm = 200;% full length in cm
corrW_pct  = 6; % corridor width in percent of corridor's length
width_pct = 4; % landmark's width in percent 

L1pos_pct = [20 60];% Landmark 1 (gratings) positions
L2pos_pct = [40 80];% Landmark 2 (plaids) positions
BGseg_pct = [(0:2:10)' (2:2:12)'];% BG segment start/end in percent
period_pct = 12;% BG periodicity in percent


%% CODE BOOK TO EASILY ACCESS THE MAPS CORRESPONDING TO SPECIFIC COMPONENTS OR TRIAL CONDITIONS
% ---------------------------------------------------------------------
% Trial-condition bookkeeping and indices into EXP.Maps (1x42 cell array)
% ---------------------------------------------------------------------
% EXP.Maps is organized as 7 consecutive blocks of 6 cells (one cell per
% condition, in the order given by condNames). Each *_idx variable below
% is a vector of 6 linear indices into EXP.Maps picking out one block; to
% get the map for condition icond within block X, index EXP.Maps{X_idx(icond)}.
% Only 4 of the 7 blocks are read by this script:
%   mapDATA_idx      block 1 (cells  1- 6): measured data
%   mapVSP_idx       block 2 (cells  7-12): full model (Vision+Speed+Position)
%   mapVS_idx        block 3 (cells 13-18): reduced model (Vision+Speed only)
%   mapVSPnoOmit_idx block 7 (cells 37-42): full model fit without the
%                    omission predictor

%condNames = {'base', 'swap23', 'swap34', 'omit2', 'omit3', 'omit4'}; %condition names (same order as condition ID)
condNames = {'base', 'swap23', 'omit2', 'omit3'}; %condition names (same order as condition ID)
ncond = numel(condNames);
mapDATA_idx = 1:ncond;
mapVS_idx = (1:ncond) + 2*ncond;
mapVSP_idx = (1:ncond) + 1*ncond;
mapVSPnoOmit_idx = (1:ncond) + 6*ncond;

%% GLM CODES FOR EACH PREDICTOR
% ---------------------------------------------------------------------
% Indices into the EXP.GLMs{1}.Tuning struct array (one element per GLM
% predictor group, see header doc for the full list/meaning)
% ---------------------------------------------------------------------
iOnOff = 1;
iSpd = 2;
itex = 3;
iBG = 4;
nBG = 6;    % number of background-texture feature columns within Tuning(iBG)
iEOC = 5;
iOmit = 6;
iPos = numel(EXP.GLMs{1}.Tuning); % spatial/position kernel is always the
                                   % LAST entry of the Tuning array

%% MODEL CODE FOR RETRIEVING LLH
% ---------------------------------------------------------------------
% Column indices into EXP.GLMs{1}.Perf.nlogL, selecting specific nested
% GLM model variants used to compute LLHi/LLHrel further below. These are
% fixed positions in this particular GLM-fitting pipeline's model
% sequence (see header doc, "GLM MODEL COLUMN INDICES") and would need
% updating if that pipeline's predictor/model ordering changes.
% "Vis*" = vision-only models (no spatial/Position predictor);
% "Pos*"  = full models (including the spatial/Position predictor);
% "no<X>" = that predictor group removed from the (vis/pos) model.
% Those indices can be retrieved from EXP_all.GLMs{1}.Perf.AllmodelsCode,
% which lists predictor combinations of all models in the same order
% ---------------------------------------------------------------------
iVismodel = 121;
iVisnoOmit = 100;
iVisnoL1L2 = 106;
iVisnoBG = 103;
iVisnoSpd = 110;
iVisnoOnOff = 115;
iPosnoOmit = 122;
iPosnoL1L2 = 125;
iPosnoBG = 124;
iPosnoSpd = 126;
iPosnoOnOff = 127;

%% RELATIVE LLH INCREASE OF THE DIFFERENT MODELS
% ---------------------------------------------------------------------
% LLHi (log-likelihood increase over the null/constant model, in
% bits/spike) and LLHrel (ratio of full-model to vision-only-model LLHi)
% for the "headline" full vs. vision-only model comparison, plus the same
% pair recomputed with one predictor group knocked out at a time
% (Omit / L1L2 / BG / Spd / OnOff), to quantify each group's contribution.
% See header doc "LLH / LLHi / LLHrel NOMENCLATURE" for full definitions.
% ---------------------------------------------------------------------
nlogL = EXP.GLMs{1}.Perf.nlogL;
nSpikes =  EXP.GLMs{1}.Perf.nSpikes(:,1);
% LLHi_vis: improvement of the full vision-only model over the null model
LLHi_vis = -(nlogL(:, iVismodel) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_vis(LLHi_vis == 0) = NaN; % avoid Inf when later dividing by this
% LLHi_pos: improvement of the full model (incl. spatial/Position
% predictor, stored in the LAST nlogL column) over the null model
LLHi_pos = -(nlogL(:, end) - nlogL(:, 1)) ./ nSpikes / log(2);
% LLHrel > 1 means adding the spatial predictor explains more than vision alone
LLHrel  = LLHi_pos ./ LLHi_vis;

% --- Same LLHi pair, but with the omission predictor group removed from
% each model, used to test whether the omission terms specifically help ---
LLHi_visnoOmit = -(nlogL(:, iVisnoOmit) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_visnoOmit(LLHi_visnoOmit == 0) = NaN; % avoid Inf
LLHi_posnoOmit = -(nlogL(:, iPosnoOmit) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_posnoOmit(LLHi_posnoOmit == 0) = NaN;

% --- Same pattern, landmarks (L1/L2) removed ---
LLHi_visnoL1L2 = -(nlogL(:, iVisnoL1L2) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_visnoL1L2(LLHi_visnoL1L2 == 0) = NaN; % avoid Inf
LLHi_posnoL1L2 = -(nlogL(:, iPosnoL1L2) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_posnoL1L2(LLHi_posnoL1L2 == 0) = NaN;

% --- Same pattern, BG textures removed ---
LLHi_visnoBG = -(nlogL(:, iVisnoBG) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_visnoBG(LLHi_visnoBG == 0) = NaN; % avoid Inf
LLHi_posnoBG = -(nlogL(:, iPosnoBG) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_posnoBG(LLHi_posnoBG == 0) = NaN;

% --- Same pattern, running speed removed ---
LLHi_visnoSpd = -(nlogL(:, iVisnoSpd) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_visnoSpd(LLHi_visnoSpd == 0) = NaN; % avoid Inf
LLHi_posnoSpd = -(nlogL(:, iPosnoSpd) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_posnoSpd(LLHi_posnoSpd == 0) = NaN;

% --- Same pattern, VR onset/offset removed ---
LLHi_visnoOnOff = -(nlogL(:, iVisnoOnOff) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_visnoOnOff(LLHi_visnoOnOff == 0) = NaN; % avoid Inf
LLHi_posnoOnOff = -(nlogL(:, iPosnoOnOff) - nlogL(:, 1)) ./ nSpikes / log(2);
LLHi_posnoOnOff(LLHi_posnoOnOff == 0) = NaN;

% Raw log-likelihoods (not normalized by spike count) for the vision
% and full models, needed below for the likelihood-ratio
% test (lratiotest requires LLH, not the normalized LLHi).
LLH_vis = -nlogL(:, iVismodel);
LLH_pos = -nlogL(:, end);
LLH_visnoOmit = -nlogL(:, iVisnoOmit);
LLH_posnoOmit = -nlogL(:, iPosnoOmit);

%% CELL SELECTION
% ---------------------------------------------------------------------
% goodcells: basic quality-control mask, independent of whether a cell
% turns out to be spatially selective. A cell must have (1) a
% significantly-better-than-constant overall response fit, (2) a
% significant GLM vision model fit, and (3) non-degenerate (non-flat)
% landmark and BG kernels (i.e. the GLM didn't collapse those predictors
% to a constant during fitting, which would indicate a failed/crashed fit
% for that cell rather than a genuine lack of texture sensitivity).
% ---------------------------------------------------------------------
goodcells = EXP.Maps{1}.Tuning.nlogL_pval <= p_th &  EXP.GLMs{1}.Tuning(itex).pval<=p_th;
% removing cells for which the GLM crashed on overall resp fit
X = squeeze(EXP.GLMs{1}.Tuning(itex).meanrespModel(:,1,:,1));
mx = max(X, [], 2);
mn = min(X, [], 2);
goodcells = goodcells & mx ~= mn; % landmark (L1) kernel must not be flat
X = squeeze(EXP.GLMs{1}.Tuning(iBG).meanrespModel(:,1,:,1));
mx = max(X, [], 2);
mn = min(X, [], 2);
goodcells = goodcells & mx ~= mn; % BG-texture kernel must not be flat

% signicells: cells whose spatial/position kernel fit is itself
% statistically significant (independent of effect size)
signicells = EXP.GLMs{1}.Tuning(iPos).pval <= p_th;
% goodLLHcells: cells for which the spatial predictor's contribution is
% large enough to matter in practice (LLHrel comfortably above 1), not
% just statistically detectable
goodLLHcells = LLHrel >= 1 + LLHrel_th;
        

% spatialcells ("spatially selective"): the conjunction of all three
% criteria above. Used below to decide, on a per-cell basis, which model
% ("with space" vs "vision-only") is the right reference model when
% computing the predictor-knockout LLHi/LLHrel variants - i.e. a cell
% that isn't spatially selective shouldn't be judged against the
% spatial model's LLHi.
spatialcells = goodcells & signicells & goodLLHcells;


% LLHi_best / LLH_best: per-cell LLHi/LLH picked from whichever model
% (full-with-space, or vision-only) is appropriate for that cell, based
% on spatialcells. This makes the downstream "knockout" LLHrel_* values
% comparable across both spatial and non-spatial cells using each cell's
% own most appropriate baseline model.
LLHi_best = NaN(size(LLHi_pos));
LLHi_best(spatialcells) = LLHi_pos(spatialcells);
LLHi_best(~spatialcells) = LLHi_vis(~spatialcells);

% LLHrel_omit: how much the omission predictor group contributes,
% relative to LLHi_best, computed per-cell using the same "appropriate
% model" logic as above.
LLHi_noomit = NaN(size(LLHi_pos));
LLHi_noomit(spatialcells) = LLHi_posnoOmit(spatialcells);
LLHi_noomit(~spatialcells) = LLHi_visnoOmit(~spatialcells);
LLHrel_omit = LLHi_best ./ LLHi_noomit;


% LLHrel_L1L2: contribution of the landmark (L1/L2) predictor group
LLHi_noL1L2 = NaN(size(LLHi_pos));
LLHi_noL1L2(spatialcells) = LLHi_posnoL1L2(spatialcells);
LLHi_noL1L2(~spatialcells) = LLHi_visnoL1L2(~spatialcells);
LLHrel_L1L2 = LLHi_best ./ LLHi_noL1L2;

% LLHrel_BG: contribution of the BG-texture predictor group
LLHi_noBG = NaN(size(LLHi_pos));
LLHi_noBG(spatialcells) = LLHi_posnoBG(spatialcells);
LLHi_noBG(~spatialcells) = LLHi_visnoBG(~spatialcells);
LLHrel_BG = LLHi_best ./ LLHi_noBG;

% LLHrel_Spd: contribution of the running-speed predictor
LLHi_noSpd = NaN(size(LLHi_pos));
LLHi_noSpd(spatialcells) = LLHi_posnoSpd(spatialcells);
LLHi_noSpd(~spatialcells) = LLHi_visnoSpd(~spatialcells);
LLHrel_Spd = LLHi_best ./ LLHi_noSpd;

% LLHrel_OnOff: contribution of the VR onset/offset predictor
LLHi_noOnOff = NaN(size(LLHi_pos));
LLHi_noOnOff(spatialcells) = LLHi_posnoOnOff(spatialcells);
LLHi_noOnOff(~spatialcells) = LLHi_visnoOnOff(~spatialcells);
LLHrel_OnOff = LLHi_best ./ LLHi_noOnOff;

% Same "pick the appropriate model per cell" logic as LLHi_best above,
% but on the raw (non-normalized) LLH, needed for the likelihood-ratio
% test just below (lratiotest compares raw log-likelihoods, not LLHi).
LLH_best = NaN(size(LLH_pos));
LLH_best(spatialcells) = LLH_pos(spatialcells);
LLH_best(~spatialcells) = LLH_vis(~spatialcells);
LLH_noomit = NaN(size(LLH_pos));
LLH_noomit(spatialcells) = LLH_posnoOmit(spatialcells);
LLH_noomit(~spatialcells) = LLH_visnoOmit(~spatialcells);

% Likelihood-ratio test: is the omission-component model significantly
% better than the same model without omission terms? 16*3=48 degrees of
% freedom (16 omission parameters in each of the three omission component).
[~, pval_omit] = lratiotest(LLH_best, LLH_noomit, 16*3);


% goodOmitcells: cells whose omission component is both statistically
% significant (LR test) AND large enough to matter in practice (LLHrel_omit
% comfortably above 1) - the "Omit"-analog of spatialcells.
goodOmitcells = pval_omit<= 0.05 & LLHrel_omit > 1 + LLHrel_th;


%% AESTHETICS
% ----- Colors (normalize 0-1) -----
% Centralized color palette so the same component (e.g. Landmark L1) is
% drawn in the same color across every figname option that shows it.
col_Spatial        = [145 145 228]/255;   % Color for spatial components, metrics, etc

col_L1        = [143 190 247]/255;   % Landmarks L1
col_L2        = [221 143 247]/255;   % Landmarks L2
col_BG_start  = [247  33  10]/255;   % BG textures start (gradient anchor 1, used with col_BG_end below)
col_BG_end    = [247 242 143]/255;   % BG textures end   (gradient anchor 2)
col_end       = [143 105 156]/255;   % End of corridor
col_omit2     = [ 12  95 196]/255;   % Skipped landmark #2
col_omit3     = [124  11 161]/255;   % Skipped landmark #3
col_omit4     = [11  124 161]/255;   % Skipped landmark #4
col_vr_on     = [ 63 123 158]/255;   % VR onset
col_vr_off    = [128 197 237]/255;   % VR offset
col_speed     = [0.45 0.45 0.45];    % Run speed
col_sim_pos   = [0.6 0.6 0.6];       % Spatial position (simulated overlay) -
                                      % NOTE: currently unused; left over
                                      % from a simulated-ground-truth
                                      % overlay that has been removed from
                                      % this version of the script

% Gradient for the BG textures: nBG+1 colors interpolated linearly from
% col_BG_start to col_BG_end (one color per background-texture segment,
% giving a visible "depth" gradient as the texture repeats down the
% corridor).
bg_cols = [linspace(col_BG_start(1), col_BG_end(1), nBG + 1)', ...
           linspace(col_BG_start(2), col_BG_end(2), nBG + 1)', ...
           linspace(col_BG_start(3), col_BG_end(3), nBG + 1)'];
       
       
%% SIMPLE FIGURE TO SHOW INCREASE IN LLH WHEN ADDING SPATIAL PREDICTOR
if strcmp(figname, 'LLHi-w/oSpace')
    % ====== LLHi-w/oSpace: paired dots with connecting lines ======  
    % Categories (same as before)
    % Categories (same as before)
    % Split goodcells into "significantly spatial" vs "not" so the two
    % groups can be visually contrasted (e.g. to check whether adding the
    % spatial predictor helps mainly cells that were already flagged as
    % significant, or also some that weren't).
    g_Sig   = goodcells & signicells;
    g_NS  = goodcells & ~signicells;
    
    % Plot
    figure; hold on;

    % Plot each category (same colors as before)
    % plot_cat draws each cell as a pair of dots (x=1 vs x=2, or x=3 vs
    % x=4) joined by a thin line, i.e. one line per cell showing how its
    % LLHi changes when the spatial predictor is added.
    plot_cat(1, LLHi_vis(g_Sig), 2, LLHi_pos(g_Sig),   [0 0.8 0], 0.9, 'p < 0.05');   % green-ish
    plot_cat(3, LLHi_vis(g_NS), 4, LLHi_pos(g_NS),   [0.5 0.5 0.5], 0.9, 'NS');
    
    % Cosmetics
    xlim([0.5 4.5]);
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'No space','Space', 'No space','Space'});
    ylabel('LLHi (bits / spike)');
    title('LLH improvement: No space vs Space');
    legend('Location','best'); legend boxoff;
    grid on; box off; set(gca,'TickDir','out');

end


%% EXAMPLE OF NATIVE PER-PREDICTOR KERNELS 
%(ONLY WORKS WITH SINGLE SESSION EXP STRUCTUREs)
if strcmp(figname,'GLM-kernels')
    % This figure shows the GLM's native per-predictor kernels, which are
    % only meaningful for a single recording session (the kernel
    % dimension structure collapses/changes shape once sessions are
    % concatenated into a batch EXP), hence this guard:
    if size(EXP.GLMs{1}.Tuning(itex).meanrespModel, 2) == 1
        error("`GLM-kernels` plots have to be called with an EXP structure from a single session, not an EXP structure obtained by concatenation of multiple sessions")
    end
    if isempty(singlecell_IDs) && nargin > 2
        error("invalid cell IDs")
    elseif isempty(singlecell_IDs)
        error("No cell Ids provided. Either provide cell indices or the full name ID of the cell (e.g. M26004_20260318_cell#281")
    end
    ncells = numel(singlecell_IDs);
    se_alpha = 0.1; % shading transparency for the SE bands (seplot)
    lwidth = 1.5;
    
    figure;
    
    % One row of 4 subplots per requested cell.
    for icell = 1:ncells
        cellidx = singlecell_IDs(icell);
        
        ax1 = subplot(ncells,4,4*(icell-1) + 1); hold on;

        % helper to flip left/right visual kernels with SE; x in [0,110]
        % All visual-space kernels (landmarks/BG/EOC/omission) are stored
        % over the full [0,110] deg visual-space axis, with 0 being lateral.
        % Since the standard is to consider 0 azimuth is front, we flip
        % them.
        get_vis_se = @(k,kf) deal( ...
            flipud(squeeze(EXP.GLMs{1}.Tuning(k).meanrespModel(cellidx,kf,:))), ...
            flipud(squeeze(EXP.GLMs{1}.Tuning(k).SErespModel(cellidx,kf,:))), ...
            linspace(0,110, size(EXP.GLMs{1}.Tuning(k).meanrespModel,3)) );

        % Landmarks (k=3): L1, L2
        [vL1,eL1,xv] = get_vis_se(itex,1);
        [vL2,eL2,~ ] = get_vis_se(itex,2);
        seplot(ax1, xv, vL1, eL1, se_alpha);
        plot(xv, vL1, 'Color', col_L1, 'LineWidth', lwidth, 'DisplayName','Landmark L1');
        seplot(ax1, xv, vL2, eL2, se_alpha);
        plot(xv, vL2, 'Color', col_L2, 'LineWidth', lwidth, 'DisplayName','Landmark L2');

        % Background textures (k=4): kf = 1..8 (colored gradient)
        % All nBG texture-segment kernels overlaid on the same axes; only
        % the first gets a legend entry (addedBG flag) to avoid a
        % cluttered legend with nBG duplicate "BG textures" labels.
        addedBG = false;
        for kf = 1:nBG
            [vbg,ebg,xbg] = get_vis_se(iBG,kf);
            seplot(ax1, xbg, vbg, ebg, se_alpha);
            h = plot(xbg, vbg, 'Color', bg_cols(kf,:), 'LineWidth', 1.0);
            if ~addedBG, set(h,'DisplayName','BG textures'); addedBG = true; else, set(h,'HandleVisibility','off'); end
        end

        % End of corridor (k=5): kf=1
        [vend,eend,xe] = get_vis_se(iEOC,1);
        seplot(ax1, xe, vend, eend, se_alpha);
        plot(xe, vend, 'Color', col_end, 'LineWidth', 1.5, 'DisplayName','End of corridor');

        % Skipped landmarks (Omit): L2, L3, L4 (dashed)
        % These are the kernels capturing the cell's response specifically
        % on trials where landmark L2, L3, or L4 was omitted (skipped) -
        % shown dashed to visually distinguish them from the "normal"
        % visual kernels above.
        [vSk1,eSk1,xs1] = get_vis_se(iOmit,1);
        [vSk2,eSk2,xs2] = get_vis_se(iOmit,2);
        [vSk3,eSk3,xs3] = get_vis_se(iOmit,3);
        seplot(ax1, xs1, vSk1, eSk1, se_alpha);
        plot(xs1, vSk1, '--', 'Color', col_omit2, 'LineWidth', lwidth, 'DisplayName','Skipped L2');
        seplot(ax1, xs2, vSk2, eSk2, se_alpha);
        plot(xs2, vSk2, '--', 'Color', col_omit3, 'LineWidth', lwidth, 'DisplayName','Skipped L3');
        seplot(ax1, xs3, vSk3, eSk3, se_alpha);
        plot(xs3, vSk3, '--', 'Color', col_omit4, 'LineWidth', lwidth, 'DisplayName','Skipped L4');

        % y=0 reference line: kernels live on the GLM's linear link-function
        % scale, so "no effect of this predictor" corresponds to 0 (unlike for
        % an log-link function for which the ref would be 1).
        hold on;yline(0, 'k--', 'HandleVisibility','off');

        xlabel('Visual space (deg)');
        ylabel(EXP.Spk.CellListString{cellidx}, 'Interpreter', 'none');
        title('Visual kernels (flipped left/right)');
        xlim([0 110]); box off; set(gca,'TickDir','out');
        legend('Location','northwest'); legend boxoff;
        pbaspect([16 9 1])


        % ---------- Subplot 2: VR onset/offset (k = 1) ----------
        ax2 = subplot(ncells,4,4*(icell-1) + 2); hold on;
        von   = squeeze(EXP.GLMs{1}.Tuning(iOnOff).meanrespModel(cellidx,1,:)); von = von(:);
        se_on = squeeze(EXP.GLMs{1}.Tuning(iOnOff).SErespModel(  cellidx,1,:)); se_on = se_on(:);
        voff   = squeeze(EXP.GLMs{1}.Tuning(iOnOff).meanrespModel(cellidx,2,:)); voff = voff(:);
        se_off = squeeze(EXP.GLMs{1}.Tuning(iOnOff).SErespModel(  cellidx,2,:)); se_off = se_off(:);
        t_on  = linspace(0,0.25,numel(von));
        t_off = linspace(0,0.25,numel(voff));
        seplot(ax2, t_on,  von,  se_on,  se_alpha);
        plot(t_on,  von,  'Color', col_vr_on,  'LineWidth',lwidth, 'DisplayName','Onset');
        seplot(ax2, t_off, voff, se_off, se_alpha);
        plot(t_off, voff, 'Color', col_vr_off, 'LineWidth',lwidth, 'DisplayName','Offset');
        
%         hold on;yline(1, 'k--', 'HandleVisibility','off');
        hold on;yline(0, 'k--', 'HandleVisibility','off');

        xlabel('Time (s)');
        title('VR onset/offset'); xlim([0 0.25]); box off; set(gca,'TickDir','out');
        legend('Location','northwest'); legend boxoff;
        pbaspect([16 9 1])


        % ---------- Subplot 3: Run speed (k = 2) ----------
        ax3 = subplot(ncells,4,4*(icell-1) + 3); hold on;
        vs   = squeeze(EXP.GLMs{1}.Tuning(iSpd).meanrespModel(cellidx,1,:)); vs = vs(:);
        se_s = squeeze(EXP.GLMs{1}.Tuning(iSpd).SErespModel(  cellidx,1,:)); se_s = se_s(:);
        xs = linspace(2,50,numel(vs));
        seplot(ax3, xs, vs, se_s, se_alpha);
        plot(xs, vs, 'Color', col_speed, 'LineWidth',lwidth);
        
%         hold on;yline(1, 'k--', 'HandleVisibility','off');
        hold on;yline(0, 'k--', 'HandleVisibility','off');

        xlabel('Run speed (cm s^{-1})');
        title('Run speed kernel'); xlim([2 50]); box off; set(gca,'TickDir','out');
        pbaspect([16 9 1])


        % ---------- Subplot 4: Spatial position kernel (k = iPos) ----------
        ax4 = subplot(ncells,4,4*(icell-1) + 4); hold on;
        vpos   = squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(cellidx,1,:)); vpos = vpos(:);
        se_pos = squeeze(EXP.GLMs{1}.Tuning(iPos).SErespModel(  cellidx,1,:)); se_pos = se_pos(:);
        xp = linspace(0,200,numel(vpos));
        seplot(ax4, xp, vpos, se_pos, se_alpha);
        plot(xp, vpos, 'k-', 'LineWidth',lwidth, 'DisplayName','Estimated');

        hold on;yline(0, 'k--', 'HandleVisibility','off');

        xlabel('Position (cm)');
        title('Spatial position kernel'); xlim([0 200]); box off; set(gca,'TickDir','out');
        pbaspect([16 9 1])
    end
end


%% TEXTURE LAYOUT IN CORRIDOR'S SPACE AND IN VISUAL SPACE FROM A SPECIFIC POSITION
if contains(figname,'Texture-layout')
    % This entire block is a static illustration of the task's fixed
    % visual layout, built purely from the geometry constants defined at
    % the top of the function - it does NOT read EXP at all.

    % viewer longitudinal position (% of corridor length)
    % Parsed from a trailing integer in figname, e.g. 'Texture-layout-50'
    % -> pos0 = 50; if no trailing number is found, default to 20%.
    pos0 = str2double(regexp(figname, '\d+$', 'match', 'once'));
    if isnan(pos0)
        pos0 = 20;
    end
    
    % subplot #1: BG textures as repeating colored patches over 0-200 cm, 
    % with landmarks overlay ---
    % Corridor and BG texture layout
    figure;
    subplot(1, 2, 1);  hold on;
    % BGseg_pct defines one repeating period's worth of segments (%
    % of corridor length); offset_pct steps through successive
    % periods (every period_pct%) to tile the pattern across the
    % full corridor.
    % Draw BG textures as contiguous patches across a 0-1 vertical band
    for offset_pct = 0:period_pct:100
        for kf = 1:size(BGseg_pct,1)
            s_pct = BGseg_pct(kf,1) + offset_pct;
            e_pct = BGseg_pct(kf,2) + offset_pct;
            if s_pct >= 100, continue; end
            e_pct = min(e_pct, 100);
            xs = s_pct * 2;
            xe = e_pct * 2;
            patch([xs xe xe xs], [0 0 1 1], bg_cols(kf,:), ...
                  'EdgeColor','none', 'FaceAlpha', 1.0);
        end
    end

    % Overlay landmarks as semi-transparent patches
    hw = width_pct/2; % half-width of a landmark patch, in % of corridor length
    % L1 at 20% and 60%
    for c = L1pos_pct
        xs = (c - hw) * 2; 
        xe = (c + hw) * 2;
        patch([xs xe xe xs], [0 0 1 1], col_L1, 'EdgeColor','none', ...
              'FaceAlpha', 1.0, 'DisplayName','L1');
    end
    % L2 at 40% and 80%
    for c = L2pos_pct
        xs = (c - hw) * 2; 
        xe = (c + hw) * 2;
        patch([xs xe xe xs], [0 0 1 1], col_L2, 'EdgeColor','none', ...
              'FaceAlpha', 1.0, 'DisplayName','L2');
    end

    % Cosmetics
    xlim([0 corridor_cm]); ylim([0 1]);
    yticks([]); yticklabels([]);
    xlabel('Position along corridor (cm)');
    title('Background textures and landmarks (0-200 cm)');
    box off; set(gca,'TickDir','out'); grid on;
    
    % --- Subplot 2: BG textures & landmarks as visual azimuth at viewer x=pos0% ---
    % Re-draws the same fixed-in-space layout, but projected through the
    % geometry of a viewer standing at corridor position pos0%, to show
    % what that viewer would actually see in visual-azimuth coordinates
    % (the coordinate system the GLM's visual kernels are defined on).
    subplot(1,2,2); hold on;
    % Geometry and mapping
    halfW_pct  = corrW_pct/2;
    % azmap: projects a corridor position (% of length) onto the visual
    % azimuth (deg) seen by a viewer standing at pos0%, using a simple
    % perspective (arctangent) projection across the corridor's width;
    % positions exactly at pos0% map to 90 deg (straight ahead).
    azmap = @(pct) 90 - atand((pct - pos0)/halfW_pct);   % pct -> azimuth (deg)

    % Repeat BG segments down the corridor and project to azimuth axis
    for offset_pct = 0:period_pct:100
        for kf = 1:size(BGseg_pct,1)
            s_pct = BGseg_pct(kf,1) + offset_pct;
            e_pct = BGseg_pct(kf,2) + offset_pct;
            if s_pct >= 100, continue; end
            e_pct = min(e_pct, 100);

            % Map endpoints to azimuth and clip to [0, 110]
            az1 = azmap(s_pct);
            az2 = azmap(e_pct);
            az_min = max(min(az1, az2), 0);
            az_max = min(max(az1, az2), 110);
            if az_max > az_min % skip segments that project entirely outside the visible field
                patch([az_min az_max az_max az_min], [0 0 1 1], bg_cols(kf,:), ...
                      'EdgeColor','none', 'FaceAlpha', 1.0);
            end
        end
    end

    % L1 patches
    for c = L1pos_pct
        s_pct = c - hw; e_pct = c + hw;
        az1 = azmap(s_pct); az2 = azmap(e_pct);
        az_min = max(min(az1, az2), 0);
        az_max = min(max(az1, az2), 110);
        if az_max > az_min
            patch([az_min az_max az_max az_min], [0 0 1 1], col_L1, ...
                  'EdgeColor','none', 'FaceAlpha', 1.0);
        end
    end

    % L2 patches
    for c = L2pos_pct
        s_pct = c - hw; e_pct = c + hw;
        az1 = azmap(s_pct); az2 = azmap(e_pct);
        az_min = max(min(az1, az2), 0);
        az_max = min(max(az1, az2), 110);
        if az_max > az_min
            patch([az_min az_max az_max az_min], [0 0 1 1], col_L2, ...
                  'EdgeColor','none', 'FaceAlpha', 1.0);
        end
    end

    % Axes cosmetics
    xlim([0 110]); ylim([0 1]);
    yticks([]); yticklabels([]);
    xlabel('Visual azimuth (deg)'); 
    title(sprintf('BG textures & landmarks in visual azimuth (viewer at x = %d %%)', pos0));
    grid on; box off; set(gca,'TickDir','out');
end


%% RESPONSE PROFILES ACROSS ALL CELLS, SPATIALLY SELECTIVE CELLS + SPATIAL KERNELS AND RESIDUALS
if contains(figname,'Resp-snake')
    figure;
    colormap(flipud(gray(256)));
    
    % Select which of the 6 conditions to display: if figname mentions
    % one or more condition names (e.g. 'Resp-snake-base-omit2'), show
    % only those; otherwise (the typical case, just 'Resp-snake') show
    % all 6.
    mask = false(size(condNames));
    for icond = 1:numel(condNames)
        if contains(figname, condNames{icond})
            mask(icond) = true;
        end
    end
    if ~any(mask)
        mask(:) = true;
    end
    
    cond2disp = find(mask);
    ncond = numel(cond2disp);
    % igroup loop: row-block 1 = all goodcells; row-block 2 = only the
    % spatially-selective subset of goodcells (signicells & goodLLHcells).
    for igroup = 1:2
        for icond = 1:ncond
            % Extract [cells x positions] matrix from EXP.Maps
            % resp_VS / resp_VSP are the GLM-predicted responses (Vision+
            % Speed only, and full Vision+Speed+Position respectively);
            % resp is the actual measured data, for this condition.
            resp_VS = EXP.Maps{mapVS_idx(cond2disp(icond))}.Tuning.meanrespModel(:,:);
            resp_VSP = EXP.Maps{mapVSP_idx(cond2disp(icond))}.Tuning.meanrespModel(:,:);
            resp = squeeze(EXP.Maps{mapDATA_idx(cond2disp(icond))}.Tuning.meanrespModel(:, 1, :));
            
            resp = resp(goodcells, :);  % keep good cells only
            resp_VS = resp_VS(goodcells,:);
            resp_VSP = resp_VSP(goodcells,:);

            % This is purely a display normalization (each cell's
            % own heatmap row is independently rescaled to [0,1]) so
            % cells with very different firing rates remain visually
            % comparable; it has no bearing on cell selection or
            % statistics.

            % Min-max normalize each cell (row) to [0,1], ignoring NaNs
            mn = min(resp, [], 2, 'omitnan');
            mx = max(resp, [], 2, 'omitnan');
            resp_norm = (resp - mn) ./ (mx - mn);
            
            mn = min(resp_VS, [], 2, 'omitnan');
            mx = max(resp_VS, [], 2, 'omitnan');
            resp_VS = (resp_VS - mn) ./ (mx - mn);
            
            mn = min(resp_VSP, [], 2, 'omitnan');
            mx = max(resp_VSP, [], 2, 'omitnan');
            resp_VSP = (resp_VSP - mn) ./ (mx - mn);
            
            % Order cells by peak position after normalization for base
            % condition. This ordering is computed ONCE (icond==1, i.e.
            % the first condition in cond2disp) and then reused for every
            % subsequent condition/row, so the same cell always appears
            % at the same row across panels - this is what makes the
            % heatmaps directly comparable column-to-column.
            if icond == 1
                [~, peakIdx] = max(resp_norm, [], 2, 'omitnan');
                [~, order] = sort(peakIdx, 'ascend', 'MissingPlacement', 'last');
            end
            resp_sorted = resp_norm(order, :);
            resp_VS = resp_VS(order,:);
            resp_VSP = resp_VSP(order,:);
            
            % Row-block 2: further restrict to the spatially-selective
            % subset (within goodcells), applied after the ordering above
            % so the subset still appears in the same relative order.
            if igroup == 2
                good_idx        = find(goodcells);
                ss_mask_good    = signicells(good_idx) & goodLLHcells(good_idx);  % spatially selective within goodcells
                ss_idx_sorted    = ss_mask_good(order); 
                resp_sorted  = resp_sorted(ss_idx_sorted, :);
                resp_VS = resp_VS(ss_idx_sorted,:);
                resp_VSP = resp_VSP(ss_idx_sorted,:);
            end
                        
            % This axis is shared by all three heatmap rows below
            % (DATA / VS / VSP) since they all span the same corridor.

            % X axis in cm (0-200)
            nBins = size(resp_sorted, 2);
            xbins = linspace(0, 200, nBins);
            nC = size(resp_sorted, 1);

            % Plot overall response for all cells 
            % column-block 1 of this igroup block: measured DATA, this condition.
            subplot(3,(ncond + 1) * 3,(igroup - 1) * (ncond + 1) * 3 + icond);
            imagesc(xbins, 1:nC, resp_sorted);
            set(gca, 'Clim', [0 1]);
            set(gca, 'YDir','normal');
            
            if icond == 1
                xlabel('Position along corridor (cm)');
                if igroup == 2
                    ylabel('Spatially selective cells (DATA)');
                else
                    ylabel('All cells (DATA)');
                end
            end
            if igroup == 1
                title(condNames{cond2disp(icond)});
            end
            box off; set(gca,'TickDir','out'); axis tight;
            
            % column-block 2: GLM-predicted response, Vision+Speed model (no
            % spatial predictor), this condition.
            subplot(3,(ncond + 1) * 3,(igroup - 1) * (ncond + 1) * 3 + icond + (ncond + 1));
            imagesc(xbins, 1:nC, resp_VS);
            set(gca, 'Clim', [0 1]);
            set(gca, 'YDir','normal');
            
            if icond == 1
                xlabel('Position along corridor (cm)');
                if igroup == 2
                    ylabel('Spatially selective cells (VS)');
                else
                    ylabel('All cells (VS)');
                end
            end
            if igroup == 1
                title(condNames{cond2disp(icond)});
            end
            box off; set(gca,'TickDir','out'); axis tight;
            
            % column-block 3: GLM-predicted response, full Vision+Speed+Position
            % model, this condition.
            subplot(3,(ncond + 1) * 3,(igroup - 1) * (ncond + 1) * 3 + icond + 2*(ncond + 1));
            imagesc(xbins, 1:nC, resp_VSP);
            set(gca, 'Clim', [0 1]);
            set(gca, 'YDir','normal');

            
            if icond == 1
                xlabel('Position along corridor (cm)');
                if igroup == 2
                    ylabel('Spatially selective cells (VSP)');
                else
                    ylabel('All cells (VSP)');
                end
            end
            if igroup == 1
                title(condNames{cond2disp(icond)});
            end
            box off; set(gca,'TickDir','out'); axis tight;
        end
    end
    

    % ---------- Spatial kernels and residuals ----------
    % Bottom block - first column block: the GLM's native spatial kernel (one value per
    % position bin, independent of trial condition) for the same cell
    % ordering used above, 
    % Bottom block - 2nd and 3rd column block: per-condition residual 
    % heatmaps showing where each model variant (VS or VSP) deviates from 
    % the actual data.
    spat_all    = squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(:, 1, :));  % [cells x posbins]
    spat_good   = spat_all(goodcells, :); % only goodcells
    spat_sorted = spat_good(order, :);
    
    good_idx        = find(goodcells);
    ss_mask_good    = signicells(good_idx) & goodLLHcells(good_idx);  % spatially selective within goodcells
    ss_idx_sorted    = ss_mask_good(order); 
    
    % Normalize kernels row-wise to [0,1]
    mn  = min(spat_sorted, [], 2, 'omitnan');
    mx  = max(spat_sorted, [], 2, 'omitnan');
    denom = mx - mn; denom(denom == 0) = NaN;
    spat_sorted_norm = (spat_sorted - mn) ./ denom;

    % xkern is reused below (shared x-axis) for all kernel/residual
    % heatmaps in this bottom block.

    xkern = linspace(0, 200, size(spat_sorted_norm,2));  % 0-200 cm axis% same row order as resp_sorted
    
    % ---------- Spatial kernels for spatially selective cells, same order as responses ----------
    axSp1 = subplot(3, (ncond + 1) * 3, 2 * (ncond + 1) * 3 + 1);
    imagesc(xkern, 1:sum(ss_idx_sorted), spat_sorted_norm(ss_idx_sorted, :), [0 1]);
    set(gca,'YDir','normal');
    xlabel('Position along corridor (cm)');
    ylabel('Spatially selective cells');
    colormap(axSp1, RedWhiteBlue);
    title('Spatial kernels (same order as responses)');
    box off; set(gca,'TickDir','out'); axis tight;
    
    for icond = 1:ncond
        VSP_all    = squeeze(EXP.Maps{mapVSP_idx(cond2disp(icond))}.Tuning.meanrespModel(:, 1, :));
        VSP_good   = VSP_all(goodcells, :);
        VSP_sorted = VSP_good(order, :);

        VS_all    = squeeze(EXP.Maps{mapVS_idx(cond2disp(icond))}.Tuning.meanrespModel(:, 1, :));
        VS_good   = VS_all(goodcells, :);
        VS_sorted = VS_good(order, :);

        DATA_all    = squeeze(EXP.Maps{mapDATA_idx(cond2disp(icond))}.Tuning.meanrespModel(:, 1, :));
        DATA_good   = DATA_all(goodcells, :);
        DATA_sorted = DATA_good(order, :);

        % denom: per-cell normalization factor (DATA's own dynamic
        % range), shared by both residual maps below, so DATA-VS and
        % DATA-VSP residuals are expressed on the same relative scale
        % per cell and remain comparable to each other.
        mn  = min(DATA_sorted(ss_idx_sorted, :), [], 2, 'omitnan');
        mx  = max(DATA_sorted(ss_idx_sorted, :), [], 2, 'omitnan');
        denom = mx - mn; denom(denom == 0) = NaN;

        % Residual map: DATA minus the Vision+Speed-only prediction. A
        % large residual here, in a cell that's spatially selective,
        % indicates the spatial predictor (absent from VS) is doing real
        % work explaining that part of the response.
        axSp1 = subplot(3, (ncond + 1) * 3, 2 * (ncond + 1) * 3 + (ncond + 1) + icond);
        VSres_norm = ((DATA_sorted(ss_idx_sorted, :) - VS_sorted(ss_idx_sorted, :))) ./ denom;

        imagesc(xkern, 1:sum(ss_idx_sorted), VSres_norm, [-1 1]);
        set(gca,'YDir','normal');
        xlabel('Position along corridor (cm)');
        ylabel('Spatially selective cells');
        colormap(axSp1, RedWhiteBlue);
        if icond == 1
            xlabel('Position along corridor (cm)');
            ylabel('Spatially selective cells (DATA - VS)');
        end
        box off; set(gca,'TickDir','out'); axis tight;
        

        % Residual map: DATA minus the full Vision+Speed+Position
        % prediction. Residuals here reflect what the full model still
        % fails to capture (model misfit / unexplained variance), since
        % unlike the VS residual above, the spatial predictor is already
        % included.
        axSp1 = subplot(3, (ncond + 1) * 3, 2 * (ncond + 1) * 3 + 2*(ncond + 1) + icond);
        VSPres_norm = ((DATA_sorted(ss_idx_sorted, :) - VSP_sorted(ss_idx_sorted, :))) ./ denom;

        imagesc(xkern, 1:sum(ss_idx_sorted), VSPres_norm, [-1 1]);
        set(gca,'YDir','normal');
        if icond == 1
            xlabel('Position along corridor (cm)');
            ylabel('Spatially selective cells (DATA - VSP)');
        end
        colormap(axSp1, RedWhiteBlue);
        box off; set(gca,'TickDir','out'); axis tight;
    end
end


%% EXAMPLE OF CELL RESPONSE WITH ACTUAL DATA AND FULL VSP PREDICTION OVERLAYED
if strcmp(figname,'Resp-singleCell')
    cellidx = singlecell_IDs;
    
    nCond = numel(condNames);
    ncells = numel(singlecell_IDs);

    figure;
    for icell = 1:ncells
        cellidx = singlecell_IDs(icell);
        % Pre-scan to set a common Y limit (start at 0). Looping once
        % over all 6 conditions first (collecting y/se/yGLM/x per
        % condition) lets every condition's subplot for this cell share
        % the same y-axis range, computed from the global max across all
        % of them, before any subplot is actually drawn.
        yMax = 0;
        yAll  = cell(1,nCond);
        seAll = cell(1,nCond);
        yGLM  = cell(1,nCond);
        xAll  = cell(1,nCond);

        for k = 1:nCond
            % y/se: measured data mean response +/- SE, this condition
            % (EXP.Maps mapDATA_idx group)
            y  = squeeze(EXP.Maps{mapDATA_idx(k)}.Tuning.meanrespModel(cellidx,1,:));   y  = y(:);
            se = squeeze(EXP.Maps{k}.Tuning.SErespModel(cellidx,1,:)); se = se(:);
            % y_g: full-model (VSP) GLM-predicted mean response, same condition
            y_g = squeeze(EXP.Maps{mapVSP_idx(k)}.Tuning.meanrespModel(cellidx,1,:)); y_g = y_g(:);

            x = linspace(0,200,numel(y));

            yAll{k}  = y;
            seAll{k} = se;
            yGLM{k}  = y_g;
            xAll{k}  = x;

            yMax = max([yMax; y + se; y_g], [], 'omitnan');
        end
        if ~isfinite(yMax) || yMax <= 0, yMax = 1; end
        yLim = [0, 1.05*yMax];

        % Colors
        col_data = [0 0 0];        % black for measured mean
        col_glm  = col_Spatial; % blue-ish for GLM prediction

        % Now draw the actual subplots, one per condition, reusing the
        % shared yLim computed above.
        for k = 1:nCond
            ax = subplot(ncells,nCond,nCond * (icell-1) + k); hold(ax,'on');

            % SE shading
            seplot(ax, xAll{k}, yAll{k}, seAll{k}, 0.25);

            % Mean response (data)
            plot(ax, xAll{k}, yAll{k}, '-', 'Color', col_data, 'LineWidth', 1.5, 'DisplayName','Data');

            % GLM predicted mean
            plot(ax, xAll{k}, yGLM{k}, '-', 'Color', col_glm, 'LineWidth', 1.5, 'DisplayName','GLM');

            % Axes / labels
            xlim(ax,[0 200]); ylim(ax,yLim);
            box(ax,'off'); set(ax,'TickDir','out');
            if icell == 1
                title(ax, sprintf('%s', condNames{k}));
            end
            if k == 1
                ylabel(ax, EXP.Spk.CellListString{cellidx}, 'Interpreter', 'none');
            end
            if icell == ncells
                xlabel(ax, 'Position (cm)');
            end

            % Legend on the last subplot only
            if k == nCond && icell == 1
                legend(ax, 'Location','northwest'); legend(ax,'boxoff');
            end
            pbaspect([16 9 1])
        end
    end

    % Overall title
    sgtitle(sprintf('Cell %d - mean responses with SE and GLM overlay', cellidx));
end


%% METRICS SUMMARY RELATING SPATIAL COMPONENT STRENGTH TO RF AND RESPONSE PROFILE PROPERTIES
if strcmp(figname,'RF-LLHrel-summary')
   
    % Masks
    % spatial_mask / nonspatial_mask partition goodcells into spatially
    % selective vs. not, used throughout this figure to contrast the two
    % populations (overlaid histograms, etc.).
    spatial_mask   = goodcells & signicells & goodLLHcells;     % spatially selective
    nonspatial_mask= goodcells & ~(signicells & goodLLHcells);    % non-spatial
    LLHplot = LLHrel;% LLHrel_omit;% the metric plotted throughout this
                      % figure; swap the active line to LLHrel_omit to
                      % instead summarize the omission-component metric

    % Spatial kernel (iPos)
    spatK = squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(:,1,:));   % [cells x posbins]
    % peak_pos: position bin where each cell's spatial GLM kernel
    % peaks; overall_peak (below) is the analogous peak position of
    % the cell's raw mean response (base condition) instead of the
    % fitted kernel - comparing the two checks whether the kernel and
    % the actual response agree on where the cell 'prefers'.

    % Peak positions (in cm, 0-200) for spatial cells
    xpos = 0:2:200;
    [~, peak_pos] = max(spatK,[],2,'omitnan');
    peak_pos(~isnan(peak_pos)) = xpos(peak_pos(~isnan(peak_pos)));
    %mean resp peak position
    overall_resp = EXP.Maps{mapDATA_idx(1)}.Tuning.meanrespModel;
    [~, overall_peak] = max(overall_resp(:,:),[],2,'omitnan');
    overall_peak(~isnan(overall_peak)) = xpos(overall_peak(~isnan(overall_peak)));

    % RF positions (visual azimuth, degrees)
    RFpos = EXP.GLMs{1}.Perf.bestRFpos(:);

    % Visual latencies (seconds)
    Vislat = EXP.GLMs{1}.Perf.bestVisdelay(:);

    % Bin definitions
    % Bin edges used by the histograms/violins/bar-plots below. RF_edges
    % is deliberately narrow (32.5-47.5 deg) rather than spanning the
    % full [0,110] visual-azimuth range, reflecting where this dataset's
    % RF positions actually concentrate (see commented-out wider range).
    % LLHrel bins
    LLH_edges = 0.9:0.01:1.5;
    % RF position bins (visual azimuth 0-110 deg, as used elsewhere)
    RF_edges = 32.5:5:47.5;% 0.5:22.5;% 
    RF_centers = (RF_edges(1:end-1)+RF_edges(2:end))/2;
    % Peak pos edges and bin size
    peak_edges = 0:2:200;% 0.5:22.5;% 
    peak_binwidth = 20;
    % Latency bins (0-0.25 s default)
    lat_min = -1000/60; lat_max = 500+1000/60;

    LAT_edges = lat_min:2000/60:lat_max; 
    LAT_centers = (LAT_edges(1:end-1)+LAT_edges(2:end))/2;

    % Figure layout: 3 rows x 6 columns of tiles; several panels below
    % span multiple tile columns ([1 2]) for a wider plot.
    figure;
    tl = tiledlayout(3,6, 'TileSpacing','compact','Padding','compact');

    % ===== (1) Dist of spatial-kernel peak positions =====
    % Compares where the cell's peak is according to two different
    % measures: the raw trial-averaged response (overall_peak) vs. the
    % fitted GLM spatial kernel (peak_pos). Large disagreement between
    % the two would suggest the kernel fit doesn't track the data well.
    ax1 = nexttile(tl); hold(ax1,'on');
    pk = overall_peak(spatial_mask);
    histogram(ax1, pk, 'BinWidth', 10, 'FaceColor',[0.3 0.3 0.3],'EdgeColor','none', 'DisplayName','Mean response peak');
    pk = peak_pos(spatial_mask);
    histogram(ax1, pk, 'BinWidth', 10, 'FaceColor',col_Spatial,'EdgeColor','none', 'DisplayName','Spatial kernel peak');
    xlim(ax1,[0 200]);
    xlabel(ax1,'Peak position (cm)'); ylabel(ax1,'# cells');
    title(ax1,'Spatial-kernel peak positions');
    box(ax1,'off'); set(ax1,'TickDir','out');
    
    % ===== (2) Scatter: spatial-kernel peak vs overall-response peak =====
    ax2 = nexttile(tl); hold(ax2,'on');
    mask_sc = spatial_mask;
    scatter(ax2, overall_peak(mask_sc), peak_pos(mask_sc), ptSize, [0 0 0], 'filled', 'MarkerFaceAlpha', ptAlpha);
    % identity line
    plot(ax2, [0 200],[0 200], 'r-', 'LineWidth',1, 'HandleVisibility','off');
    xlim(ax2,[0 200]); ylim(ax2,[0 200]);
    xlabel(ax2,'Overall response peak (cm)'); ylabel(ax2,'Spatial-kernel peak (cm)');
    title(ax2,'Kernel peak vs response peak');
    box(ax2,'off'); set(ax2,'TickDir','out');
    pbaspect([1, 1, 1])
    
    % ===== (3) LLHrel vs peak position (violin) =====
    % Two violin plots back to back: LLHrel binned by overall-response
    % peak, then LLHrel binned by spatial-kernel peak, so it's easy to
    % see whether spatial coding strength varies systematically across
    % the corridor (e.g. weaker near the start/end) for either measure.
    ax3 = nexttile(tl); hold(ax3,'on');
    msk = spatial_mask;%goodcells;
    plotViolin(ax3, overall_peak(msk), LLHplot(msk), peak_edges, peak_binwidth);
    xlim(ax3,[peak_edges(1) peak_edges(end)]);
    xlabel(ax3,'Overall peak position (cm)'); ylabel(ax3,'LLH_{rel}');
    title(ax3,'LLH_{rel} vs Overall response peak');
    legend(ax3,'Location','northwest'); legend(ax3,'boxoff');
    box(ax3,'off'); set(ax3,'TickDir','out');
    
    ax3 = nexttile(tl); hold(ax3,'on');
    msk = spatial_mask;
    plotViolin(ax3, peak_pos(msk), LLHplot(msk), peak_edges, peak_binwidth);
    xlim(ax3,[peak_edges(1) peak_edges(end)]);
    xlabel(ax3,'Space peak position (cm)'); ylabel(ax3,'LLH_{rel}');
    title(ax3,'LLH_{rel} vs Spatial kernel Peak');
    legend(ax3,'Location','northwest'); legend(ax3,'boxoff');
    box(ax3,'off'); set(ax3,'TickDir','out');

    % ===== (4) Dist of LLHrel - overlayed non-spatial vs spatial =====
    ax3 = nexttile(tl, [1 2]); hold(ax3,'on');
    h_ns = histogram(ax3, LLHplot(goodcells), LLH_edges, ...
        'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none', 'FaceAlpha',0.6, 'DisplayName',sprintf('Non-spatial (n=%d)',sum(nonspatial_mask)));
    h_sp = histogram(ax3, LLHplot(spatial_mask),    LLH_edges, ...
        'FaceColor',[0.2 0.2 0.2], 'EdgeColor','none', 'FaceColor',col_Spatial, 'FaceAlpha',1.0, 'DisplayName',sprintf('Spatial (n=%d)',sum(spatial_mask)));
    xlabel(ax3,'LLH_{rel}'); ylabel(ax3,'# cells');
    title(ax3,'LLH_{rel} distribution (overlay)');
    legend(ax3,'Location','northwest'); legend(ax3,'boxoff');
    box(ax3,'off'); set(ax3,'TickDir','out');

    % ===== (5) Fraction spatially selective vs RF position =====
    % N_all==0 -> NaN guards against a 0/0 division producing a
    % misleading 0% bar for RF-position bins with no goodcells at all
    % (NaN bars are simply omitted by bar(), rather than drawn as 0%).
    ax3 = nexttile(tl, [1 2]); hold(ax3,'on');
    [N_all,~] = histcounts(RFpos(goodcells), RF_edges);
    [N_spa,~] = histcounts(RFpos(spatial_mask), RF_edges);
    N_all(N_all == 0) = NaN;
    frac_spa = N_spa ./ N_all;
    bar(ax3, RF_centers, frac_spa, 'stacked', 'BarWidth', 1, 'EdgeColor', 'none');
    ylim(ax3,[0 1]); xlim(ax3,[RF_edges(1) RF_edges(end)]);
    xlabel(ax3,'RF position (deg)');
    ylabel(ax3,'Fraction spatial');
    title(ax3,'Spatial fraction vs RF position');
    box(ax3,'off'); set(ax3,'TickDir','out');

    % ===== (6) LLHrel vs RF position (violin) =====
    ax4 = nexttile(tl, [1 2]); hold(ax4,'on');
    msk = goodcells;
    plotViolin(ax4, RFpos(msk), LLHplot(msk), RF_edges);
    xlim(ax4,[RF_edges(1) RF_edges(end)]);
    xlabel(ax4,'RF position (deg)'); ylabel(ax4,'LLH_{rel}');
    title(ax4,'LLH_{rel} vs RF position');
    legend(ax4,'Location','northwest'); legend(ax4,'boxoff');
    box(ax4,'off'); set(ax4,'TickDir','out');

    % ===== (7) Dist of RF positions - overlayed non-spatial vs spatial =====
    ax5 = nexttile(tl, [1 2]); hold(ax5,'on');
    histogram(ax5, RFpos(goodcells), RF_edges, 'FaceColor',[0.3 0.3 0.3], 'EdgeColor','none', 'DisplayName',sprintf('Non-spatial (n=%d)',sum(nonspatial_mask)));
    histogram(ax5, RFpos(spatial_mask), RF_edges, 'EdgeColor','none', 'FaceColor',col_Spatial, 'FaceAlpha',1.0, 'DisplayName',sprintf('Spatial (n=%d)',sum(spatial_mask)));
    xlim(ax5,[RF_edges(1) RF_edges(end)]);
    xlabel(ax5,'RF position (deg)'); ylabel(ax5,'# cells');
    title(ax5,'RF position distribution');
    legend(ax5,'Location','northwest'); legend(ax5,'boxoff');
    box(ax5,'off'); set(ax5,'TickDir','out');

    % ===== (8) Fraction spatially selective vs visual latency =====
    ax6 = nexttile(tl, [1 2]); hold(ax6,'on');
    [N_allL,~] = histcounts(Vislat(goodcells), LAT_edges);
    [N_spaL,~] = histcounts(Vislat(spatial_mask & goodcells), LAT_edges);
    frac_spa_L = N_spaL ./ max(N_allL,1);
    bar(ax6, LAT_centers, frac_spa_L, 'stacked', 'BarWidth', 1, 'EdgeColor', 'none');
    xlim(ax6,[LAT_edges(1) LAT_edges(end)]); ylim(ax6,[0 1]);
    xlabel(ax6,'Visual latency (s)');
    ylabel(ax6,'Fraction spatial');
    title(ax6,'Spatial fraction vs RF latency');
    box(ax6,'off'); set(ax6,'TickDir','out');

    % ===== (9) LLHrel vs visual latency (violin) =====
    ax7 = nexttile(tl, [1 2]); hold(ax7,'on');
    mskL = goodcells;
    plotViolin(ax7, Vislat(mskL), LLHplot(mskL), LAT_edges);
    xlim(ax7,[LAT_edges(1) LAT_edges(end)]);
    xlabel(ax7,'Visual latency (s)'); ylabel(ax7,'LLH_{rel}');
    title(ax7,'LLH_{rel} vs RF latency');
    legend(ax7,'Location','northwest'); legend(ax7,'boxoff');
    box(ax7,'off'); set(ax7,'TickDir','out');

    % ===== (10) Dist of visual latency =====
    ax8 = nexttile(tl, [1 2]); hold(ax8,'on');
    histogram(ax8, Vislat(goodcells), LAT_edges, 'FaceColor',[0.3 0.3 0.3], 'EdgeColor','none', 'DisplayName',sprintf('Non spatial (n=%d)',sum(nonspatial_mask)));
    histogram(ax8, Vislat(spatial_mask), LAT_edges, 'EdgeColor','none', 'FaceColor',col_Spatial, 'FaceAlpha',1.0, 'DisplayName',sprintf('Spatial (n=%d)',sum(spatial_mask)));
    xlim(ax8,[LAT_edges(1) LAT_edges(end)]);
    xlabel(ax8,'Visual latency (s)'); ylabel(ax8,'# cells');
    title(ax8,'Latency distribution');
    legend(ax8,'Location','northwest'); legend(ax8,'boxoff');
    box(ax8,'off'); set(ax8,'TickDir','out');

    % Overall title
    sgtitle(tl, 'Spatial kernel peaks, LLH_{rel}, RF position & latency summaries');
end


%% 
if strcmp(figname,'Cluster-SpaceKernels')
    % spaceK: GLM spatial kernel per cell (the thing being clustered)
    % resp: raw mean response per cell (base condition), used only for
    % a side-by-side visual comparison panel further below, not in the
    % clustering itself.
    spaceK = squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(:,1,:));
    resp  = squeeze(EXP.Maps{mapDATA_idx(1)}.Tuning.meanrespModel(:,1,:));
    x_pos = linspace(0,200,size(spaceK,2));
    % Combine series + animal into a single unique session identifier
    % (only meaningful when EXP.Spk.animal exists, i.e. data spans
    % multiple animals using overlapping series numbers).
    sessionid = EXP.Spk.series;
    if isfield(EXP.Spk, 'animal')
        sessionid = sessionid + 10*EXP.Spk.animal;
    end

    % valid "good" cells
    % Restrict to the canonical spatially-selective population before
    % clustering - clustering non-spatial cells' kernels would just be
    % clustering noise.
    spatial_mask = goodcells & signicells & goodLLHcells;
    spaceK = spaceK(spatial_mask,:);
    resp  = resp(spatial_mask,:);
    sess_id = sessionid(spatial_mask);
    [sess_vals, ~, sess_uniq] = unique(sess_id);

    % --- Cluster ---
    % See the ClusterSpaceKernels local function (end of file) for the
    % full clustering pipeline: z-score, optional smoothing, optional
    % PCA, hierarchical linkage for the dendrogram, then k-means for the
    % final cluster assignment T.
    [T, Z, leafOrder, colorCutoff] = ClusterSpaceKernels(spaceK, n_Clusters, clust_nPCs, clust_smthwin, dist_metric, dist_method);

    % --- Smoothing window; 0 = no smoothing ---
%     smth_win = 0;
    if clust_smthwin > 0 && exist('SpecialSmooth','file')==2
        spaceK_smth = SpecialSmooth(spaceK, [0 clust_smthwin/100], [1, 100]);
    else
        spaceK_smth = spaceK;
    end

    % --- Order clusters by peak position of their (smoothed) mean kernel ---
    % cluster_ids: cluster labels in the order they first appear along
    % the dendrogram's leaf order, so clusters are drawn left-to-right in
    % roughly the same order the dendrogram visually groups them.
    cluster_ids = unique(T(leafOrder), 'stable');
    cluster_ids = cluster_ids(:)';

    % --- Precompute per-cluster scaled matrices + mean/SE for plotting ---
    Nmax_cells = sum(T == mode(T)); % size of the largest cluster, used to
                                     % give all per-cluster imagesc panels
                                     % the same y-axis range
    cmap = lines(numel(cluster_ids));

    % row-wise scaler: min-max normalize each row (cell) to [0,1], purely
    % for display so cells with different firing-rate scales remain
    % visually comparable within and across cluster panels.
    minmax_scale = @(M) ( ...
        (M - min(M,[],2,'omitnan')) ./ max(min( max(M,[],2,'omitnan') - min(M,[],2,'omitnan'), [], 2), eps) );

    spaceK_smth = minmax_scale(spaceK_smth);
    spaceK      = minmax_scale(spaceK);
    resp_scaled     = minmax_scale(resp);

    % ===================== PLOTTING =====================
    figure;
    colormap(flipud(gray(256)));
    ncol = 6;
    
    % Use a single, consistent grid: 5 rows x Nclusters columns
    % Col 1: dendrogram spanning all columns
    % Col 2: per-cluster mean (+/-SE) of smoothed, scaled kernels
    % Col 3: per-cluster imagesc of smoothed, scaled kernels
    % Col 4: per-cluster imagesc of raw, scaled kernels
    % Col 5: per-cluster imagesc of response, scaled
    % Col 6: distribution of each cluster across recording sessions

    % ----- Row 1: dendrogram spanning all columns (consistent subplot grid) -----
    subplot(numel(cluster_ids), ncol, 1:ncol:ncol*numel(cluster_ids));
%     dendrogram(Z, size(spaceK,1), ...
%         'Reorder', leafOrder, ...
%         'Labels',  num2str(T), ...
%         'ColorThreshold', colorCutoff, ...
%         'CheckCrossing', false,...
%         'Orientation','left');
    dendrogram(Z, 4*n_Clusters, ...
        'Reorder', leafOrder, ...
        'Labels',  num2str(T), ...
        'ColorThreshold', colorCutoff, ...
        'CheckCrossing', false,...
        'Orientation','left');
    set(gca,'TickDir','out'); box off;
    title('Dendrogram (labels = cluster IDs)');

    % Loop clusters in the ordered list
    for p = 1:numel(cluster_ids)
        c = cluster_ids(p);
        col_c = cmap(p,:);
        rows = (T==c); % logical mask selecting this cluster's cells

        % ----- Row 2: mean (+/-SE) of smoothed, scaled kernels -----
        ax = subplot(numel(cluster_ids), ncol, (p - 1) * ncol + 2); hold(ax,'on');
        M = spaceK(rows,:);           % [n_c x posbins]
        mu = mean(M,1,'omitnan');
        se = std(M,0,1,'omitnan');% / sum(rows).^0.5;
        % SE (gray)
        seplot(ax, x_pos, mu, se, 0.2);
        % mean line (cluster color)
        plot(ax, x_pos, mu, 'Color', col_c, 'LineWidth', 1.6);
        ylim(ax,[0 1]); xlim(ax,[0 200]);
        title(ax, sprintf('C%d mean (n=%d)', c, sum(rows)));
        if p==1, ylabel(ax,'Scaled'); end
        set(ax,'TickDir','out'); box(ax,'off');

        % ----- Row 3: smoothed, scaled kernels (imagesc) -----
        % Same data as row 4 below but with clust_smthwin smoothing
        % applied (if enabled); compare the two rows to judge how much
        % smoothing changes the apparent kernel shape.
        ax = subplot(numel(cluster_ids),ncol, (p - 1) * ncol + 3);
        Kc = spaceK_smth(rows,:);
        imagesc(ax, x_pos, 1:size(Kc,1), Kc, [0 1]);
        ylim(ax,[0.5 Nmax_cells+0.5]); xlim(ax,[0 200]);
        if p==1, ylabel(ax,'Smoothed kernels'); end
        title(ax, sprintf('C%d', c));
        set(ax,'TickDir','out'); box(ax,'off');

        % ----- Row 4: raw, scaled kernels (imagesc) -----
        ax = subplot(numel(cluster_ids), ncol, (p - 1) * ncol + 4);
        Kc = spaceK(rows,:);
        imagesc(ax, x_pos, 1:size(Kc,1), Kc, [0 1]);
        ylim(ax,[0.5 Nmax_cells+0.5]); xlim(ax,[0 200]);
        if p==1, ylabel(ax,'Raw kernels'); end
        set(ax,'TickDir','out'); box(ax,'off');

        % ----- Row 5: responses, scaled (imagesc) -----
        % The actual mean response (base condition), for the same cells
        % in the same row order as rows 3/4, to visually cross-check that
        % each cluster's GLM kernel shape is reflected in the real data.
        ax = subplot(numel(cluster_ids), ncol, (p - 1) * ncol + 5);
        Rc = resp_scaled(rows,:);
        imagesc(ax, x_pos, 1:size(Rc,1), Rc, [0 1]);
        ylim(ax,[0.5 Nmax_cells+0.5]); xlim(ax,[0 200]);
        if p==1, ylabel(ax,'Responses'); end
        xlabel(ax,'Position (cm)');
        set(ax,'TickDir','out'); box(ax,'off');
        
        % ----- Row 6: distribution across sessions or simulated shapes -----
        % Sanity check: if a cluster is dominated by cells from a single
        % recording session (or animal, via the sessionid offset above),
        % that cluster might reflect a session-specific artifact rather
        % than a genuine, reproducible response category.
        ax = subplot(numel(cluster_ids), ncol, (p - 1) * ncol + 6);
        Sc = sess_uniq(rows,:);
        histogram(Sc, 0.5:1:numel(sess_vals)+0.5, 'EdgeColor','none', 'FaceColor',col_Spatial, 'FaceAlpha',1.0);
        set(gca, 'XTick', 1:numel(sess_vals), 'XTickLabel', sess_vals);
        if p==1, ylabel(ax,'# cells'); end
        xlabel(ax,'Session id');
        set(ax,'TickDir','out'); box(ax,'off');
    end
end


%% KERNELS AND MAPS ORDERED BY SIMILARITY
if contains(figname,'Ordered-Kernels')
    % Parse figname into a list of hyphen-separated tags, dropping the
    % leading 'Ordered' and 'Kernels' tokens. S(1) (the FIRST remaining
    % tag) is always treated as the reference kernel/map used for cell
    % ordering/grouping below; every tag matching a known kernel/map type
    % (handled in the X.* blocks further down) becomes an extra heatmap
    % column, in addition to whichever modifier tags (all/omission/hist/
    % peak) are also present anywhere in the list.
    S = strsplit(figname, '-');
    S = S(3:end);
    
    if any(contains(S, 'all'))
        show_all = true;
    else
        show_all = false;
    end
    if any(contains(S, 'omission'))
        show_omit = true;
    else
        show_omit = false;
    end
    
    if any(contains(S, 'hist'))
        show_hist = true;
    else
        show_hist = false;
    end
    if any(contains(S, 'peak'))
        orderby = 'peak';
    else
        orderby = 'similarity';
    end
    ref = S(1); % reference tag used for ordering/grouping (see below)
    
    % X: struct of fieldname -> {maps, scale, clims, colormap}, one entry
    % per requested kernel/map tag found in S. Built incrementally below;
    % only tags actually present in S get an entry, so X may end up with
    % anywhere from 1 to many fields depending on figname.
    X = [];
    if any(strcmp(S, 'Space'))
        X.Space.maps = squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(:,1,:));
        X.Space.scale = 'minmax';
        X.Space.clims = [0 1]; 
        X.Space.colormap = RedWhiteBlue;
    end
    if any(strcmp(S, 'Landmarks'))
        X.Landmarks.maps = squeeze(EXP.GLMs{1}.Tuning(itex).meanrespModel(:,1,:,1));
        X.Landmarks.scale = 'minmax';
        X.Landmarks.clims = [0 1];
        X.Landmarks.colormap = RedWhiteBlue;
    end
    if any(strcmp(S, 'BG'))
        X.BG.maps = squeeze(EXP.GLMs{1}.Tuning(iBG).meanrespModel(:,1,:,1));
        X.BG.scale = 'minmax';
        X.BG.clims = [0 1];
        X.BG.colormap = RedWhiteBlue;
    end
    if any(strcmp(S, 'EOC'))
        X.Tails.maps = squeeze(EXP.GLMs{1}.Tuning(iEOC).meanrespModel(:,1,:,1));
        X.Tails.scale = 'minmax';
        X.Tails.clims = [0 1];
        X.Tails.colormap = RedWhiteBlue;
    end
    if any(strcmp(S, 'Vis'))
        % 'Vis' = combined visual drive: landmark kernel + BG-texture
        % kernel summed together (both feature 1 only), as a single
        % proxy for "how much does vision alone drive this cell".
        X.Vis.maps = squeeze(EXP.GLMs{1}.Tuning(itex).meanrespModel(:,1,:,1))...
            + squeeze(EXP.GLMs{1}.Tuning(iBG).meanrespModel(:,1,:,1));
        X.Vis.scale = 'minmax';
        X.Vis.clims = [0 1];
        X.Vis.colormap = RedWhiteBlue;
    end
    if any(strcmp(S, 'FullMap'))
        X.Full.maps = squeeze(EXP.Maps{mapVSP_idx(1)}.Tuning.meanrespModel(:,1,:));
        X.Full.scale = 'minmax';
        X.Full.clims = [0 1];
        X.Full.colormap = flipud(gray(256));%RedWhiteBlue;
    end
    
    % Condition IDs for base, swap and Omit trials
    % Indices into condNames (1-6), used below to pick out specific
    % conditions from the various EXP.Maps / Tuning feature
    % dimensions for the diff/comparison tags (Omit, OmitMap,
    % DataOmitMapDiff, etc.).
    baseidx = find(strcmp('base', condNames));
    swap23idx = find(strcmp('swap23', condNames));
    swap34idx = find(strcmp('swap34', condNames));
    omit2idx = find(strcmp('omit2', condNames));
    omit3idx = find(strcmp('omit3', condNames));
    omit4idx = find(strcmp('omit4', condNames));
    
    if any(strcmp(S, 'Omit'))
        % GLM native omission KERNELS (Tuning(iOmit), feature dim indexed
        % by condition: omit2idx/3/4 here are used to index the FEATURE
        % dimension, not a separate trial-condition selection - i.e. each
        % of the 3 omission kernels is the cell's fitted response
        % specifically to that landmark being skipped, on its own native
        % visual-space axis. Contrast with 'OmitMap' below, which uses
        % the actual reconstructed MAPS for the omission conditions
        % instead, on the corridor-position axis.
        X.Omit2.maps = squeeze(EXP.GLMs{1}.Tuning(iOmit).meanrespModel(:,1,:,omit2idx));
        X.Omit2.scale = 'minmax';
        X.Omit2.clims = [0 1];
        X.Omit2.colormap = RedWhiteBlue;
        
        X.Omit3.maps = squeeze(EXP.GLMs{1}.Tuning(iOmit).meanrespModel(:,1,:,omit3idx));
        X.Omit3.scale = 'minmax';
        X.Omit3.clims = [0 1];
        X.Omit3.colormap = RedWhiteBlue;
        
        X.Omit4.maps = squeeze(EXP.GLMs{1}.Tuning(iOmit).meanrespModel(:,1,:,omit4idx));
        X.Omit4.scale = 'minmax';
        X.Omit4.clims = [0 1];
        X.Omit4.colormap = RedWhiteBlue;
    end
    
    if any(strcmp(S, 'OmitMap'))
        % Measured mean response (EXP.Maps mapDATA_idx group) for each
        % omission TRIAL CONDITION, i.e. the actual recorded data on
        % trials where that landmark was skipped - distinct from the
        % 'Omit' kernels above.
        X.omit2map.maps = squeeze(EXP.Maps{omit2idx}.Tuning.meanrespModel(:,1,:));
        X.omit2map.scale = 'minmax';
        X.omit2map.clims = [0 1];
        X.omit2map.colormap = flipud(gray(256));
        
        X.omit3map.maps = squeeze(EXP.Maps{omit3idx}.Tuning.meanrespModel(:,1,:));
        X.omit3map.scale = 'minmax';
        X.omit3map.clims = [0 1];
        X.omit3map.colormap = flipud(gray(256));
        
        X.omit4map.maps = squeeze(EXP.Maps{omit4idx}.Tuning.meanrespModel(:,1,:));
        X.omit4map.scale = 'minmax';
        X.omit4map.clims = [0 1];
        X.omit4map.colormap = flipud(gray(256));
    end
    
    if any(strcmp(S, 'DataOmitMapDiff'))
        % Each omission-condition diff = (omission trials) - (base
        % trials), computed on the MEASURED DATA. A positive diff right
        % around where the skipped landmark would have appeared suggests
        % the cell's responds positively to the background being removed.
        base_maps = squeeze(EXP.Maps{mapDATA_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit2_maps = squeeze(EXP.Maps{mapDATA_idx(omit2idx)}.Tuning.meanrespModel(:,1,:));
        X.omit2diff.maps = omit2_maps - base_maps;
        X.omit2diff.scale = 'minmax';
        X.omit2diff.clims = [0 1];
        X.omit2diff.colormap = RedWhiteBlue;
        
        base_maps = squeeze(EXP.Maps{mapDATA_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit3_maps = squeeze(EXP.Maps{mapDATA_idx(omit3idx)}.Tuning.meanrespModel(:,1,:));
        X.omit3diff.maps = omit3_maps - base_maps;
        X.omit3diff.scale = 'minmax';
        X.omit3diff.clims = [0 1];
        X.omit3diff.colormap = RedWhiteBlue;
        
        base_maps = squeeze(EXP.Maps{mapDATA_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit4_maps = squeeze(EXP.Maps{mapDATA_idx(omit4idx)}.Tuning.meanrespModel(:,1,:));
        X.omit4diff.maps = omit4_maps - base_maps;
        X.omit4diff.scale = 'minmax';
        X.omit4diff.clims = [0 1];
        X.omit4diff.colormap = RedWhiteBlue;
    end
    if any(strcmp(S, 'VisOmitMapDiff'))
        % Same omission-minus-base diff as above, but computed on the
        % Vision+Speed-only model PREDICTION instead of the data: shows
        % what the model would predict for the omission effect using
        % only visual-drive + speed, with no spatial/Position predictor.
        % Comparing this to DataOmitMapDiff tests whether the cell's
        % real omission response goes beyond what vision alone predicts.
        base_maps = squeeze(EXP.Maps{mapVS_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit2_maps = squeeze(EXP.Maps{mapVS_idx(omit2idx)}.Tuning.meanrespModel(:,1,:));
        X.VisOmit2diff.maps = omit2_maps - base_maps;
        X.VisOmit2diff.scale = 'minmax';
        X.VisOmit2diff.clims = [0 1];
        X.VisOmit2diff.colormap = RedWhiteBlue;
        
        base_maps = squeeze(EXP.Maps{mapVS_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit3_maps = squeeze(EXP.Maps{mapVS_idx(omit3idx)}.Tuning.meanrespModel(:,1,:));
        X.VisOmit3diff.maps = omit3_maps - base_maps;
        X.VisOmit3diff.scale = 'minmax';
        X.VisOmit3diff.clims = [0 1];
        X.VisOmit3diff.colormap = RedWhiteBlue;
        
        base_maps = squeeze(EXP.Maps{mapVS_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit4_maps = squeeze(EXP.Maps{mapVS_idx(omit4idx)}.Tuning.meanrespModel(:,1,:));
        X.VisOmit4diff.maps = omit4_maps - base_maps;
        X.VisOmit4diff.scale = 'minmax';
        X.VisOmit4diff.clims = [0 1];
        X.VisOmit4diff.colormap = RedWhiteBlue;
    end
    
    if any(strcmp(S, 'noOmitFullOmitMapdiff'))
        % Same omission-minus-base diff again, but on the full model's
        % prediction when fit WITHOUT the omission predictor
        % (mapVSPnoOmit_idx): isolates whether the full model, even
        % lacking explicit omission terms, still "leaks" some omission-
        % like signal through its other predictors (e.g. via the spatial
        % kernel itself), versus needing dedicated omission terms to
        % reproduce the effect.
        base_maps = squeeze(EXP.Maps{mapVSPnoOmit_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit2_maps = squeeze(EXP.Maps{mapVSPnoOmit_idx(omit2idx)}.Tuning.meanrespModel(:,1,:));
        X.noOmitVisOmit2diff.maps =omit2_maps - base_maps;
        X.noOmitVisOmit2diff.scale = 'minmax';
        X.noOmitVisOmit2diff.clims = [0 1];
        X.noOmitVisOmit2diff.colormap = RedWhiteBlue;
        
        base_maps = squeeze(EXP.Maps{mapVSPnoOmit_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit3_maps = squeeze(EXP.Maps{mapVSPnoOmit_idx(omit3idx)}.Tuning.meanrespModel(:,1,:));
        X.noOmitVisOmit3diff.maps =omit3_maps - base_maps;
        X.noOmitVisOmit3diff.scale = 'minmax';
        X.noOmitVisOmit3diff.clims = [0 1];
        X.noOmitVisOmit3diff.colormap = RedWhiteBlue;
        
        base_maps = squeeze(EXP.Maps{mapVSPnoOmit_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        omit4_maps = squeeze(EXP.Maps{mapVSPnoOmit_idx(omit4idx)}.Tuning.meanrespModel(:,1,:));
        X.noOmitVisOmit4diff.maps =omit4_maps - base_maps;
        X.noOmitVisOmit4diff.scale = 'minmax';
        X.noOmitVisOmit4diff.clims = [0 1];
        X.noOmitVisOmit4diff.colormap = RedWhiteBlue;
    end
    
    if any(strcmp(S, 'swapMap'))
        % Measured mean response (data) on the landmark-swap conditions
        % (the swapped-order trial variants, distinct from the omission
        % conditions above).
        X.swap23.maps = squeeze(EXP.Maps{mapDATA_idx(swap23idx)}.Tuning.meanrespModel(:,1,:));
        X.swap23.scale = 'minmax';
        X.swap23.clims = [0 1];
        X.swap23.colormap = flipud(gray(256));
        
        X.swap34.maps = squeeze(EXP.Maps{mapDATA_idx(swap34idx)}.Tuning.meanrespModel(:,1,:));
        X.swap34.scale = 'minmax';
        X.swap34.clims = [0 1];
        X.swap34.colormap = flipud(gray(256));
    end
    
    if any(strcmp(S, 'DataSwapMapDiff'))
        % swap-condition minus base, on the data, same logic as
        % DataOmitMapDiff above but for the swap conditions.
        X.swap23diff.maps = squeeze(EXP.Maps{mapDATA_idx(swap23idx)}.Tuning.meanrespModel(:,1,:))...
            - squeeze(EXP.Maps{mapDATA_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        X.swap23diff.scale = 'minmax';
        X.swap23diff.clims = [0 1];
        X.swap23diff.colormap = RedWhiteBlue;
        
        X.swap34diff.maps = squeeze(EXP.Maps{mapDATA_idx(swap34idx)}.Tuning.meanrespModel(:,1,:))...
            - squeeze(EXP.Maps{mapDATA_idx(baseidx)}.Tuning.meanrespModel(:,1,:));
        X.swap34diff.scale = 'minmax';
        X.swap34diff.clims = [0 1];
        X.swap34diff.colormap = RedWhiteBlue;
    end
    
    Nbins = size(squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(:,1,:)), 2);
    x_pos = linspace(0,200,Nbins);

    % Cell selection: default restricts to spatially-selective cells
    % (the 'all' and 'omission' modifier tags override this; see header
    % doc). show_omit takes priority over show_all when both are absent.
    spatial_mask = goodcells;
    if ~show_all && ~show_omit %showing spatially selective cells only by default
        spatial_mask = spatial_mask & signicells & goodLLHcells;
    elseif show_omit
        spatial_mask = spatial_mask & goodOmitcells;
    end
    N = sum(spatial_mask);
    
    % Apply the cell mask to every requested X.* map. The try/catch
    % handles two possible orientations of X.(fnames{k}).maps: most maps
    % are [cells x bins] (first branch), but a few computed via squeeze()
    % above can come out as [bins x cells] for single-cell-dimension
    % cases, hence the fallback transpose-and-index in the catch branch.
    fnames = fieldnames(X);
    for k = 1:numel(fnames)
        try
            X.(fnames{k}).maps = X.(fnames{k}).maps(spatial_mask,:);
        catch
            X.(fnames{k}).maps = X.(fnames{k}).maps(:, spatial_mask)';
        end
    end
    
    % SpatialWeights: the metric shown in the leftmost "spatial weight"
    % sidebar of the figure (red/black ticks). Defaults to LLHrel; when
    % 'omission' was requested (show_omit), switches to LLHrel_omit so
    % the sidebar reflects omission-component strength instead of
    % overall spatial-coding strength.
    if ~show_omit
        SpatialWeights = LLHrel(spatial_mask);
    else
        SpatialWeights = LLHrel_omit(spatial_mask);
    end
    Spatialcells = signicells(spatial_mask) & goodLLHcells(spatial_mask);
    
    sessionid = EXP.Spk.series;
    if isfield(EXP.Spk, 'animal')
        sessionid = sessionid + 10*EXP.Spk.animal;
    end
    
    [sess_vals, ~, sess_uniq] = unique(sessionid);
    sess_sel = sess_uniq(spatial_mask);
    
    % === Ordering rows so successive kernels are most similar ===
    % Xref: the reference kernel/map (possibly several tags concatenated
    % side by side, if ref has more than one entry - though typically
    % ref is a single tag, the first one in S) used purely to compute the
    % cell ordering; not itself drawn differently from other columns.
    Xref = [];
    for k = 1:numel(ref)
        Xref = cat(2, Xref, X.(ref{k}).maps);
    end
    
    % Row-wise z-score so the ordering reflects kernel SHAPE, not
    % overall gain/offset differences between cells.
    mu = mean(Xref, 2, 'omitnan');
    sd = std(Xref, 0, 2, 'omitnan');
    Xref = (Xref - mu) ./ sd;
    if strcmp(orderby, 'similarity')
        optiorder = orderKernelsBySimilarity(Xref, clust_nPCs, dist_metric, dist_method);
    elseif strcmp(orderby, 'peak')
        [~, peakIdx] = max(Xref, [], 2, 'omitnan');
        [~, optiorder] = sort(peakIdx, 'ascend', 'MissingPlacement', 'last');
    end
    % Apply the computed ordering to every map/kernel column and to the
    % auxiliary per-cell vectors (weights, spatial-cell flag, session id)
    % so everything stays aligned to the same row order.
    for k = 1:numel(fnames)
        X.(fnames{k}).maps = X.(fnames{k}).maps(optiorder,:);
    end
    SpatialWeights = SpatialWeights(optiorder);
    Spatialcells = Spatialcells(optiorder);
    sess_ord = sess_sel(optiorder);

    % From here on, ref collapses to a single string (the first/primary
    % reference tag) for indexing X.(ref) directly.
    ref = ref{1};
    mu = mean(X.(ref).maps, 2, 'omitnan');
    sd = std(X.(ref).maps, 0, 2, 'omitnan');
    Xref_norm = (X.(ref).maps - mu) ./ sd;
    
    % === Build similarity groups of K_n cells by lowest variance (non-overlapping) ===
    % Greedy algorithm: repeatedly find the contiguous (in the now-
    % ordered list) window of Kn(1)-Kn(2) cells with the lowest internal
    % variance (i.e. the most mutually similar adjacent cells), remove
    % those cells from the pool, and repeat - up to max_examples groups,
    % or until fewer than Kn(1) cells remain ungrouped.
    remaining = 1:N;
    keep_on = true;
    groups = {};
    while numel(remaining) >= min(Kn) && keep_on && numel(groups)<=max_examples
        best_sse = inf; best_window = [];
        % slide a contiguous window over the remaining index list (contiguity in the ordered list)
        mu = [];
        for k = Kn(1):Kn(2)
            for s = 1:(numel(remaining)-k+1)
                win_idx = remaining(s:s+k-1);
                if sum(diff(win_idx) > 1) > 0
                    continue;
                end
                M = Xref_norm(win_idx,:);
                mu = mean(M,1,'omitnan');
                sse = mean(sum((M - mu).^2, 2, 'omitnan'), 'omitnan') / numel(win_idx)^0.5;
                if sse < best_sse
                    best_sse = sse;
                    best_window = win_idx;
                end
            end
        end
        if isempty(mu)
            keep_on = false; % no valid window found this pass; stop grouping
        else
            groups{end+1} = best_window; %#ok<AGROW>
            % remove selected indices from remaining
            mask_rem = true(1,numel(remaining));
            mask_rem(ismember(remaining, best_window)) = false;
            remaining = remaining(mask_rem);
        end
    end
    Ng = numel(groups);

    figure;
    
    % Colors for groups
    % Trim the hsv colormap's saturated extremes (full red wrap-around)
    % so adjacent groups remain visually distinguishable.
    cmap = colormap(hsv);
    cmap = cmap(32:224,:);
    group_cmap = cmap(round(linspace(1, size(cmap, 1), Ng)),:);
    
    sess_cmap = cmap(round(linspace(1, size(cmap, 1), numel(sess_vals))),:);

    % === Layout: first column = heatmap (spans all rows),
    %     last columns = per-group mean+SE subplots (arranged in grid with nrow<=5),
    % Tile grid columns: 1 "spatial weight" sidebar + Nmaps heatmap
    % columns (col_offset = Nmaps+1), followed by ncol_groups columns of
    % per-group mean(+hist) subplots, arranged in up to 5 rows so the
    % group panels don't get too short when there are many groups.
    nrow = min(5, Ng);
    ncol_groups = ceil(Ng / nrow);
    Nmaps = numel(fnames);
    col_offset = Nmaps + 1;
    
    ncol_total  = col_offset + ncol_groups; % 1 heatmap col + (mean) per group-column
    if show_hist
        ncol_total = ncol_total + ncol_groups; %(mean,hist) per group-column
    end
    tl = tiledlayout(nrow, ncol_total, 'TileSpacing','compact','Padding','compact');

    % ---- Weights and heatmaps of reordered kernels ----
    % Spatial weights sidebar: one tick mark per cell (row), at its
    % SpatialWeights value (LLHrel or LLHrel_omit); red = spatially
    % selective (Spatialcells), black = not. The gray line is a local
    % running median (+/-5 neighboring cells in the ordered list) to show
    % the overall trend along the ordering; the vertical line at x=1
    % marks the "no improvement" reference value.
    axW = nexttile(tl, [nrow 1]); hold(axW,'on');
    sw_median = NaN(size(SpatialWeights));
    for j = 1:numel(SpatialWeights)
        if Spatialcells(j)
            c = 'r';
        else
            c = 'k';
        end
        plot(SpatialWeights(j)*[1 1], j + [-0.5 0.5], c);
        if j > 5 && j < numel(SpatialWeights) - 5
            sw_median(j) = median(SpatialWeights(j-5:j+5));
        end
    end
    plot(sw_median, 1:numel(SpatialWeights), 'Color', [.5 .5 .5]);
    hold on;plot([1 1], [1 numel(SpatialWeights)])
    
    % kernels heatmaps
    for k = 1:numel(fnames)
        axH = nexttile(tl, [nrow 1]); hold(axH,'on');
        if strcmp(X.(fnames{k}).scale,'log')
            imagesc(axH, x_pos, 1:N, log(X.(fnames{k}).maps), X.(fnames{k}).clims);
        elseif strcmp(X.(fnames{k}).scale,'minmax')
            mx = max(X.(fnames{k}).maps, [], 2);
            mn = min(X.(fnames{k}).maps, [], 2);
            maps = (X.(fnames{k}).maps - mn) ./ (mx - mn);
            imagesc(axH, x_pos, 1:N, maps, X.(fnames{k}).clims);
        else
            imagesc(axH, x_pos, 1:N, X.(fnames{k}).maps, X.(fnames{k}).clims);
        end
    
        if strcmp(fnames{k}, ref)
            % These two annotations only get drawn once, on the
            % reference column's heatmap (since the group/session
            % structure was computed from the reference kernel).
            % Group brackets: a colored vertical bar + label to the
            % right of the heatmap, marking the row-range belonging to
            % each contiguous group formed above.
            for g = 1:Ng
                idx = groups{g};
                y1 = min(idx); y2 = max(idx);
                x0 = x_pos(end) + 2;  % place each group at distinct x for visibility
                plot(axH, [x0 x0], [y1 y2], '-', 'Color', group_cmap(g,:), 'LineWidth', 3);
                txt = sprintf('group #%.0f', g);
                text(x0 + 2, (y1 + y2) / 2, txt, ...
                 'HorizontalAlignment','left','VerticalAlignment','middle','FontWeight','bold');
            end
    
            % Session scatter strip: a column of colored dots to the left
            % of the heatmap, one color per recording session, so it's
            % visible at a glance whether the ordering/grouping above
            % clusters cells by session (a potential confound) or mixes
            % sessions freely.
            for s = 1:numel(sess_vals)
                vals = find(sess_vals(sess_ord) == sess_vals(s));
                scatter((x_pos(1) - 5 - 5*s)*ones(size(vals)), vals,...
                    'filled', 'MarkerFaceColor', sess_cmap(s,:), 'MarkerFaceAlpha', ptAlpha,...
                    'MarkerEdgeColor', 'none');
            end
        end
    
        set(axH,'YDir','normal');
        colormap(axH, X.(fnames{k}).colormap);
        colorbar(axH);
        xlabel(axH,'Position (cm)'); ylabel(axH,'Cells (reordered)');
        title(axH,fnames{k});
        box(axH,'off'); set(axH,'TickDir','out');
    end

    % Precompute y-limits shared across mean+SE plots
    % Computing this up front (before drawing any of the per-group
    % subplots below) means every group's mean+SE subplot, across every
    % kernel/map tag, shares one common y-axis range - necessary for the
    % per-group panels to be visually comparable to each other.
    yMin = inf; yMax = -inf;
    mu = []; se = [];
    for k = 1:numel(fnames)
        mu.(fnames{k}) = cell(Ng,1); se.(fnames{k}) = cell(Ng,1);
        for g = 1:Ng
            idx = groups{g};
            M = X.(fnames{k}).maps(idx,:);%     Xord_scaled(idx,:);
            mu_temp = mean(M,1,'omitnan');
            se_temp = std(M,0,1,'omitnan') / max(numel(idx),1)^0.5;
            mu.(fnames{k}){g} = mu_temp;
            se.(fnames{k}){g} = se_temp;
            yMin = min(yMin, min(mu_temp - se_temp, [], 'omitnan'));
            yMax = max(yMax, max(mu_temp + se_temp, [], 'omitnan'));
        end
    end
    
    if ~isfinite(yMin), yMin = 0; end
    if ~isfinite(yMax) || yMax<=yMin, yMax = yMin + 1; end
    yLimMeans = [yMin yMax]*1.05;

    % ---- Per-group mean (+SE) subplots ----
    for g = 1:Ng
        % tile row/column positions
        rowIdx      = mod(g-1, nrow) + 1;
        groupColIdx = ceil(g / nrow);
        if show_hist
            colMean = col_offset + 2*(groupColIdx-1) + 1;    % mean plot column
            colHist = col_offset + 2*(groupColIdx-1) + 2;    % histogram column
        else
            colMean = col_offset + (groupColIdx-1) + 1;    % mean plot column
            colHist = 0;
        end

        % indices for this group (in reordered space)
        idx = groups{g};
        color_g = group_cmap(g,:);

        % ---- Mean (+SE) subplot ----
        % All requested kernel/map tags are overlaid on the same axes for
        % this group: the reference tag gets the group's own color plus
        % an SE shaded band (it's the tag the group was formed from);
        % every other tag is drawn as a thin colored line (RGB cycled via
        % mod(k,3)) purely so its SHAPE can be visually compared to the
        % reference, without claiming the same statistical weight.
        axG = nexttile(tl, (rowIdx-1)*ncol_total + colMean); hold(axG,'on');
        for k = 1:numel(fnames)
            if ~strcmp(fnames{k}, ref)
                c = zeros(1,3);
                c(mod(k, 3)+1) = 1;
                plot(axG, x_pos, mu.(fnames{k}){g}, 'LineWidth',1.6, 'Color',c);
            end
        end
        seplot(axG, x_pos, mu.(ref){g}, se.(ref){g}, 0.2);       % gray SE
        plot(axG, x_pos, mu.(ref){g}, 'LineWidth',1.6, 'Color', color_g);
        xlim(axG,[0 200]); ylim(axG, yLimMeans);
        if rowIdx == nrow, xlabel(axG,'Position (cm)'); end
        ylabel(axG, sprintf('Group %d (n=%d)', g, numel(idx)));
        title(axG, 'Mean kernel');
        box(axG,'off'); set(axG,'TickDir','out');

        % ---- Session ID histogram subplot (next to the mean) ----
        % Only drawn if the 'hist' modifier tag was requested; same
        % "check for session confounds" purpose as the scatter strip on
        % the reference heatmap column above, but per-group instead of
        % across all cells at once.
        if colHist > 0
            axS = nexttile(tl, (rowIdx-1)*ncol_total + colHist); hold(axS,'on');
            sess_g = sess_ord(idx); % session ids for this group
            % histogram with session bins
            if ~isempty(sess_g)
                histogram(axS, sess_g, 'BinMethod','integers', ...
                          'FaceColor', color_g, 'EdgeColor','none');
            else
                histogram(axS, nan, 'BinMethod','integers', 'FaceColor', color_g, 'EdgeColor','none');
            end
            title(axS, 'Sessions');
            xlabel(axS, 'session id'); ylabel(axS, '# cells');
            xtickangle(axS, 45);
            box(axS,'off'); set(axS,'TickDir','out');
            grid(axS,'on');
        end
    end

    % Overall title
    sgtitle(tl, sprintf('SpaceKernels seriation (N=%d, K_n=%d): heatmap, group means, ranges', N, Kn));
end

% =============================================================================
% LOCAL HELPER FUNCTIONS
% =============================================================================

% Helper to plot one category
% Used by the 'LLHi-w/oSpace' figure: draws one category of cells as a
% pair of x-positions (x1, x2) connected by a thin line per cell, with
% dots at each end - i.e. a "paired observations" plot showing how each
% cell's value changes between the two conditions/models being compared.
function hLine = plot_cat(x1, y1, x2, y2, lineColor, alpha, dispName)
    % X positions (two columns) with tiny jitter for visibility
    if isempty(y1), hLine = gobjects(1); return; end
    % Lines
    for i = 1:numel(y1)
        plot([x1 x2], [y1(i) y2(i)], '-', ...
             'Color', lineColor, 'LineWidth', 0.1, ...
             'HandleVisibility','off', ...
             'AlignVertexCenters','on', 'Clipping','on', 'LineStyle','-');
    end
    % Dots
    h1 = scatter(x1 + 0*y1, y1, ptSize, 'o', ...
        'MarkerEdgeColor', lineColor, 'MarkerFaceColor', lineColor, ...
        'MarkerFaceAlpha', alpha, 'DisplayName', dispName);
    scatter(x2 + 0*y2, y2,   ptSize, 'o', ...
        'MarkerEdgeColor', lineColor, 'MarkerFaceColor', lineColor, ...
        'MarkerFaceAlpha', alpha, 'HandleVisibility','off');
    hLine = h1; % return handle for legend
end
    
end

function plotViolin(ax, x, y, edges, binw, support, c)
% plotViolin  Draw a violin plot of y, binned by x into the bins given by
% edges, on existing axes ax.
%   x, y      paired data (e.g. a position value and an LLHrel value, one
%             per cell)
%   edges     bin edges along x
%   binw      (optional) bin width to actually use for grouping points
%             into each violin; if smaller than diff(edges), bins can
%             overlap less than their edges would suggest (effectively
%             narrowing the data window per violin while keeping violins
%             spaced at the edges' centers). Defaults to diff(edges).
%   support   (optional) [min max] y-range passed to ksdensity, fixing
%             the density estimate's support across all bins for
%             comparability; computed from data range with a 5% margin
%             if not given.
%   c         (optional) violin fill/line color, default black.
% Bins with >=3 points get a real kernel-density violin (mirrored
% left/right) plus a median tick; bins with 1-2 points fall back to
% plotting raw points instead, since a density estimate from so few
% points would be meaningless.
    if nargin < 7
        c = [0, 0, 0];
    end
    if nargin < 6
        support = [];
    end
    if nargin < 5
        binw = [];
    end
    nbins = numel(edges) - 1;
    if isempty(binw)
        binw    = diff(edges);
    else
        binw = binw * ones(size(edges(1:end-1)));
    end
    centers = edges(1:end-1) + diff(edges) / 2;
    w_plt = min(binw, diff(edges));
    if isempty(support)
        mi = min(y)+eps;
        ma = max(y)+eps;
        ra = ma - mi;
        mi = mi - 0.05*ra;
        ma = ma + 0.05*ra;
        support = [mi, ma];
    end
    for b = 1:nbins
        inb = (x >= centers(b) - binw(b) / 2) & (x < centers(b) + binw(b) / 2);
        yb  = y(inb);
        if numel(yb) >= 3
            x0 = centers(b);% + ((j - (numel(shape_vals)+1)/2) / max(numel(shape_vals),1)) * 0.6*binw(b);
            % kernel density of r (limit support to [0,1] to match axis)
            [f, yi] = ksdensity(yb, 'Support', support);
            if all(~isfinite(f)), continue; end
            f  = f / max(f);                          % normalize width
            xi = x0 + 0.3*w_plt(b) * f;                % half-width violin
            % Mirror the density estimate (xi) around x0 to form the
            % violin's left and right halves (2*x0-xi reflects xi).
            patch(ax, [xi, fliplr(2*x0 - xi)], [yi, fliplr(yi)], c, ...
                  'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility','off');
            % median line
            med = median(yb,'omitnan');
            plot(ax, [x0-0.25*w_plt(b), x0+0.25*w_plt(b)], [med med], '-', ...
                 'Color', c, 'LineWidth', 1.1, 'HandleVisibility','off');
        elseif numel(yb) > 0
            % fallback for tiny samples: show points
            x0 = centers(b);
            plot(ax, repmat(x0, size(yb)), yb, '.', 'Color', c, 'MarkerSize', 6, 'HandleVisibility','off');
        end
    end
end

function [T, Z, leafOrder, colorCutoff] = ClusterSpaceKernels(spaceKernels, n_Clusters, n_PCs, smth_win, cluster_distance, cluster_method)
% Cluster spatial kernels with hierarchical clustering.
% Returns T (cluster labels), Z (linkage), leafOrder (optimal dendrogram order),
% and colorCutoff (dendrogram cutoff level).
%
% IMPORTANT: this function actually uses TWO different clustering
% algorithms for two different purposes. Hierarchical clustering
% (linkage/cluster) is used only to build the dendrogram structure (Z,
% leafOrder, colorCutoff) for visualization; the cluster LABELS actually
% returned in T are computed separately by k-means on the very last line.
% This means the colors/grouping implied by the dendrogram's cut height
% (colorCutoff) do not necessarily match the k-means group assignment T
% used elsewhere in the figure - the dendrogram is illustrative of
% similarity structure, not a literal visualization of how T was formed.

% ----- defaults -----
if nargin < 2 || isempty(n_Clusters), n_Clusters = 8; end
if nargin < 3 || isempty(n_PCs),      n_PCs      = 20; end
if nargin < 4 || isempty(smth_win),   smth_win   = 0; end

% ----- row-wise z-score with guards -----
% Normalizes each cell's kernel to zero mean / unit variance across
% position bins, so clustering groups cells by kernel SHAPE rather than
% by overall response magnitude.
mu = mean(spaceKernels, 2, 'omitnan');
sd = std(spaceKernels, 0, 2, 'omitnan');          % N-1 normalization
spaceKernels = (spaceKernels - mu) ./ sd;

% optional smoothing (omit if win==0 or function missing)
if smth_win > 0 && exist('SpecialSmooth', 'file') == 2
    maps_space = SpecialSmooth(spaceKernels, [0 smth_win/100], [1, 100]);
else
    maps_space = spaceKernels;
end
% Re-z-score after smoothing (smoothing can shrink the variance, so this
% restores unit variance per row before the PCA/clustering step).
% NOTE: re-derives mu/sd from maps_space but applies them to the
% un-smoothed spaceKernels on the right-hand side - if smth_win>0 this is
% deliberately mixing the smoothed normalization with the raw values;
% double-check this is intended if you change the smoothing behavior.
mu = mean(maps_space, 2, 'omitnan');
sd = std(maps_space, 0, 2, 'omitnan');          % N-1 normalization
maps_space = (spaceKernels - mu) ./ sd;

% maps_space = cat(2, spaceKernels, maps_space);

% PCA (cap components to valid maximum)
% Reduces dimensionality before clustering/distance computation; ncomp is
% capped so it never exceeds the number of position bins or (cells-1),
% both of which would make the requested PCA invalid.
if n_PCs > 0
    ncomp = min([n_PCs, size(maps_space,2), size(maps_space,1)-1]);
    if ncomp > 0
        [~, score] = pca(maps_space, 'NumComponents', ncomp);
        maps_clust = score;
    else
        maps_clust = maps_space;
    end
else
    maps_clust = maps_space;
end
mu = mean(maps_clust, 2, 'omitnan');
sd = std(maps_clust, 0, 2, 'omitnan');          % N-1 normalization
maps_clust = (maps_clust - mu) ./ sd;

% distances + linkage
% Y: pairwise distance vector (for optimalleaforder); Z: the hierarchical
% linkage tree used purely to build the dendrogram and its leaf ordering.
Y = pdist(maps_clust, cluster_distance);
Z = linkage(maps_clust, cluster_method, cluster_distance);

leafOrder = optimalleaforder(Z, Y);
% T = cluster(Z, 'MaxClust', n_Clusters);
% This hierarchical cluster() call's result is used ONLY to derive
% colorCutoff (an approximate dendrogram coloring threshold below); its
% actual cluster labels T are immediately discarded and overwritten by
% k-means at the end of this function.
T = cluster(Z, 'Cutoff', 2.5, 'Depth', 100);

% color cutoff
% Picks a linkage-distance threshold that would produce exactly K=max(T)
% clusters if the dendrogram were cut there - purely a cosmetic choice
% for where to color-split the dendrogram in the figure; has no effect
% on the actual k-means cluster assignment below.
K = max(T);
Zs = sort(Z(:,3), 'ascend');
idx = max(numel(Zs) - (K-1), 1);
colorCutoff = Zs(idx);

% Final cluster assignment actually used everywhere else in the
% 'Cluster-SpaceKernels' figure: k-means on the (possibly PCA-reduced,
% z-scored) kernels, requesting exactly n_Clusters groups.
T = kmeans(maps_clust,n_Clusters);
end

function optimalOrder = orderKernelsBySimilarity(spaceKernels, n_PCs, metric, method)
% Returns a row permutation that makes successive kernels as close as possible.
% Used by the 'Ordered-Kernels-*' figures (orderby='similarity', the
% default) to seriate cells so that adjacent rows in the heatmap have
% maximally similar kernel shapes - this is what makes the
% sliding-window grouping (in the main figname block above) effective at
% finding genuinely similar contiguous runs of cells.

% method   = 'farthest';%'ward';
% metric = 'cosine';   % ward/centroid/median require Euclidean

% Row-wise z-score (first pass, on the raw kernel) so the optional PCA
% below operates on shape-normalized data rather than raw gain/offset
% differences between cells.
mu = mean(spaceKernels, 2, 'omitnan');
sd = std(spaceKernels, 0, 2, 'omitnan');          % N-1 normalization
spaceKernels = (spaceKernels - mu) ./ sd;

X = spaceKernels;
if n_PCs > 0
    ncomp = min([n_PCs, size(X,2), size(X,1)-1]);
    if ncomp > 0
        [~, score] = pca(X, 'NumComponents', ncomp);
        X = score;
    end
end

% Row-wise zscore (robust) so distance focuses on SHAPE, not gain/offset
mu = mean(X, 2, 'omitnan');
sd = std( X, 0, 2, 'omitnan'); sd(~isfinite(sd) | sd==0) = eps;
X = (X - mu) ./ sd;
% Cells whose kernel z-scores to all-NaN (e.g. a fully flat/degenerate
% kernel after PCA) can't be placed in the distance-based ordering;
% valid_idx are clustered/ordered normally, nan_idx are simply appended
% at the end of the final order (rather than dropped entirely).
valid_idx = find(~isnan(mu));
nan_idx = find(isnan(mu));

% Distances + linkage
Y = pdist(X(valid_idx,:), metric);
Z = linkage(X(valid_idx,:), method, metric);

% Seriation that minimizes sum of adjacent distances
% optimalleaforder finds the leaf ordering of the hierarchical tree Z
% that minimizes the sum of distances between adjacent leaves - this is
% the actual "similarity ordering" returned by this function (the
% clustering tree Z/Y themselves are discarded after this point, only
% the resulting order survives).
optimalOrder = optimalleaforder(Z, Y);        % indices into idx_keep
optimalOrder = valid_idx(optimalOrder);
optimalOrder = cat(1, optimalOrder, nan_idx);

end