% RFMapping_Gaussian2DFit_current.m
%
% Fits an anisotropic 2D Gaussian (no rotation, no baseline term) to each
% bouton's mean response across the RF mapping grid -- matching
% Timplalexi et al.:
%
%   G(x,y) = A * exp( -( (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2) ) )
%
% Trust criteria (REVISED):
%   1. r2Floor -- a LOW screening floor only (0.1), used purely to skip catastrophic non-fits
%      before spending 100 bootstrap refits on them. This is NOT the trust threshold.
%   2. trustR2Threshold + bootstrap percentile CI -- the REAL trust criterion, set and applied
%      entirely AFTER fitting (see the sweep section near the end), from the raw per-bouton
%      bootstrap R^2 distribution stored during fitting. A bouton is "robust at threshold th" if
%      the (1-ciLevel) lower percentile of its bootstrap R^2 distribution is itself >= th -- i.e.
%      we are ciLevel-confident (e.g. 95%) that the TRUE R^2 clears th, not merely that "most"
%      resamples happened to. robustCIThreshold is MATCHED to trustR2Threshold (not decoupled to
%      a separately-set lower value) -- both the point estimate AND the CI lower bound must clear
%      the SAME threshold, which is what was actually validated (n=82, sigAz=15.34+/-9.17,
%      sigEl=8.59+/-3.54; a decoupled/looser CI check gave a different, less-validated result).
%
%   3. isDegenerateSigma -- a DIRECT, per-axis identifiability check: is the fitted sigma >= a
%      FRACTION (degenerateSigmaFrac = 1/turnoverWidths, turnoverWidths = 2.5) of the actual
%      tested span (range(azRange) / range(elRange)) on that axis? This fraction is derived from
%      Gaussian geometry -- seeing a Gaussian turn over on both sides of its peak requires
%      roughly 2-3 sigma-widths of coverage, not just 1 -- and is applied IDENTICALLY to both
%      axes, since it is a property of the Gaussian shape, not of grid density. The fitting
%      bound (ub(4)/ub(5) = 3*range) is deliberately loose so the optimizer never hits a wall --
%      it is NOT a claim about plausible RF size, and a fit can sit nowhere near that bound while
%      still being unidentifiable. Visual inspection of fits with sigma at ~45-60% of the tested
%      range showed clear monotonic, edge-going responses along that axis -- never turning over
%      anywhere in the grid -- confirming the coverage logic above. This check GATES isTrusted
%      alongside the R^2/bootstrap-CI criteria (previously this exclusion was accidentally left
%      commented out in the commit block -- FIXED; without it, the trusted count was 114 instead
%      of the validated 82, silently re-including the monotonic/degenerate fits this check exists
%      to catch).
%
%      Axis-specific grid SPARSITY (e.g. few elevation positions) is a separate concern from
%      coverage and is NOT folded into degenerateSigmaFrac. It is tracked instead as a
%      reporting-only diagnostic -- nPointsInWindowAz/El and isSparseWindowAz/El, stored on both
%      gaussFitResults and allRFMapping -- so sparse-axis fits remain visible as a caveat rather
%      than silently passing a loosened threshold.
%
%   4. isBeyondRange / isAtBoundary -- per-axis REPORTING-ONLY checks on the fitted CENTER (not
%      sigma), using the stimulus-footprint-extended range (grid centers + stimulus radius,
%      clipped to actual screen limits -- see screenAzLimits/screenElLimits below) rather than
%      the bare grid-center coordinates:
%        - isBeyondRange*: fitted center is outside the visual space the outermost stimulus
%          actually covered. NOT GATED on isTrusted: empirically, observed overshoots were only
%          ~5-10 deg past the extended edge -- within roughly one grid-step, i.e. within the
%          fit's own center-estimate precision rather than genuine extrapolation. Gating on this
%          excluded plausible boundary-adjacent preferences for reasons that don't track
%          reliability (the same failure mode as the old width-biased center-coverage bootstrap).
%        - isAtBoundary*: fitted center sits AT/NEAR the outermost TESTED GRID position, but
%          still within the stimulated range. Also reporting-only -- may reflect a genuine
%          preference for the most peripheral position actually tested.
%      Both flags remain available on gaussFitResults/allRFMapping for inspection; only
%      isDegenerateSigma (a genuine identifiability problem) and the R^2/bootstrap-CI criteria
%      gate isTrusted.
%
% RESPONSE WINDOW CONSISTENCY: both the point estimate and every
% bootstrap resample now call the SAME local function
% (computeGridRespVec, defined at the end of this file) to build the
% response vector from allRFMapping.rawAlignedTrials, restricted to
% respIdx (respWin = [0.1 3], defined explicitly below -- keep in sync
% with the main pooling script). This guarantees the point estimate and
% its own bootstrap are always computed identically, with no possibility
% of drifting onto different windows.
%
% REQUIRES (already in workspace from the main pooling script):
%   allRFMapping   - pooled bouton struct array, with .rawAlignedTrials
%                    ([nEl x nAz] cell of [Trials x Time] matrices)
%   uAz, uEl_plot  - grid position vectors (uEl_plot rows, uAz columns)
%   timeVector     - time vector matching the Time dimension above, used
%                    to compute respIdx from respWin
%
% REQUIRES: Optimization Toolbox (lsqcurvefit). Statistics Toolbox (prctile).

%% params
rng(1, 'twister'); % FIXED SEED -- without this, every run draws different bootstrap resamples,
                    % which can change WHICH boutons pass a given criterion.

r2Floor        = 0.1;  % LOW screening floor only -- excludes catastrophic non-fits before
                        % spending 100 bootstrap refits on them. NOT the trust threshold; do not
                        % raise this to tighten your trusted set -- use trustR2Threshold below
                        % (set post-hoc, in the commit block) instead.
nBootstrap     = 1000;  % number of within-position trial resamples for the R^2 bootstrap distribution

% NOTE: azRange/elRange are the raw grid CENTER coordinates. But each stimulus
% is a stimDiameterDeg-diameter grating, so a position centered at e.g. -70 deg already visually
% stimulates out to roughly -70 - 15 = -85 deg. The tested visual field is therefore wider than
% the array of grid centers implies. stimRadius is added when deciding isBeyondRange below, so a
% fitted center isn't wrongly flagged as "extrapolated" when it actually falls within visual
% space the outermost stimulus already covered.
stimDiameterDeg = 30; % set to match your actual grating diameter
stimRadius = stimDiameterDeg / 2;

% Screen physically clips the stimulus near the edges -- a grating centered at the
% most-negative azimuth position does NOT fully extend by stimRadius, because part of it is
% drawn off-screen.
screenAzLimits = [-80, 20];   % NaN = no clipping on this side (extends past the stimulus footprint)
screenElLimits = [-25, 40];
edgeTolFrac = 0.5; % fraction of one grid-step counted as "at/near boundary" for isAtBoundary
                    % (reporting only -- a real preference for the most peripheral tested
                    % position should not be excluded just for landing there)

respWin = [0.1 3];
if ~exist('timeVector', 'var')
    error('timeVector not found -- run the main pooling script first (needed to compute respIdx from respWin).');
end
respIdx = timeVector >= respWin(1) & timeVector <= respWin(2);

azStep  = mean(diff(sort(uAz)));
elStep  = mean(diff(sort(uEl_plot)));
azRange = [min(uAz) max(uAz)];
elRange = [min(uEl_plot) max(uEl_plot)];

% Stimulus-extended range: a fitted center is only genuine extrapolation
% (isBeyondRange) if it falls outside the visual space the outermost stimulus actually covered,
% not just outside the bare grid-center coordinates. Clipped by the actual screen limits above.
% Used for REPORTING/inspection only (isBeyondRange is not gated -- see header note #4).
azExtLow  = max(azRange(1)-stimRadius, screenAzLimits(1));
azExtHigh = azRange(2)+stimRadius;
if ~isnan(screenAzLimits(2)), azExtHigh = min(azExtHigh, screenAzLimits(2)); end

elExtLow  = elRange(1)-stimRadius;
if ~isnan(screenElLimits(1)), elExtLow = max(elExtLow, screenElLimits(1)); end
elExtHigh = elRange(2)+stimRadius;
if ~isnan(screenElLimits(2)), elExtHigh = min(elExtHigh, screenElLimits(2)); end

azRangeExt = [azExtLow, azExtHigh];
elRangeExt = [elExtLow, elExtHigh];

[AzGrid, ElGrid] = meshgrid(uAz, uEl_plot);   % [nEl x nAz]
xyList = [AzGrid(:), ElGrid(:)];

gaussFit2D = @(p, xy) p(1) * exp( -( (xy(:,1)-p(2)).^2 ./ (2*p(4)^2) + (xy(:,2)-p(3)).^2 ./ (2*p(5)^2) ) );
% p = [A, x0, y0, sigmaX, sigmaY]

fitOpts = optimoptions('lsqcurvefit', 'Display', 'off');
lb = [0,   azRange(1)-azStep, elRange(1)-elStep, azStep/4,        elStep/4];
ub = [Inf, azRange(2)+azStep, elRange(2)+elStep, 3*range(azRange), 3*range(elRange)];

% fit only the responsive rois
% fitBoutonIdx = 1:numel(allRFMapping);
fitBoutonIdx = respIdxList;
numBoutons   = numel(fitBoutonIdx);

staleFields = {'gaussFit_A','gaussFit_Az0','gaussFit_El0','gaussFit_sigmaAz','gaussFit_sigmaEl', ...
    'gaussFit_R2', 'gaussFit_isTrusted', ...
    'gaussFit_bootR2LowerCI', 'gaussFit_isRobust', ...
    'gaussFit_isDegenerateSigmaAz', 'gaussFit_isDegenerateSigmaEl', 'gaussFit_isDegenerateSigma', ...
    'gaussFit_nPointsInWindowAz', 'gaussFit_nPointsInWindowEl', ...
    'gaussFit_isSparseWindowAz', 'gaussFit_isSparseWindowEl', ...
    'gaussFit_isBeyondRangeAz', 'gaussFit_isBeyondRangeEl', 'gaussFit_isBeyondRange', ...
    'gaussFit_isAtBoundaryAz', 'gaussFit_isAtBoundaryEl', 'gaussFit_isAtBoundary', ...
    'gaussFit_pShuffle', 'gaussFit_isSignificant', ... % old fields, cleared if left over from a prior run (shuffle test removed)
    'gaussFit_fracR2AboveFloor', ... % old field name, ditto
    'gaussFit_bootCenterStdAz','gaussFit_bootCenterStdEl', ... % older field names, ditto
    'gaussFit_pctWithinStepAz', 'gaussFit_pctWithinStepEl'};
for iROI = fitBoutonIdx(:)'
    for f = staleFields
        if isfield(allRFMapping, f{1})
            allRFMapping(iROI).(f{1}) = [];
        end
    end
end

% bootR2 is stored RAW (full 1000-value distribution) per bouton, not just a summary at one
% fixed threshold. This is what lets trustR2Threshold and ciLevel be set/swept freely afterward,
% with NO REFITTING -- see the sweep section near the end of this script.
gaussFitResults = struct('iROI', {}, 'isResponsive', {}, 'A', {}, 'x0', {}, 'y0', {}, 'sigmaX', {}, 'sigmaY', {}, 'R2', {}, ...
    'isDegenerateSigmaAz', {}, 'isDegenerateSigmaEl', {}, 'isDegenerateSigma', {}, ...
    'nPointsInWindowAz', {}, 'nPointsInWindowEl', {}, 'isSparseWindowAz', {}, 'isSparseWindowEl', {}, ...
    'isBeyondRangeAz', {}, 'isBeyondRangeEl', {}, 'isBeyondRange', {}, ...
    'isAtBoundaryAz', {}, 'isAtBoundaryEl', {}, 'isAtBoundary', {}, ...
    'bootR2', {}, 'nBootValid', {});

fprintf('Fitting 2D Gaussian RFs for %d ROIs (visually responsive rois; screening floor R^2 >= %.2f)...\n', numBoutons, r2Floor);
fprintf('Running bootstrap R^2 distribution (%d resamples) on all fits above the screening floor -- this may take a while.\n', nBootstrap);
fprintf('Trust threshold (trustR2Threshold) and bootstrap CI level are NOT set here -- see the commit block at the end of this script.\n');

allR2 = nan(numBoutons, 1);

for k = 1:numBoutons
    iROI = fitBoutonIdx(k);
    b = allRFMapping(iROI);
    trialMatrix = b.baselineSubtracted; % cell array [nEl x nAz], GRID POSITIONS ONLY
    if isempty(trialMatrix)
        continue;
    end
    gridDims = size(trialMatrix);

    % point estimate: SAME function used for the bootstrap below, called here with the
    % untouched trial set (no resampling)
    respVec = computeGridRespVec(trialMatrix, respIdx);
    if any(isnan(respVec))
        continue; % incomplete grid coverage for this bouton
    end

    [maxVal, maxIdx] = max(respVec);
    x0_init = xyList(maxIdx, 1);
    y0_init = xyList(maxIdx, 2);
    p0 = [maxVal, x0_init, y0_init, azStep, elStep];

    try
        [pFit, resnorm] = lsqcurvefit(gaussFit2D, p0, xyList, respVec, lb, ub, fitOpts);
    catch ME
        warning('Fit failed for bouton %d: %s', iROI, ME.message);
        continue;
    end

    ssTot = sum((respVec - mean(respVec)).^2);
    R2    = 1 - resnorm / ssTot;
    allR2(k) = R2;

    if R2 < r2Floor
        continue; % excludes only catastrophic non-fits -- NOT the trust decision
    end

    %% degenerate-sigma check: is sigma >= a FRACTION of the ACTUAL tested span on that axis?
    % Coverage requirement, from Gaussian geometry (NOT tuned per axis, NOT set post-hoc):
    % resolving BOTH flanks of a Gaussian -- i.e. seeing it turn over on both sides of its
    % peak, rather than a monotonic edge-going ramp -- requires roughly 2-3 sigma-widths of
    % tested span (see header note #3). turnoverWidths is fixed ONCE here at the midpoint of
    % that range and applied IDENTICALLY to both axes; the fraction is a property of the
    % Gaussian shape, not of how many grid points either axis happens to have.
    turnoverWidths = 2.5; % 2 sigma = bare minimum to see both flanks; 3 sigma = comfortable margin
    degenerateSigmaFrac = 1 / turnoverWidths; % ~0.40, SAME value for az and el
    isDegenerateSigmaAz = pFit(4) >= degenerateSigmaFrac * range(azRange);
    isDegenerateSigmaEl = pFit(5) >= degenerateSigmaFrac * range(elRange);
    isDegenerateSigma   = isDegenerateSigmaAz || isDegenerateSigmaEl; % too large on either axis

    % Sampling-density diagnostic (REPORTING ONLY -- does NOT gate isDegenerateSigma/isTrusted).
    % This is where axis-specific grid sparsity (e.g. few elevation positions) belongs: as a
    % visible per-fit caveat, rather than folded silently into a looser coverage threshold for
    % that axis. Counts how many grid points actually fall within the resolving window
    % (+/- turnoverWidths/2 * sigma) around the fitted peak on each axis.
    minPointsInWindow = 3; % need at least 3 points to distinguish a curve from a line
    winAz = turnoverWidths/2 * pFit(4);
    winEl = turnoverWidths/2 * pFit(5);
    nPointsInWindowAz = sum(uAz(:)      >= pFit(2)-winAz & uAz(:)      <= pFit(2)+winAz);
    nPointsInWindowEl = sum(uEl_plot(:) >= pFit(3)-winEl & uEl_plot(:) <= pFit(3)+winEl);
    isSparseWindowAz  = nPointsInWindowAz < minPointsInWindow;
    isSparseWindowEl  = nPointsInWindowEl < minPointsInWindow;
    
    %% isBeyondRange / isAtBoundary: two DIFFERENT questions about the fitted center
    isBeyondRangeAz = (pFit(2) < azRangeExt(1)) || (pFit(2) > azRangeExt(2));
    isBeyondRangeEl = (pFit(3) < elRangeExt(1)) || (pFit(3) > elRangeExt(2));
    isBeyondRange   = isBeyondRangeAz || isBeyondRangeEl;

    isAtBoundaryAz = ~isBeyondRangeAz && ...
                     ((pFit(2) <= azRange(1) + edgeTolFrac*azStep) || (pFit(2) >= azRange(2) - edgeTolFrac*azStep));
    isAtBoundaryEl = ~isBeyondRangeEl && ...
                     ((pFit(3) <= elRange(1) + edgeTolFrac*elStep) || (pFit(3) >= elRange(2) - edgeTolFrac*elStep));
    isAtBoundary   = isAtBoundaryAz || isAtBoundaryEl;

    %% bootstrap R^2 distribution (evaluated post-hoc): trial-resampling WITH REPLACEMENT,
    % WITHIN each position's correct assignment (positions are NOT scrambled). The full
    % distribution of nBootstrap R^2 values is stored raw -- no threshold or CI level applied
    % at fit time. Robustness at any given trustR2Threshold/ciLevel is computed afterward.
    bootR2 = nan(nBootstrap, 1);
    for bIter = 1:nBootstrap
        resampledTrialMatrix = cell(gridDims);
        for posIdx = 1:numel(trialMatrix)
            trials = trialMatrix{posIdx};
            if isempty(trials)
                continue;
            end
            trials = double(trials);
            nTrialsHere = size(trials, 1);
            if nTrialsHere == 0
                continue;
            end
            sampIdx = randi(nTrialsHere, nTrialsHere, 1);
            resampledTrialMatrix{posIdx} = trials(sampIdx, :);
        end
        bootRespVec = computeGridRespVec(resampledTrialMatrix, respIdx);
        if any(isnan(bootRespVec))
            continue;
        end
        try
            [~, resnormBoot] = lsqcurvefit(gaussFit2D, pFit, xyList, bootRespVec, lb, ub, fitOpts);
            ssTotBoot = sum((bootRespVec - mean(bootRespVec)).^2);
            bootR2(bIter) = 1 - resnormBoot / ssTotBoot;
        catch
            continue;
        end
    end

    nBootValid = sum(~isnan(bootR2));

    gaussFitResults(end+1) = struct('iROI', iROI, 'isResponsive', b.isResponsive, ...
        'A', pFit(1), 'x0', pFit(2), 'y0', pFit(3), 'sigmaX', pFit(4), 'sigmaY', pFit(5), 'R2', R2, ...
        'isDegenerateSigmaAz', isDegenerateSigmaAz, 'isDegenerateSigmaEl', isDegenerateSigmaEl, 'isDegenerateSigma', isDegenerateSigma, ...
        'nPointsInWindowAz', nPointsInWindowAz, 'nPointsInWindowEl', nPointsInWindowEl, ...
        'isSparseWindowAz', isSparseWindowAz, 'isSparseWindowEl', isSparseWindowEl, ...
        'isBeyondRangeAz', isBeyondRangeAz, 'isBeyondRangeEl', isBeyondRangeEl, 'isBeyondRange', isBeyondRange, ...
        'isAtBoundaryAz', isAtBoundaryAz, 'isAtBoundaryEl', isAtBoundaryEl, 'isAtBoundary', isAtBoundary, ...
        'bootR2', bootR2, 'nBootValid', nBootValid); 

    allRFMapping(iROI).gaussFit_A       = pFit(1);
    allRFMapping(iROI).gaussFit_Az0     = pFit(2);
    allRFMapping(iROI).gaussFit_El0     = pFit(3);
    allRFMapping(iROI).gaussFit_sigmaAz = pFit(4);
    allRFMapping(iROI).gaussFit_sigmaEl = pFit(5);
    allRFMapping(iROI).gaussFit_R2      = R2;
    allRFMapping(iROI).gaussFit_isDegenerateSigmaAz = isDegenerateSigmaAz;
    allRFMapping(iROI).gaussFit_isDegenerateSigmaEl = isDegenerateSigmaEl;
    allRFMapping(iROI).gaussFit_isDegenerateSigma   = isDegenerateSigma;
    allRFMapping(iROI).gaussFit_nPointsInWindowAz = nPointsInWindowAz; % reporting only, not a gate
    allRFMapping(iROI).gaussFit_nPointsInWindowEl = nPointsInWindowEl; % reporting only, not a gate
    allRFMapping(iROI).gaussFit_isSparseWindowAz  = isSparseWindowAz;  % reporting only, not a gate
    allRFMapping(iROI).gaussFit_isSparseWindowEl  = isSparseWindowEl;  % reporting only, not a gate
    allRFMapping(iROI).gaussFit_isBeyondRangeAz = isBeyondRangeAz;
    allRFMapping(iROI).gaussFit_isBeyondRangeEl = isBeyondRangeEl;
    allRFMapping(iROI).gaussFit_isBeyondRange   = isBeyondRange;
    allRFMapping(iROI).gaussFit_isAtBoundaryAz  = isAtBoundaryAz; % reporting only, not a gate
    allRFMapping(iROI).gaussFit_isAtBoundaryEl  = isAtBoundaryEl; % reporting only, not a gate
    allRFMapping(iROI).gaussFit_isAtBoundary    = isAtBoundary;   % reporting only, not a gate
    % NOTE: gaussFit_isTrusted / gaussFit_isRobust / gaussFit_bootR2LowerCI are NOT set here --
    % they depend on trustR2Threshold/ciLevel, chosen in the commit block below.
end


nAboveFloor = numel(gaussFitResults);
fprintf('\n%d / %d ROIs passed the screening floor R^2 >= %.2f.\n', nAboveFloor, numBoutons, r2Floor);
fprintf('No isTrusted decision has been made yet -- see the commit block below.\n');

%% R^2 threshold x bootstrap-CI sensitivity sweep (uses stored bootR2 -- NO REFITTING)
R2vals      = [gaussFitResults.R2];
sigAzAll    = [gaussFitResults.sigmaX];
sigElAll    = [gaussFitResults.sigmaY];
isDegenAll  = [gaussFitResults.isDegenerateSigma];
isBeyondAll = [gaussFitResults.isBeyondRange];
isAtBoundAll = [gaussFitResults.isAtBoundary];

% Wall-clipping check: was the fitted center pinned against the lsqcurvefit box constraint
% (lb/ub), rather than landing at a genuine interior optimum? If so, the "true" unconstrained
% center/overshoot is unknown -- this is what we gate on, NOT mere isBeyondRange (which, per the
% header notes, only reflects centers landing ~5-20 deg past the extended range -- about one
% grid step -- consistent with normal center-estimate imprecision, not genuine extrapolation).
tol = 1e-6;
x0All = [gaussFitResults.x0];
y0All = [gaussFitResults.y0];
isWallClippedAz = (abs(x0All - ub(2)) < tol) | (abs(x0All - lb(2)) < tol);
isWallClippedEl = (abs(y0All - ub(3)) < tol) | (abs(y0All - lb(3)) < tol);
isWallClipped   = isWallClippedAz | isWallClippedEl;

ciLevel = 0.95; % one-sided confidence level for the bootstrap lower bound
lowerPct = 100 * (1 - ciLevel);
fprintf('\n=== R^2 threshold x %.0f%% bootstrap-CI x non-degenerate-sigma x non-wall-clipped sensitivity sweep (n = %d fitted-above-screening-floor ROIs) ===\n', 100*ciLevel, numBoutons);
threshList = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
for th = threshList
    passR2 = R2vals >= th;
    isRobustHere = false(size(gaussFitResults));
    for gi = 1:numel(gaussFitResults)
        validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
        if numel(validBoot) >= 800
            lowerBound = prctile(validBoot, lowerPct);
            isRobustHere(gi) = lowerBound >= th;
        end
    end
    % GATED on: R2 threshold, bootstrap-CI robustness, non-degenerate-sigma, non-wall-clipped.
    % isBeyondRange itself is NOT gated (reporting only -- see header note #4 and isBeyondAll below).
    passAll  = passR2 & isRobustHere & ~isDegenAll & ~isWallClipped;
    nR2Only  = sum(passR2);
    nTrusted = sum(passAll);
    if nTrusted > 0
        sAz = mean(sigAzAll(passAll));
        sEl = mean(sigElAll(passAll));
        sdAz = std(sigAzAll(passAll));
        sdEl = std(sigElAll(passAll));
    else
        sAz = NaN; sEl = NaN; sdAz = NaN; sdEl = NaN;
    end
    fprintf('  R^2 >= %.1f : %d pass R^2 alone | %d pass R^2 + %.0f%% bootstrap-CI + non-degenerate-sigma + non-wall-clipped (%.1f%% of R^2-passing) | sigAz=%.2f+/-%.2f sigEl=%.2f+/-%.2f\n', ...
        th, nR2Only, nTrusted, 100*ciLevel, 100*nTrusted/max(nR2Only,1), sAz, sdAz, sEl, sdEl);
end
fprintf('  For reference, degenerate-sigma alone excludes %d / %d above-floor fits (sigma >= %.0f%% of tested span on at least one axis, or peak-at-edge for the coarse elevation axis).\n', ...
    sum(isDegenAll), numel(gaussFitResults), 100*degenerateSigmaFrac);
fprintf('  Beyond-range (reporting only, NOT gated): %d / %d above-floor fits have a center outside the stimulus-footprint-extended range.\n', ...
    sum(isBeyondAll), numel(gaussFitResults));
fprintf('  Wall-clipped (GATED): %d / %d above-floor fits have a center pinned against the fitting bound (lb/ub) -- true optimum unknown.\n', ...
    sum(isWallClipped), numel(gaussFitResults));
fprintf('  At-boundary (KEPT, reporting only): %d / %d fits have a center at/near the outermost tested position but within range.\n', ...
    sum(isAtBoundAll), numel(gaussFitResults));

%% COMMIT: final inclusion criteria
% trustR2Threshold and robustCIThreshold are MATCHED (both set to the same value).
% NOTE: this commit block now gates on isWallClipped, not isBeyondRange -- see header note #4
% and the sweep section above for why (isBeyondRange alone was excluding legitimate
% boundary-adjacent fits whose overshoot was within ~1 grid step of the extended range; only
% fits actually pinned against the optimizer's lb/ub -- where the true center is unknown --
% are excluded here).
trustR2Threshold  = 0.2;
robustCIThreshold = 0.2;

for gi = 1:numel(gaussFitResults)
    iROI = gaussFitResults(gi).iROI;
    validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));

    if numel(validBoot) >= 800
        lowerBound = prctile(validBoot, lowerPct);
    else
        lowerBound = NaN;
    end

    isRobust = ~isnan(lowerBound) && (lowerBound >= robustCIThreshold);

    r2here = gaussFitResults(gi).R2;
    if isempty(r2here), r2here = -Inf; end

    degenFlag = gaussFitResults(gi).isDegenerateSigma;
    if isempty(degenFlag), degenFlag = false; end
    degenFlag = logical(degenFlag(1));

    isTrusted = (r2here >= trustR2Threshold) && isRobust && ~degenFlag && ~isWallClipped(gi);

    gaussFitResults(gi).bootR2LowerCI = lowerBound;
    gaussFitResults(gi).isRobust      = isRobust;
    gaussFitResults(gi).isTrusted     = isTrusted;
    gaussFitResults(gi).isWallClipped = isWallClipped(gi); % stored for inspection/plotting

    allRFMapping(iROI).gaussFit_bootR2LowerCI = lowerBound;
    allRFMapping(iROI).gaussFit_isRobust      = isRobust;
    allRFMapping(iROI).gaussFit_isTrusted     = isTrusted;
    allRFMapping(iROI).gaussFit_isWallClipped = isWallClipped(gi);
end

isTrustedAll = [gaussFitResults.isTrusted];
fprintf('\nCOMMITTED: trustR2Threshold=%.2f, robustCIThreshold=%.2f (matched), ciLevel=%.0f%% -> %d / %d ROIs trusted.\n', ...
    trustR2Threshold, robustCIThreshold, 100*ciLevel, sum(isTrustedAll), nAboveFloor);

nResponsiveTrusted = sum(isTrustedAll & [gaussFitResults.isResponsive]);
fprintf('%d / %d of the trusted set are also isResponsive -- remember to AND with isResponsive\n', ...
    nResponsiveTrusted, sum(isTrustedAll));
fprintf('downstream, since this script fits ALL ROIs, not just responsive ones.\n');

trustedIdx = [gaussFitResults.isTrusted];
azVals = [gaussFitResults(trustedIdx).sigmaX];
elVals = [gaussFitResults(trustedIdx).sigmaY];
fprintf('n = %d\n', sum(trustedIdx));
fprintf('sigAz = %.2f +/- %.2f (SD), median = %.2f, IQR = [%.2f, %.2f]\n', ...
    mean(azVals), std(azVals), median(azVals), prctile(azVals,25), prctile(azVals,75));
fprintf('sigEl = %.2f +/- %.2f (SD), median = %.2f, IQR = [%.2f, %.2f]\n', ...
    mean(elVals), std(elVals), median(elVals), prctile(elVals,25), prctile(elVals,75));


%% plot and view
% ---- SELECT WHICH SET OF FITS TO VIEW HERE ----
% Exactly one of the blocks below should be uncommented at a time. All of them populate
% `trustedIdxList` (indices into gaussFitResults) and `numFitted`, which the viewer loop below
% uses generically -- so the loop itself never needs to change, only this selection block.

% Plot limits: two DIFFERENT ranges for two DIFFERENT panels.
% dataLims = tight grid range -- used for the raw-data heatmap (ax1), since that panel only
% ever has real measured tiles within the grid itself; extending past this leaves a blank/white
% gap where no data exists.
% modelLims = stimulus-extended, screen-clipped range (azRangeExt/elRangeExt) -- used for the
% fitted-model heatmap (ax2) only, so a fitted peak near the tested edge isn't visually clipped
% against the plot wall. This does NOT affect isDegenerateSigma (uses azRange/elRange) or
% isBeyondRange (already uses azRangeExt/elRangeExt) -- purely a display choice for ax2.
dataLims.az = [min(uAz) - azStep/2, max(uAz) + azStep/2];
dataLims.el = [min(uEl_plot) - elStep/2, max(uEl_plot) + elStep/2];
modelLims.az = azRangeExt;
modelLims.el = elRangeExt;

% --- Option 1: the final COMMITTED trusted set (default) ---
trustedIdxList = find([gaussFitResults.isTrusted]);

% --- Option 2: fits flagged degenerate on elevation but NOT on azimuth (inspect for
% peak-inside-grid vs monotonic/edge-going elevation profiles) ---
% trustedIdxList = find([gaussFitResults.isDegenerateSigmaEl] & ~[gaussFitResults.isDegenerateSigmaAz]);

% --- Option 3: fits beyond the extended range but NOT wall-clipped (the "recovered" set --
% overshoot due to normal center-estimate imprecision, not genuine extrapolation) ---
% isDegenAzOnly_forPlot = [gaussFitResults.isDegenerateSigmaAz];
% passR2_02 = R2vals >= 0.2;
% isRobustHere_02 = false(size(gaussFitResults));
% for gi = 1:numel(gaussFitResults)
%     validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
%     if numel(validBoot) >= 800
%         isRobustHere_02(gi) = prctile(validBoot, lowerPct) >= 0.2;
%     end
% end
% passAll_current_forPlot = passR2_02 & isRobustHere_02 & ~isDegenAll & ~isWallClipped;
% passAll_azOnly_forPlot  = passR2_02 & isRobustHere_02 & ~isDegenAzOnly_forPlot & ~isWallClipped;
% trustedIdxList = find(passAll_azOnly_forPlot & ~passAll_current_forPlot);

numFitted = numel(trustedIdxList);
fprintf('\nViewer: showing %d fits (see "plot and view" section to change selection).\n', numFitted);

hFig = figure('Position', [50, 100, 1250, 400], 'Name', 'RF Comprehensive Diagnostic Viewer', 'NumberTitle', 'off');

for gi_i = 1:numFitted
    if ~ishandle(hFig)
        break;
    end

    gi = trustedIdxList(gi_i);   % map back to the real index in gaussFitResults
    res = gaussFitResults(gi);
    b = allRFMapping(res.iROI);

    % Reconstruct raw response grid
    trialMatrix = b.baselineSubtracted;
    respVec = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        trials = trialMatrix{posIdx};
        if ~isempty(trials)
            trials = double(trials);
            respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
        end
    end
    respGrid = reshape(respVec, size(AzGrid));

    % Generate high-resolution model fit grid
    [xq, yq] = meshgrid(linspace(modelLims.az(1), modelLims.az(2), 100), linspace(modelLims.el(1), modelLims.el(2), 100));
    pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
    fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));

    % Prepare boutonData struct dynamically for the trace-plotting function
    boutonData = b;
    boutonData.meanGridResponse = respGrid;
    boutonData.meanTemporalResponse = b.meanTemporalResponse;
    boutonData.meanBlankResponse = b.meanBlankResponse;
    boutonData.peakAmplitude = max(respGrid(:));

    clf(hFig);

    sgtitle(sprintf('ROI %d (%d of %d shown) | Raw R^2 = %.2f | Boot Lower Bound = %.2f | Trusted = %d | WallClipped = %d | Degen = %d', ...
    res.iROI, gi_i, numFitted, res.R2, res.bootR2LowerCI, res.isTrusted, isWallClipped(gi), res.isDegenerateSigma), 'FontWeight', 'bold');

    % Subplot 1: Trace-Overlaid Spatial Heatmap -- uses dataLims (tight grid range, no white gap)
    ax1 = subplot(1, 3, 1);
    plotRFHeatmapWithTraces(ax1, boutonData, uAz, uEl_plot, timeVector, ...
        'Colormap', 'bone', 'Smooth', false);
    xlim(ax1, dataLims.az);
    ylim(ax1, dataLims.el);
    set(ax1, 'YDir', 'normal');
    title(ax1, 'Mean Grid Response & Traces');

    % Subplot 2: Fitted 2D Gaussian Model -- uses modelLims (extended range, avoids clipping
    % boundary-adjacent peaks against the plot wall)
    ax2 = subplot(1, 3, 2);
    imagesc(linspace(modelLims.az(1), modelLims.az(2), 100), linspace(modelLims.el(1), modelLims.el(2), 100), fittedGrid);
    axis image;
    colormap(ax2, 'bone');
    colorbar;
    xlim(ax2, modelLims.az);
    ylim(ax2, modelLims.el);
    set(ax2, 'YDir', 'normal');
    hold on;
    plot(ax2, res.x0, res.y0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    t_ellipse = linspace(0, 2*pi, 100);
    plot(ax2, res.x0 + res.sigmaX * cos(t_ellipse), res.y0 + res.sigmaY * sin(t_ellipse), 'r--', 'LineWidth', 1.5);
    title(ax2, 'Fitted 2D Gaussian Model');
    xlabel(ax2, 'Azimuth (deg)');
    ylabel(ax2, 'Elevation (deg)');

    % Subplot 3: Bootstrap R^2 Distribution
    ax3 = subplot(1, 3, 3);
    validBoot = res.bootR2(~isnan(res.bootR2));
    histogram(ax3, validBoot, linspace(-0.2, 1, 40), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on;
    yl = ylim(ax3);
    line([res.R2, res.R2], yl, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-');
    line([res.bootR2LowerCI, res.bootR2LowerCI], yl, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    line([trustR2Threshold, trustR2Threshold], yl, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
    legend(ax3, 'Bootstrap R^2', 'Raw R^2', sprintf('%.0fth Pct Lower Bound', lowerPct), sprintf('Threshold (%.1f)', trustR2Threshold), 'Location', 'northwest');
    title(ax3, sprintf('%d Bootstrap R^2 Distribution', nBootstrap));
    xlabel(ax3, 'R^2');
    ylabel(ax3, 'Count');
    grid on;

    drawnow;

    if gi_i < numFitted
        waitforbuttonpress;
    end
end

%% local function
function respVec = computeGridRespVec(trialMatrix, respIdx)
% Computes the per-grid-position mean response, restricted to respIdx, from a cell array of
% [Trials x Time] matrices. Used IDENTICALLY for the point estimate (full trial set) and every
% bootstrap resample (resampled trial set), so both are guaranteed to use the same operation and
% the same response window.
respVec = nan(numel(trialMatrix), 1);
for posIdx = 1:numel(trialMatrix)
    trials = trialMatrix{posIdx};
    if isempty(trials)
        continue;
    end
    trials = double(trials);
    respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
end
end


% nAboveFloor = numel(gaussFitResults);
% fprintf('\n%d / %d ROIs passed the screening floor R^2 >= %.2f.\n', nAboveFloor, numBoutons, r2Floor);
% fprintf('No isTrusted decision has been made yet -- see the commit block below.\n');
% 
% %% R^2 threshold x bootstrap-CI sensitivity sweep (uses stored bootR2 -- NO REFITTING)
% R2vals      = [gaussFitResults.R2];
% sigAzAll    = [gaussFitResults.sigmaX];
% sigElAll    = [gaussFitResults.sigmaY];
% isDegenAll  = [gaussFitResults.isDegenerateSigma];
% isBeyondAll = [gaussFitResults.isBeyondRange];
% isAtBoundAll = [gaussFitResults.isAtBoundary];
% 
% tol = 1e-6;
% x0All = [gaussFitResults.x0];
% y0All = [gaussFitResults.y0];
% isWallClippedAz = (abs(x0All - ub(2)) < tol) | (abs(x0All - lb(2)) < tol);
% isWallClippedEl = (abs(y0All - ub(3)) < tol) | (abs(y0All - lb(3)) < tol);
% isWallClipped   = isWallClippedAz | isWallClippedEl;
% 
% ciLevel = 0.95; % one-sided confidence level for the bootstrap lower bound
% lowerPct = 100 * (1 - ciLevel);
% fprintf('\n=== R^2 threshold x %.0f%% bootstrap-CI x non-degenerate-sigma x in-range sensitivity sweep (n = %d fitted-above-screening-floor ROIs) ===\n', 100*ciLevel, numBoutons);
% threshList = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
% for th = threshList
%     passR2 = R2vals >= th;
%     isRobustHere = false(size(gaussFitResults));
%     for gi = 1:numel(gaussFitResults)
%         validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
%         if numel(validBoot) >= 800
%             lowerBound = prctile(validBoot, lowerPct);
%             isRobustHere(gi) = lowerBound >= th;
%         end
%     end
%     %passAll  = passR2 & isRobustHere & ~isDegenAll & ~isBeyondAll; % isBeyondRange now gated -- center must be on-screen
%     passAll  = passR2 & isRobustHere & ~isDegenAll & ~isWallClipped;
%     nR2Only  = sum(passR2);
%     nTrusted = sum(passAll);
%     if nTrusted > 0
%         sAz = mean(sigAzAll(passAll));
%         sEl = mean(sigElAll(passAll));
%     else
%         sAz = NaN; sEl = NaN;
%     end
%     fprintf('  R^2 >= %.1f : %d pass R^2 alone | %d pass R^2 + %.0f%% bootstrap-CI + non-degenerate-sigma + in-range (%.1f%% of R^2-passing) | sigAz=%.2f sigEl=%.2f\n', ...
%         th, nR2Only, nTrusted, 100*ciLevel, 100*nTrusted/max(nR2Only,1), sAz, sEl);
% end
% fprintf('  For reference, degenerate-sigma alone excludes %d / %d above-floor fits (sigma >= %.0f%% of tested span on at least one axis).\n', ...
%     sum(isDegenAll), numel(gaussFitResults), 100*degenerateSigmaFrac);
% fprintf('  Beyond-range (GATED): %d / %d above-floor fits have a center outside the stimulus-footprint-extended range.\n', ...
%     sum(isBeyondAll), numel(gaussFitResults));
% fprintf('  At-boundary (KEPT, reporting only): %d / %d fits have a center at/near the outermost tested position but within range.\n', ...
%     sum(isAtBoundAll), numel(gaussFitResults));
% %% COMMIT: final inclusion criteria 
% %trustR2Threshold and robustCIThreshold are MATCHED (both set to the same value) -- this is
% % what was actually validated (n=82, sigAz=15.34+/-9.17, sigEl=8.59+/-3.54). 
% trustR2Threshold  = 0.2;
% robustCIThreshold = 0.2;
% 
% for gi = 1:numel(gaussFitResults)
%     iROI = gaussFitResults(gi).iROI;
%     validBoot = gaussFitResults(gi).bootR2(~isnan(gaussFitResults(gi).bootR2));
% 
%     if numel(validBoot) >= 800
%         lowerBound = prctile(validBoot, lowerPct);
%     else
%         lowerBound = NaN;
%     end
% 
%     isRobust = ~isnan(lowerBound) && (lowerBound >= robustCIThreshold);
% 
%     r2here = gaussFitResults(gi).R2;
%     if isempty(r2here), r2here = -Inf; end
% 
%     degenFlag = gaussFitResults(gi).isDegenerateSigma;
%     if isempty(degenFlag), degenFlag = false; end
%     degenFlag = logical(degenFlag(1));
% %     isTrusted = (lowerBound >= trustR2Threshold) && isRobust && ~degenFlag && ~gaussFitResults(gi).isBeyondRange; 
%     isTrusted = (r2here >= trustR2Threshold) && isRobust && ~degenFlag && ~gaussFitResults(gi).isBeyondRange;
% 
%     gaussFitResults(gi).bootR2LowerCI = lowerBound;
%     gaussFitResults(gi).isRobust      = isRobust;
%     gaussFitResults(gi).isTrusted     = isTrusted;
% 
%     allRFMapping(iROI).gaussFit_bootR2LowerCI = lowerBound;
%     allRFMapping(iROI).gaussFit_isRobust      = isRobust;
%     allRFMapping(iROI).gaussFit_isTrusted     = isTrusted;
% end
% 
% fprintf('n = %d\n', sum([gaussFitResults.isTrusted]));
% 
% isTrustedAll = [gaussFitResults.isTrusted];
% fprintf('\nCOMMITTED: trustR2Threshold=%.2f, robustCIThreshold=%.2f (matched), ciLevel=%.0f%% -> %d / %d ROIs trusted.\n', ...
%     trustR2Threshold, robustCIThreshold, 100*ciLevel, sum(isTrustedAll), nAboveFloor);
% 
% nResponsiveTrusted = sum(isTrustedAll & [gaussFitResults.isResponsive]);
% fprintf('%d / %d of the trusted set are also isResponsive -- remember to AND with isResponsive\n', ...
%     nResponsiveTrusted, sum(isTrustedAll));
% fprintf('downstream, since this script fits ALL ROIs, not just responsive ones.\n');
% 
% trustedIdx = [gaussFitResults.isTrusted];
% azVals = [gaussFitResults(trustedIdx).sigmaX];
% elVals = [gaussFitResults(trustedIdx).sigmaY];
% fprintf('n = %d\n', sum(trustedIdx));
% fprintf('sigAz = %.2f +/- %.2f (SD), median = %.2f, IQR = [%.2f, %.2f]\n', ...
%     mean(azVals), std(azVals), median(azVals), prctile(azVals,25), prctile(azVals,75));
% fprintf('sigEl = %.2f +/- %.2f (SD), median = %.2f, IQR = [%.2f, %.2f]\n', ...
%     mean(elVals), std(elVals), median(elVals), prctile(elVals,25), prctile(elVals,75));
% 
% 
% %% plot and view
% % azLims = [min(uAz) - azStep/2, max(uAz) + azStep/2];
% % elLims = [min(uEl_plot) - elStep/2, max(uEl_plot) + elStep/2];
% % 
% % trustedIdxList = find([gaussFitResults.isTrusted]); % trusted-fit indices only
% % numFitted = numel(trustedIdxList);
% 
% 
% % degenElIdx = find([gaussFitResults.isDegenerateSigmaEl] & ~[gaussFitResults.isDegenerateSigmaAz]);
% % numFitted = numel(degenElIdx);
% % trustedIdxList = degenElIdx; % feed straight into your existing plotting loop
% hFig = figure('Position', [50, 100, 1250, 400], 'Name', 'RF Comprehensive Diagnostic Viewer (Trusted Only)', 'NumberTitle', 'off');
% 
% for gi_i = 1:numFitted
%     if ~ishandle(hFig)
%         break;
%     end
% 
%     gi = trustedIdxList(gi_i);   % map back to the real index in gaussFitResults
%     res = gaussFitResults(gi);
%     b = allRFMapping(res.iROI);
% 
%     % Reconstruct raw response grid
%     trialMatrix = b.baselineSubtracted;
%     respVec = nan(numel(trialMatrix), 1);
%     for posIdx = 1:numel(trialMatrix)
%         trials = trialMatrix{posIdx};
%         if ~isempty(trials)
%             trials = double(trials);
%             respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
%         end
%     end
%     respGrid = reshape(respVec, size(AzGrid));
% 
%     % Generate high-resolution model fit grid
%     [xq, yq] = meshgrid(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100));
%     pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
%     fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));
% 
%     % Prepare boutonData struct dynamically for the trace-plotting function
%     boutonData = b;
%     boutonData.meanGridResponse = respGrid;
%     boutonData.meanTemporalResponse = b.meanTemporalResponse;
%     boutonData.meanBlankResponse = b.meanBlankResponse;
%     boutonData.peakAmplitude = max(respGrid(:));
% 
%     clf(hFig);
% 
%     sgtitle(sprintf('ROI %d (%d of %d TRUSTED) | Raw R^2 = %.2f | Boot Lower Bound = %.2f | Trusted = %d', ...
%     res.iROI, gi_i, numFitted, res.R2, res.bootR2LowerCI, res.isTrusted), 'FontWeight', 'bold');
% 
%     % Subplot 1: Trace-Overlaid Spatial Heatmap
%     ax1 = subplot(1, 3, 1);
%     plotRFHeatmapWithTraces(ax1, boutonData, uAz, uEl_plot, timeVector, ...
%         'Colormap', 'bone', 'Smooth', false);
%     xlim(ax1, azLims);
%     ylim(ax1, elLims);
%     set(ax1, 'YDir', 'normal');
%     title(ax1, 'Mean Grid Response & Traces');
% 
%     % Subplot 2: Fitted 2D Gaussian Model
%     ax2 = subplot(1, 3, 2);
%     imagesc(linspace(azLims(1), azLims(2), 100), linspace(elLims(1), elLims(2), 100), fittedGrid);
%     axis image;
%     colormap(ax2, 'bone');
%     colorbar;
%     xlim(ax2, azLims);
%     ylim(ax2, elLims);
%     set(ax2, 'YDir', 'normal');
%     hold on;
%     plot(ax2, res.x0, res.y0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
%     t_ellipse = linspace(0, 2*pi, 100);
%     plot(ax2, res.x0 + res.sigmaX * cos(t_ellipse), res.y0 + res.sigmaY * sin(t_ellipse), 'r--', 'LineWidth', 1.5);
%     title(ax2, 'Fitted 2D Gaussian Model');
%     xlabel(ax2, 'Azimuth (deg)');
%     ylabel(ax2, 'Elevation (deg)');
% 
%     % Subplot 3: Bootstrap R^2 Distribution
%     ax3 = subplot(1, 3, 3);
%     validBoot = res.bootR2(~isnan(res.bootR2));
%     histogram(ax3, validBoot, linspace(-0.2, 1, 40), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
%     hold on;
%     yl = ylim(ax3);
%     line([res.R2, res.R2], yl, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '-');
%     line([res.bootR2LowerCI, res.bootR2LowerCI], yl, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
%     line([trustR2Threshold, trustR2Threshold], yl, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':');
%     legend(ax3, 'Bootstrap R^2', 'Raw R^2', sprintf('%.0fth Pct Lower Bound', lowerPct), sprintf('Threshold (%.1f)', trustR2Threshold), 'Location', 'northwest');
%     title(ax3, sprintf('%d Bootstrap R^2 Distribution', nBootstrap));
%     xlabel(ax3, 'R^2');
%     ylabel(ax3, 'Count');
%     grid on;
% 
%     drawnow;
% 
%     if gi_i < numFitted
%         waitforbuttonpress;
%     end
% end
% 
% %% local function
% function respVec = computeGridRespVec(trialMatrix, respIdx)
% % Computes the per-grid-position mean response, restricted to respIdx, from a cell array of
% % [Trials x Time] matrices. Used IDENTICALLY for the point estimate (full trial set) and every
% % bootstrap resample (resampled trial set), so both are guaranteed to use the same operation and
% % the same response window.
% respVec = nan(numel(trialMatrix), 1);
% for posIdx = 1:numel(trialMatrix)
%     trials = trialMatrix{posIdx};
%     if isempty(trials)
%         continue;
%     end
%     trials = double(trials);
%     respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
% end
% end