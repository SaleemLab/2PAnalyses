% RFMapping_Gaussian2DFit.m
%
% Fits an anisotropic 2D Gaussian (no rotation, no baseline term) to each
% RESPONSIVE bouton's mean response across the RF mapping grid:
%
%   G(x,y) = A * exp( -( (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2) ) )
%
% QC APPROACH: R^2 alone cannot detect a key failure mode in a sparse/
% broad-RF-relative-to-grid design like this one -- if the true RF is as
% wide as (or wider than) the sampled range on an axis, the response never
% turns over, and a rising/falling partial slope can look like a good
% Gaussian fit (high R^2) for many different, mutually incompatible
% (center, sigma) combinations. So R^2 here is used only as a permissive
% floor to exclude catastrophic non-fits, NOT as the primary reliability
% criterion. Reliability is instead assessed via 500-resample bootstrapping
% (trials resampled with replacement, per grid position, refit each time)
% for EVERY fit above the floor, tracking:
%   - stability of the fitted CENTER (x0, y0) across resamples
%   - stability of the fitted WIDTH (sigmaX, sigmaY) across resamples
%   - a DIRECT degenerate-sigma check (sigma >= the actual tested span,
%     not the loose fitting bound) -- catches flat/noise fits that a
%     loose "sigma pinned at bound" check would miss entirely
% Fitted centers are also checked against the tested grid range, but split
% into two DIFFERENT categories rather than one blanket "edge" flag:
%   - isBeyondRange: center is genuinely OUTSIDE the min/max tested values
%     on that axis -- this is true extrapolation and is excluded from the
%     trusted set.
%   - isAtBoundary: center sits AT/NEAR the outermost tested position, but
%     still WITHIN range -- this is NOT excluded, since it may reflect a
%     real preference for the most peripheral position actually tested
%     (confirmed necessary after inspecting the preferred-position
%     histograms, which showed a large, likely genuine cluster of boutons
%     preferring the most lateral/uppermost tested positions).
% A combined R^2-floor x bootstrap-stability sensitivity table is printed
% so the final inclusion criteria can be chosen based on how your actual
% data behaves, rather than importing a single number from a paper with a
% different grid density (7x3=21 positions there vs. 16 here).
%
% REQUIRES (already in workspace from the main pooling script):
%   allRFMapping   - pooled bouton struct array, with .meanGridResponse
%                    ([nEl x nAz] matrix) and .baselineSubtracted
%                    ([nEl x nAz] cell of [Trials x Time] matrices)
%   uAz, uEl_plot  - grid position vectors (uEl_plot rows, uAz columns,
%                    matching the orientation of meanGridResponse)
%   respIdx        - response-window time mask (from main script)
%
% REQUIRES: Optimization Toolbox (lsqcurvefit).
%
% COST WARNING: bootstrapping now runs for every responsive bouton above
% the R^2 floor (not just edge-adjacent ones), so this is meaningfully
% slower than the previous version. Consider nBootstrap = 100 for a first
% pass before committing to the full 500.

%% params
rng(1, 'twister'); % FIXED SEED -- without this, every run draws different bootstrap resamples,
                    % which can change WHICH boutons pass the coverage criterion. This was
                    % discovered directly: n went from 87->62 trusted and the sigmaAz-vs-sigmaEl
                    % asymmetry flipped from p<0.0001 to p=0.15 between two runs of this exact
                    % script on the same data, purely due to random resampling differences.
r2Floor      = 0.2;   % primary filter, closer to the paper's approach -- NOT just a permissive floor
                       % anymore. Bootstrap stability is no longer a hard gate (see note below);
                       % R^2 is doing the main work here, same as Timplalexi et al., just tuned to
                       % your grid via the earlier sensitivity check rather than importing 0.7 blindly.
edgeTolFrac       = 0.5;  % fraction of one grid-step counted as "edge-adjacent" (reporting only)
coverageThreshold = 0.95; % matches the published convention EXACTLY (95% within one grid-step).
                          % Empirically verified on this dataset: all 87 boutons passing the looser
                          % 68% version also pass 95% (74/74 boundary-adjacent, 13/13 interior) --
                          % so using the literal published threshold costs nothing here and lets the
                          % methods section cite the convention directly rather than justify a deviation.
coverageWindowSteps = 1;  % literal "one step size," matching the published wording exactly
sigmaBoundTolFrac = 0.9;  % fraction of the upper sigma bound counted as "pinned" (signature of a
                          % receptive field broader than the sampled range -- never turns over)
nBootstrap        = 500;  % consider starting with e.g. 100 for a faster first pass -- this now runs
                          % for every responsive bouton, not just edge-adjacent ones

if ~exist('respIdx', 'var')
    error('respIdx not found -- run the main pooling script first.');
end

azStep  = mean(diff(sort(uAz)));
elStep  = mean(diff(sort(uEl_plot)));

% NOTE: azRange/elRange below are the raw grid CENTER coordinates. But each
% stimulus is a 30-degree diameter grating, so the position centered at
% e.g. -70 deg already visually stimulates out to roughly -70 - 15 = -85
% deg. The tested visual field is therefore wider than the array of grid
% centers implies. stimRadius is added when deciding isBeyondRange, so a
% fitted center isn't wrongly flagged as "extrapolated" when it actually
% falls within visual space the outermost stimulus already covered.
stimDiameterDeg = 30; % set to match your actual grating diameter
stimRadius = stimDiameterDeg / 2;

% Screen physically clips the stimulus footprint near the edges -- a
% grating centered at your most-negative azimuth position does NOT fully
% extend by stimRadius, because part of it is drawn off-screen. Confirmed:
%   azimuth lower screen limit  = -80 deg (clipped -- deep temporal/lateral
%                                 periphery, near the physical monitor edge)
%   azimuth upper screen limit  = no clipping (this side sits toward central/
%                                 nasal visual space, well within the monitor;
%                                 the flat stimRadius extension applies as-is)
%   elevation lower screen limit = -25 deg (clipped)
%   elevation upper screen limit = +45 deg (clipped)
screenAzLimits = [-80, NaN];   % NaN here means "no clipping on this side" -- flat extension applies
screenElLimits = [-25, 45];    % both confirmed, both clipped

azRange = [min(uAz) max(uAz)];
elRange = [min(uEl_plot) max(uEl_plot)];

azExtLow  = max(azRange(1)-stimRadius, screenAzLimits(1));
azExtHigh = azRange(2)+stimRadius; % placeholder -- confirm actual screen limit
if ~isnan(screenAzLimits(2)), azExtHigh = min(azExtHigh, screenAzLimits(2)); end

elExtLow  = elRange(1)-stimRadius; % placeholder -- confirm actual screen limit
if ~isnan(screenElLimits(1)), elExtLow = max(elExtLow, screenElLimits(1)); end
elExtHigh = elRange(2)+stimRadius; % placeholder -- confirm actual screen limit
if ~isnan(screenElLimits(2)), elExtHigh = min(elExtHigh, screenElLimits(2)); end

azRangeExt = [azExtLow, azExtHigh];
elRangeExt = [elExtLow, elExtHigh];

[AzGrid, ElGrid] = meshgrid(uAz, uEl_plot);   % [nEl x nAz], matches meanGridResponse orientation
xyList = [AzGrid(:), ElGrid(:)];

gaussFit2D = @(p, xy) p(1) * exp( -( (xy(:,1)-p(2)).^2 ./ (2*p(4)^2) + (xy(:,2)-p(3)).^2 ./ (2*p(5)^2) ) );
% p = [A, x0, y0, sigmaX, sigmaY]

fitOpts = optimoptions('lsqcurvefit', 'Display', 'off');

respBoutonIdx = find([allRFMapping.isResponsive]);
numBoutons    = numel(respBoutonIdx);

% Clear any gaussFit_* fields left over from a PREVIOUS run with different criteria
% (e.g. a different r2Floor) -- otherwise boutons that don't clear this run's floor
% keep stale field values from before, silently contaminating downstream scripts.
staleFields = {'gaussFit_A','gaussFit_Az0','gaussFit_El0','gaussFit_sigmaAz','gaussFit_sigmaEl', ...
    'gaussFit_R2','gaussFit_isBeyondRange','gaussFit_isAtBoundary','gaussFit_isDegenerateSigma', ...
    'gaussFit_isStable','gaussFit_isTrusted','gaussFit_bootCenterStdAz','gaussFit_bootCenterStdEl'};
for iROI = respBoutonIdx(:)'
    for f = staleFields
        if isfield(allRFMapping, f{1})
            allRFMapping(iROI).(f{1}) = [];
        end
    end
end

gaussFitResults = struct('iROI', {}, 'A', {}, 'x0', {}, 'y0', {}, 'sigmaX', {}, 'sigmaY', {}, 'R2', {}, ...
    'isBeyondRange', {}, 'isAtBoundary', {}, 'isDegenerateSigma', {}, 'bootCenterStdAz', {}, 'bootCenterStdEl', {}, ...
    'bootSigmaStdAz', {}, 'bootSigmaStdEl', {}, ...
    'pctWithinStepAz', {}, 'pctWithinStepEl', {}, ...
    'fracSigmaPinnedAz', {}, 'fracSigmaPinnedEl', {}, ...
    'isStable', {}, 'nBootValid', {});

fprintf('Fitting 2D Gaussian RFs for %d RESPONSIVE boutons (R^2 floor = %.2f)...\n', numBoutons, r2Floor);
fprintf('Bootstrapping ALL fits (%d resamples each) -- this may take a while.\n', nBootstrap);

allR2 = nan(numBoutons, 1);
sigmaUB = 3*max(range(azRange), range(elRange)); % matches ub(4)/ub(5) below -- kept as one variable for clarity

for k = 1:numBoutons
    iROI = respBoutonIdx(k);
    b = allRFMapping(iROI);
    meanGridResponse = b.meanGridResponse;
    if isempty(meanGridResponse)
        continue;
    end
    respVec = double(meanGridResponse(:));

    [maxVal, maxIdx] = max(respVec);
    x0_init = xyList(maxIdx, 1);
    y0_init = xyList(maxIdx, 2);

    p0 = [maxVal, x0_init, y0_init, azStep, elStep];
    lb = [0,        azRange(1)-azStep, elRange(1)-elStep, azStep/4,        elStep/4];
    ub = [Inf,       azRange(2)+azStep, elRange(2)+elStep, 3*range(azRange), 3*range(elRange)];

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
        continue; % excludes only catastrophic non-fits
    end

    % Split into two distinct cases, since they mean very different things:
    %  - isBeyondRange: fitted center is genuinely OUTSIDE the tested min/max on that axis.
    %    This is the suspicious case (extrapolation beyond what was ever measured) and should
    %    be excluded from the trusted set.
    %  - isAtBoundary: fitted center sits AT or NEAR the outermost tested position, but not
    %    beyond it. This can be a completely real preference (the bouton may genuinely prefer
    %    the most peripheral position you tested) and should NOT be excluded just for that.
    % isBeyondRange now uses the STIMULUS-FOOTPRINT-extended range: a
    % center is only flagged as genuine extrapolation if it falls outside
    % the visual space the outermost stimulus actually covered, not just
    % outside the bare grid-center coordinates.
    isBeyondRange = (pFit(2) < azRangeExt(1)) || (pFit(2) > azRangeExt(2)) || ...
                    (pFit(3) < elRangeExt(1)) || (pFit(3) > elRangeExt(2));
    % isAtBoundary still refers to the original grid-center range -- this
    % is reporting "did the fit land at/near a position you actually
    % tested," which is a different question from whether it's extrapolated.
    isAtBoundary  = ~isBeyondRange && ( ...
                    (pFit(2) <= azRange(1) + edgeTolFrac*azStep) || (pFit(2) >= azRange(2) - edgeTolFrac*azStep) || ...
                    (pFit(3) <= elRange(1) + edgeTolFrac*elStep) || (pFit(3) >= elRange(2) - edgeTolFrac*elStep) );

    % NOTE: the fitting bound on sigma (ub(4)/ub(5) = 3x the tested range) was set loose
    % deliberately, to let a genuinely broad RF be estimated rather than artificially capped.
    % But that means the earlier "fraction pinned at upper bound" check could never catch a
    % degenerate fit (sigma ~60-100 deg on a ~100 deg grid is clearly not a real RF estimate,
    % but is nowhere near 3x the range). This is a DIFFERENT, more honest check: is sigma bigger
    % than the actual tested span itself? If so, the fit cannot be distinguishing "real broad RF"
    % from "flat noise the optimizer gave up on" -- exclude it regardless of R^2 or stability.
    isDegenerateSigma = (pFit(4) >= range(azRange)) || (pFit(5) >= range(elRange));

    %% bootstrap: now runs for every fit above the floor, not just edge-adjacent ones
    bootStdAz = NaN; bootStdEl = NaN; bootSigStdAz = NaN; bootSigStdEl = NaN;
    fracPinAz = NaN; fracPinEl = NaN; pctWithinAz = NaN; pctWithinEl = NaN; nBootValid = 0; isStable = false;

    trialMatrix = b.baselineSubtracted;
    if ~isempty(trialMatrix)
        bootCenters = nan(nBootstrap, 2);
        bootSigmas  = nan(nBootstrap, 2);
        for bIter = 1:nBootstrap
            bootRespVec = nan(size(respVec));
            for posIdx = 1:numel(respVec)
                [rIdx, cIdx] = ind2sub(size(meanGridResponse), posIdx);
                trials = trialMatrix{rIdx, cIdx};
                if isempty(trials)
                    continue;
                end
                trials = double(trials);
                nTrialsHere = size(trials, 1);
                if nTrialsHere == 0
                    continue;
                end
                sampIdx = randi(nTrialsHere, nTrialsHere, 1);
                resampledTrials = trials(sampIdx, :);
                bootRespVec(posIdx) = mean(mean(resampledTrials(:, respIdx), 2, 'omitnan'), 'omitnan');
            end
            if any(isnan(bootRespVec))
                continue;
            end
            try
                pBoot = lsqcurvefit(gaussFit2D, pFit, xyList, bootRespVec, lb, ub, fitOpts);
                bootCenters(bIter, :) = [pBoot(2), pBoot(3)];
                bootSigmas(bIter, :)  = [pBoot(4), pBoot(5)];
                nBootValid = nBootValid + 1;
            catch
                continue;
            end
        end
        validRows = ~any(isnan(bootCenters), 2);
        if sum(validRows) >= 20
            bootStdAz    = std(bootCenters(validRows, 1));
            bootStdEl    = std(bootCenters(validRows, 2));
            bootSigStdAz = std(bootSigmas(validRows, 1));
            bootSigStdEl = std(bootSigmas(validRows, 2));
            % "pinned" = sigma landed near its upper bound in this resample -- the signature
            % of a receptive field broader than the sampled range (never turns over)
            fracPinAz = mean(bootSigmas(validRows, 1) >= sigmaBoundTolFrac * ub(4));
            fracPinEl = mean(bootSigmas(validRows, 2) >= sigmaBoundTolFrac * ub(5));

            % COVERAGE-BASED criterion, adapted for grid coarseness (see param note above): what
            % fraction of bootstrap-resampled centers fall within coverageWindowSteps grid-steps
            % of the ORIGINAL fit's center, per axis? Require >= coverageThreshold (95%) on BOTH axes.
            pctWithinAz = mean(abs(bootCenters(validRows, 1) - pFit(2)) <= coverageWindowSteps*azStep);
            pctWithinEl = mean(abs(bootCenters(validRows, 2) - pFit(3)) <= coverageWindowSteps*elStep);

            isStable = (pctWithinAz >= coverageThreshold) && (pctWithinEl >= coverageThreshold);
        end
    end
    % a fit can't be considered trustworthy if its sigma is degenerate, regardless of bootstrap
    % center stability -- a wide flat "hill" can still refit to a similar (degenerate) shape
    % consistently across resamples without that meaning the fit is meaningful
    isStable = isStable && ~isDegenerateSigma;

    gaussFitResults(end+1) = struct('iROI', iROI, 'A', pFit(1), 'x0', pFit(2), 'y0', pFit(3), ...
        'sigmaX', pFit(4), 'sigmaY', pFit(5), 'R2', R2, 'isBeyondRange', isBeyondRange, 'isAtBoundary', isAtBoundary, ...
        'isDegenerateSigma', isDegenerateSigma, ...
        'bootCenterStdAz', bootStdAz, 'bootCenterStdEl', bootStdEl, ...
        'bootSigmaStdAz', bootSigStdAz, 'bootSigmaStdEl', bootSigStdEl, ...
        'pctWithinStepAz', pctWithinAz, 'pctWithinStepEl', pctWithinEl, ...
        'fracSigmaPinnedAz', fracPinAz, 'fracSigmaPinnedEl', fracPinEl, ...
        'isStable', isStable, 'nBootValid', nBootValid); %#ok<SAGROW>

    allRFMapping(iROI).gaussFit_A       = pFit(1);
    allRFMapping(iROI).gaussFit_Az0     = pFit(2);
    allRFMapping(iROI).gaussFit_El0     = pFit(3);
    allRFMapping(iROI).gaussFit_sigmaAz = pFit(4);
    allRFMapping(iROI).gaussFit_sigmaEl = pFit(5);
    allRFMapping(iROI).gaussFit_R2      = R2;
    allRFMapping(iROI).gaussFit_isBeyondRange = isBeyondRange;
    allRFMapping(iROI).gaussFit_isAtBoundary  = isAtBoundary;
    allRFMapping(iROI).gaussFit_isDegenerateSigma = isDegenerateSigma;
    allRFMapping(iROI).gaussFit_isStable = isStable;
    % final combined trust flag: R^2 floor + bootstrap coverage (95% of resamples within one
    % grid-step, matching published convention) + the two visually-validated checks. All four
    % are hard gates.
    allRFMapping(iROI).gaussFit_isTrusted = isStable && ~isBeyondRange && ~isDegenerateSigma;
    allRFMapping(iROI).gaussFit_bootCenterStdAz = bootStdAz;
    allRFMapping(iROI).gaussFit_bootCenterStdEl = bootStdEl;
end

nAboveFloor = numel(gaussFitResults);
fprintf('\n%d / %d RESPONSIVE boutons passed the R^2 >= %.2f floor.\n', nAboveFloor, numBoutons, r2Floor);

isStableAll        = [gaussFitResults.isStable];   % now coverage-based, and IS the gate
isBeyondRangeAll    = [gaussFitResults.isBeyondRange];
isAtBoundaryAll     = [gaussFitResults.isAtBoundary];
isDegenAll          = [gaussFitResults.isDegenerateSigma];
isTrustedAll        = isStableAll & ~isBeyondRangeAll & ~isDegenAll;

fprintf('%d / %d have centers genuinely OUTSIDE the tested range (excluded from trusted set).\n', sum(isBeyondRangeAll), nAboveFloor);
fprintf('%d / %d have centers AT/NEAR the tested boundary but within range (KEPT -- may be real peripheral preference).\n', ...
    sum(isAtBoundaryAll), nAboveFloor);
fprintf('%d / %d have a degenerate (span-exceeding) sigma on at least one axis (excluded from trusted set).\n', sum(isDegenAll), nAboveFloor);
fprintf('%d / %d pass bootstrap coverage (>=%.0f%% of resamples within one grid-step, both axes) (excluded from trusted set if false).\n', ...
    sum(isStableAll), nAboveFloor, 100*coverageThreshold);
fprintf('%d / %d pass the trusted set (R^2 floor AND coverage-stable AND non-degenerate AND not-beyond-range).\n', ...
    sum(isTrustedAll), nAboveFloor);

fracPinAzAll = [gaussFitResults.fracSigmaPinnedAz];
fracPinElAll = [gaussFitResults.fracSigmaPinnedEl];
fprintf('\nMedian fraction of bootstrap resamples with sigma pinned at its upper bound:\n');
fprintf('  Azimuth:   %.1f%%\n', 100*median(fracPinAzAll, 'omitnan'));
fprintf('  Elevation: %.1f%%\n', 100*median(fracPinElAll, 'omitnan'));
fprintf(['NOTE: a high pinned fraction means the bootstrap keeps failing to see the response turn back down --\n' ...
    'that axis'' true RF width likely exceeds your sampled range for a meaningful fraction of boutons,\n' ...
    'regardless of R^2. Treat sigma (and possibly center) on that axis with real caution for those boutons.\n']);

%% ===================== combined R^2-floor x (stability + non-degenerate + not-beyond-range) sensitivity =====================
fprintf('\n=== Combined sensitivity: R^2 floor x [stable AND non-degenerate AND not-beyond-range] (n = %d responsive boutons) ===\n', numBoutons);
threshList = [0.3 0.4 0.5 0.6 0.7];
for th = threshList
    passR2 = allR2 >= th;
    trustedForThresh = false(numBoutons, 1);
    for gi = 1:numel(gaussFitResults)
        kMatch = find(respBoutonIdx == gaussFitResults(gi).iROI, 1);
        if ~isempty(kMatch) && gaussFitResults(gi).R2 >= th
            trustedForThresh(kMatch) = gaussFitResults(gi).isStable && ~gaussFitResults(gi).isBeyondRange && ~gaussFitResults(gi).isDegenerateSigma;
        end
    end
    nR2Only = sum(passR2);
    nTrusted = sum(passR2 & trustedForThresh);
    fprintf('  R^2 >= %.1f : %d pass R^2 alone | %d pass ALL checks (%.1f%% of R^2-passing)\n', ...
        th, nR2Only, nTrusted, 100*nTrusted/max(nR2Only,1));
end

nPassed = numel(gaussFitResults);
fprintf('\n%d / %d RESPONSIVE boutons passed the R^2 floor (%.2f).\n', nPassed, numBoutons, r2Floor);
fprintf('Use the combined-sensitivity table above to choose your actual final R^2 x stability criteria.\n');

%% ===================== coverage summary (uses stored per-bouton coverage, no re-run needed) =====================
pctWithinAzAll = [gaussFitResults.pctWithinStepAz];
pctWithinElAll = [gaussFitResults.pctWithinStepEl];
fprintf('\n=== Coverage summary (n = %d fitted boutons, R^2 floor only) ===\n', nPassed);
fprintf('Median %% of resamples within one grid-step: azimuth = %.1f%%, elevation = %.1f%%\n', ...
    100*median(pctWithinAzAll, 'omitnan'), 100*median(pctWithinElAll, 'omitnan'));
fprintf('Coverage threshold in use: %.0f%% (matches published convention); %d / %d pass this + non-degenerate + not-beyond-range.\n', ...
    100*coverageThreshold, sum(isTrustedAll), nPassed);
