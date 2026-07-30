% RFMapping_Gaussian2DFit.m
%
% Fits an anisotropic 2D Gaussian (no rotation, no baseline term) to each
% bouton's mean response across the RF mapping grid:
%
%   G(x,y) = A * exp( -( (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2) ) )
%
% NOTE ON SCOPE: this version fits ALL ROIs in allRFMapping, not just
% responsive ones (previous versions restricted fitting to
% allRFMapping.isResponsive). This means gaussFit_isTrusted alone is NOT
% a safe stand-in for "this is a real RF" -- a non-responsive bouton's
% flat/noisy response can still clear the permissive R^2 floor by chance,
% especially with only 16 grid positions. Downstream analyses should
% combine isResponsive AND gaussFit_isTrusted (or the per-axis trusted
% flags below), not gaussFit_isTrusted alone. The isResponsive field is
% left untouched by this script so that filter stays available.
%
% NOTE ON RESPONSE WINDOW CONSISTENCY: the response vector fed to the
% ORIGINAL (non-bootstrap) fit is now built directly from
% allRFMapping.baselineSubtracted, averaged over respIdx, exactly matching
% how the bootstrap below has always computed its resampled response
% vectors. Previously the original fit used allRFMapping.meanGridResponse
% as-is, which risked being computed over a DIFFERENT (e.g. full-trial,
% including baseline/ISI) window upstream -- diluting peak-vs-surround
% contrast and biasing fitted sigma downward, while the bootstrap
% resampled a properly-windowed signal around that mismatched point
% estimate. Deriving both from the same respIdx-restricted window removes
% this inconsistency at the source. baselineSubtracted contains grid
% positions ONLY (blank trials live separately in baselineSubtractedBlank
% and are never referenced by this script), so blank is correctly
% excluded from every fit.
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
%
% PER-AXIS GATING (NEW): isDegenerateSigma, isBeyondRange, isAtBoundary,
% and the final trusted flag are now computed BOTH per-axis (...Az /
% ...El) AND as a combined OR/AND across axes, kept for backward
% compatibility. Previously a bouton with a perfectly good, narrow
% azimuth RF could be excluded entirely just because elevation coverage
% happened to be too narrow to resolve width on that one axis -- a
% granularity problem, not a real quality problem. Now:
%   - gaussFit_isTrustedAz / gaussFit_isTrustedEl: trust per individual axis
%   - gaussFit_isTrusted: TRUE only if BOTH axes are trusted (equivalent
%     to the old combined behavior) -- use this for analyses that need a
%     fully-trusted 2D fit (e.g. sigmaAz vs sigmaEl comparisons, 2D center
%     plots). Use the per-axis flags for analyses that only need one axis
%     to be reliable (e.g. an azimuth-only tuning summary).
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
% COST WARNING: bootstrapping now runs for every fit above the R^2 floor
% (not just edge-adjacent ones), AND this version fits every ROI (not
% just responsive ones), so this is meaningfully slower than previous
% versions. Consider nBootstrap = 100 for a first pass before committing
% to the full 500.

%% params
rng(1, 'twister'); % FIXED SEED -- without this, every run draws different bootstrap resamples,
                    % which can change WHICH boutons pass the coverage criterion. This was
                    % discovered directly: n went from 87->62 trusted and the sigmaAz-vs-sigmaEl
                    % asymmetry flipped from p<0.0001 to p=0.15 between two runs of this exact
                    % script on the same data, purely due to random resampling differences.
r2Floor      = 0.1;   % primary filter, closer to the paper's approach -- NOT just a permissive floor
                       % anymore. Bootstrap stability is no longer a hard gate (see note below);
                       % R^2 is doing the main work here, same as Timplalexi et al., just tuned to
                       % your grid via the earlier sensitivity check rather than importing 0.7 blindly.
                       % [NOTE: this comment is stale relative to the actual code below -- bootstrap
                       % coverage (isStable) IS still a hard gate on the trusted flags. Kept here
                       % verbatim pending a decision on which description to keep; see chat history.]
edgeTolFrac       = 0.5;  % fraction of one grid-step counted as "edge-adjacent" (reporting only)
coverageThreshold = 0.95; % matches the published convention EXACTLY (95% within one grid-step).
                          % Empirically verified on this dataset: all 87 boutons passing the looser
                          % 68% version also pass 95% (74/74 boundary-adjacent, 13/13 interior) --
                          % so using the literal published threshold costs nothing here and lets the
                          % methods section cite the convention directly rather than justify a deviation.
coverageWindowSteps = 1;  % literal "one step size," matching the published wording exactly
sigmaBoundTolFrac = 0.9;  % fraction of the upper sigma bound counted as "pinned" (signature of a
                          % receptive field broader than the sampled range -- never turns over)
nBootstrap        = 50;  % consider starting with e.g. 100 for a faster first pass -- this now runs
                          % for every fit above the R^2 floor, not just edge-adjacent ones

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
screenAzLimits = [-80, 25];   % NaN here means "no clipping on this side" -- flat extension applies
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

% CHANGED: fit ALL ROIs, not just responsive ones. isResponsive is left
% untouched on allRFMapping so it remains available as a separate
% downstream filter -- see the top-of-file NOTE ON SCOPE.
% fitBoutonIdx = 1:numel(allRFMapping);
fitBoutonIdx = respIdxList; % the ACTUAL responsive bouton indices, not 1:numel(respIdxList)
                             % (that generated 1,2,3,...N -- an arbitrary early slice of
                             % allRFMapping, NOT the responsive population)
numBoutons   = numel(fitBoutonIdx);

% Clear any gaussFit_* fields left over from a PREVIOUS run with different criteria
% (e.g. a different r2Floor) -- otherwise boutons that don't clear this run's floor
% keep stale field values from before, silently contaminating downstream scripts.
staleFields = {'gaussFit_A','gaussFit_Az0','gaussFit_El0','gaussFit_sigmaAz','gaussFit_sigmaEl', ...
    'gaussFit_R2', ...
    'gaussFit_isBeyondRange','gaussFit_isBeyondRangeAz','gaussFit_isBeyondRangeEl', ...
    'gaussFit_isAtBoundary','gaussFit_isAtBoundaryAz','gaussFit_isAtBoundaryEl', ...
    'gaussFit_isDegenerateSigma','gaussFit_isDegenerateSigmaAz','gaussFit_isDegenerateSigmaEl', ...
    'gaussFit_isStable','gaussFit_isTrusted','gaussFit_isTrustedAz','gaussFit_isTrustedEl', ...
    'gaussFit_bootCenterStdAz','gaussFit_bootCenterStdEl'};
for iROI = fitBoutonIdx(:)'
    for f = staleFields
        if isfield(allRFMapping, f{1})
            allRFMapping(iROI).(f{1}) = [];
        end
    end
end

gaussFitResults = struct('iROI', {}, 'isResponsive', {}, 'A', {}, 'x0', {}, 'y0', {}, 'sigmaX', {}, 'sigmaY', {}, 'R2', {}, ...
    'isBeyondRangeAz', {}, 'isBeyondRangeEl', {}, 'isBeyondRange', {}, ...
    'isAtBoundaryAz', {}, 'isAtBoundaryEl', {}, 'isAtBoundary', {}, ...
    'isDegenerateSigmaAz', {}, 'isDegenerateSigmaEl', {}, 'isDegenerateSigma', {}, ...
    'bootCenterStdAz', {}, 'bootCenterStdEl', {}, ...
    'bootSigmaStdAz', {}, 'bootSigmaStdEl', {}, ...
    'pctWithinStepAz', {}, 'pctWithinStepEl', {}, ...
    'fracSigmaPinnedAz', {}, 'fracSigmaPinnedEl', {}, ...
    'isTrustedAz', {}, 'isTrustedEl', {}, 'isTrusted', {}, 'nBootValid', {});

fprintf('Fitting 2D Gaussian RFs for %d ROIs (ALL ROIs, regardless of responsiveness; R^2 floor = %.2f)...\n', numBoutons, r2Floor);
fprintf('Bootstrapping ALL fits above the floor (%d resamples each) -- this may take a while.\n', nBootstrap);

allR2 = nan(numBoutons, 1);
sigmaUB = 3*max(range(azRange), range(elRange)); % matches ub(4)/ub(5) below -- kept as one variable for clarity

for k = 1:numBoutons
    iROI = fitBoutonIdx(k);
    b = allRFMapping(iROI);

    % Build the response vector DIRECTLY from baselineSubtracted using respIdx, rather than
    % trusting allRFMapping.meanGridResponse (which may have been computed over a different,
    % e.g. full-trial, window upstream). This guarantees the ORIGINAL point-estimate fit uses
    % the EXACT SAME response window as the bootstrap resampling below (which has always used
    % respIdx), eliminating any window mismatch between the point estimate and its own
    % bootstrap. baselineSubtracted is grid-positions-only (blank trials live separately in
    % baselineSubtractedBlank and are never touched here), so blank is correctly excluded.
    trialMatrix = b.baselineSubtracted; % cell array [nEl x nAz], GRID POSITIONS ONLY
    if isempty(trialMatrix)
        continue;
    end
    gridDims = size(trialMatrix);
    respVec  = nan(numel(trialMatrix), 1);
    for posIdx = 1:numel(trialMatrix)
        [rIdx, cIdx] = ind2sub(gridDims, posIdx);
        trials = trialMatrix{rIdx, cIdx};
        if isempty(trials)
            continue;
        end
        trials = double(trials);
        respVec(posIdx) = mean(mean(trials(:, respIdx), 2, 'omitnan'), 'omitnan');
    end
    if any(isnan(respVec))
        warning('Bouton %d has incomplete grid coverage (NaN in respVec) -- skipping.', iROI);
        continue;
    end
    meanGridResponse = reshape(respVec, gridDims); % kept only for downstream size/ind2sub reference below

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
    %  - isBeyondRange*: fitted center is genuinely OUTSIDE the tested min/max on that axis.
    %    This is the suspicious case (extrapolation beyond what was ever measured) and should
    %    be excluded from the trusted set.
    %  - isAtBoundary*: fitted center sits AT or NEAR the outermost tested position, but not
    %    beyond it. This can be a completely real preference (the bouton may genuinely prefer
    %    the most peripheral position you tested) and should NOT be excluded just for that.
    % isBeyondRange* uses the STIMULUS-FOOTPRINT-extended range: a center is only flagged as
    % genuine extrapolation if it falls outside the visual space the outermost stimulus
    % actually covered, not just outside the bare grid-center coordinates.
    % NEW: computed per-axis so a bad fit on one axis doesn't disqualify a good fit on the other.
    isBeyondRangeAz = (pFit(2) < azRangeExt(1)) || (pFit(2) > azRangeExt(2));
    isBeyondRangeEl = (pFit(3) < elRangeExt(1)) || (pFit(3) > elRangeExt(2));
    isBeyondRange   = isBeyondRangeAz || isBeyondRangeEl; % kept for backward compatibility / reporting

    % isAtBoundary* still refers to the original grid-center range -- this is reporting "did
    % the fit land at/near a position you actually tested," a different question from whether
    % it's extrapolated. Reporting-only, not a gate, per-axis kept for symmetry with the above.
    isAtBoundaryAz = ~isBeyondRangeAz && ...
                     ((pFit(2) <= azRange(1) + edgeTolFrac*azStep) || (pFit(2) >= azRange(2) - edgeTolFrac*azStep));
    isAtBoundaryEl = ~isBeyondRangeEl && ...
                     ((pFit(3) <= elRange(1) + edgeTolFrac*elStep) || (pFit(3) >= elRange(2) - edgeTolFrac*elStep));
    isAtBoundary   = isAtBoundaryAz || isAtBoundaryEl;

    % NOTE: the fitting bound on sigma (ub(4)/ub(5) = 3x the tested range) was set loose
    % deliberately, to let a genuinely broad RF be estimated rather than artificially capped.
    % But that means the earlier "fraction pinned at upper bound" check could never catch a
    % degenerate fit (sigma ~60-100 deg on a ~100 deg grid is clearly not a real RF estimate,
    % but is nowhere near 3x the range). This is a DIFFERENT, more honest check: is sigma bigger
    % than the actual tested span itself? If so, the fit cannot be distinguishing "real broad RF"
    % from "flat noise the optimizer gave up on" -- exclude it regardless of R^2 or stability.
    % NEW: computed per-axis. A bouton can have a perfectly good, narrow azimuth sigma but a
    % degenerate elevation sigma (e.g. because the elevation grid has a narrower tested range) --
    % previously this excluded the WHOLE fit; now only the affected axis is excluded.
    isDegenerateSigmaAz = pFit(4) >= range(azRange);
    isDegenerateSigmaEl = pFit(5) >= range(elRange);
    isDegenerateSigma   = isDegenerateSigmaAz || isDegenerateSigmaEl; % kept for backward compatibility / reporting

    %% bootstrap: runs for every fit above the floor
    bootStdAz = NaN; bootStdEl = NaN; bootSigStdAz = NaN; bootSigStdEl = NaN;
    fracPinAz = NaN; fracPinEl = NaN; pctWithinAz = NaN; pctWithinEl = NaN; nBootValid = 0;
    isCoverageOKAz = false; isCoverageOKEl = false;

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
            % of the ORIGINAL fit's center, per axis? Require >= coverageThreshold (95%) per axis.
            pctWithinAz = mean(abs(bootCenters(validRows, 1) - pFit(2)) <= coverageWindowSteps*azStep);
            pctWithinEl = mean(abs(bootCenters(validRows, 2) - pFit(3)) <= coverageWindowSteps*elStep);

            isCoverageOKAz = pctWithinAz >= coverageThreshold;
            isCoverageOKEl = pctWithinEl >= coverageThreshold;
        end
    end

    % Per-axis trust: coverage-stable on that axis AND not beyond-range on that axis AND
    % sigma not degenerate on that axis. A fit can't be trustworthy on an axis whose sigma is
    % degenerate, regardless of bootstrap center stability -- a wide flat "hill" can still
    % refit to a similar (degenerate) shape consistently across resamples without that meaning
    % the fit is meaningful on that axis.
    isTrustedAz = isCoverageOKAz && ~isBeyondRangeAz && ~isDegenerateSigmaAz;
    isTrustedEl = isCoverageOKEl && ~isBeyondRangeEl && ~isDegenerateSigmaEl;
    % Combined flag, TRUE only if BOTH axes are trusted -- equivalent to the old all-or-nothing
    % behavior. Use this for analyses that need a fully-trusted 2D fit (e.g. sigmaAz vs sigmaEl
    % comparisons); use the per-axis flags above for analyses that only need one axis reliable.
    isTrusted = isTrustedAz && isTrustedEl;

    gaussFitResults(end+1) = struct('iROI', iROI, 'isResponsive', b.isResponsive, ...
        'A', pFit(1), 'x0', pFit(2), 'y0', pFit(3), 'sigmaX', pFit(4), 'sigmaY', pFit(5), 'R2', R2, ...
        'isBeyondRangeAz', isBeyondRangeAz, 'isBeyondRangeEl', isBeyondRangeEl, 'isBeyondRange', isBeyondRange, ...
        'isAtBoundaryAz', isAtBoundaryAz, 'isAtBoundaryEl', isAtBoundaryEl, 'isAtBoundary', isAtBoundary, ...
        'isDegenerateSigmaAz', isDegenerateSigmaAz, 'isDegenerateSigmaEl', isDegenerateSigmaEl, 'isDegenerateSigma', isDegenerateSigma, ...
        'bootCenterStdAz', bootStdAz, 'bootCenterStdEl', bootStdEl, ...
        'bootSigmaStdAz', bootSigStdAz, 'bootSigmaStdEl', bootSigStdEl, ...
        'pctWithinStepAz', pctWithinAz, 'pctWithinStepEl', pctWithinEl, ...
        'fracSigmaPinnedAz', fracPinAz, 'fracSigmaPinnedEl', fracPinEl, ...
        'isTrustedAz', isTrustedAz, 'isTrustedEl', isTrustedEl, 'isTrusted', isTrusted, ...
        'nBootValid', nBootValid); %#ok<SAGROW>

    allRFMapping(iROI).gaussFit_A       = pFit(1);
    allRFMapping(iROI).gaussFit_Az0     = pFit(2);
    allRFMapping(iROI).gaussFit_El0     = pFit(3);
    allRFMapping(iROI).gaussFit_sigmaAz = pFit(4);
    allRFMapping(iROI).gaussFit_sigmaEl = pFit(5);
    allRFMapping(iROI).gaussFit_R2      = R2;
    allRFMapping(iROI).gaussFit_isBeyondRangeAz = isBeyondRangeAz;
    allRFMapping(iROI).gaussFit_isBeyondRangeEl = isBeyondRangeEl;
    allRFMapping(iROI).gaussFit_isBeyondRange   = isBeyondRange;
    allRFMapping(iROI).gaussFit_isAtBoundaryAz  = isAtBoundaryAz;
    allRFMapping(iROI).gaussFit_isAtBoundaryEl  = isAtBoundaryEl;
    allRFMapping(iROI).gaussFit_isAtBoundary    = isAtBoundary;
    allRFMapping(iROI).gaussFit_isDegenerateSigmaAz = isDegenerateSigmaAz;
    allRFMapping(iROI).gaussFit_isDegenerateSigmaEl = isDegenerateSigmaEl;
    allRFMapping(iROI).gaussFit_isDegenerateSigma   = isDegenerateSigma;
    % Per-axis trusted flags (NEW) -- prefer these for axis-specific analyses.
    allRFMapping(iROI).gaussFit_isTrustedAz = isTrustedAz;
    allRFMapping(iROI).gaussFit_isTrustedEl = isTrustedEl;
    % Combined trusted flag: TRUE only if both axes pass. Equivalent to previous all-or-nothing
    % behavior -- R^2 floor + bootstrap coverage (95% within one grid-step, both axes) + non-
    % degenerate sigma (both axes) + not-beyond-range (both axes). Remember to additionally
    % filter by allRFMapping.isResponsive downstream, since this script now fits ALL ROIs and a
    % non-responsive bouton's noise can still occasionally clear these gates by chance.
    allRFMapping(iROI).gaussFit_isTrusted = isTrusted;
    allRFMapping(iROI).gaussFit_bootCenterStdAz = bootStdAz;
    allRFMapping(iROI).gaussFit_bootCenterStdEl = bootStdEl;
end

nAboveFloor = numel(gaussFitResults);
fprintf('\n%d / %d ROIs passed the R^2 >= %.2f floor.\n', nAboveFloor, numBoutons, r2Floor);

nResponsiveAboveFloor = sum([gaussFitResults.isResponsive]);
fprintf('Of these, %d / %d are flagged isResponsive (the rest passed the R^2 floor on non-responsive data --\n', ...
    nResponsiveAboveFloor, nAboveFloor);
fprintf('remember to AND with isResponsive downstream, not rely on gaussFit_isTrusted alone).\n');

isBeyondRangeAzAll = [gaussFitResults.isBeyondRangeAz];
isBeyondRangeElAll = [gaussFitResults.isBeyondRangeEl];
isAtBoundaryAzAll  = [gaussFitResults.isAtBoundaryAz];
isAtBoundaryElAll  = [gaussFitResults.isAtBoundaryEl];
isDegenAzAll       = [gaussFitResults.isDegenerateSigmaAz];
isDegenElAll       = [gaussFitResults.isDegenerateSigmaEl];
isTrustedAzAll     = [gaussFitResults.isTrustedAz];
isTrustedElAll     = [gaussFitResults.isTrustedEl];
isTrustedAll       = [gaussFitResults.isTrusted];

fprintf('\n--- Azimuth axis ---\n');
fprintf('%d / %d centers genuinely OUTSIDE the tested range (excluded on this axis).\n', sum(isBeyondRangeAzAll), nAboveFloor);
fprintf('%d / %d centers AT/NEAR the tested boundary but within range (KEPT).\n', sum(isAtBoundaryAzAll), nAboveFloor);
fprintf('%d / %d have a degenerate (span-exceeding) sigma (excluded on this axis).\n', sum(isDegenAzAll), nAboveFloor);
fprintf('%d / %d pass full azimuth trust (isTrustedAz).\n', sum(isTrustedAzAll), nAboveFloor);

fprintf('\n--- Elevation axis ---\n');
fprintf('%d / %d centers genuinely OUTSIDE the tested range (excluded on this axis).\n', sum(isBeyondRangeElAll), nAboveFloor);
fprintf('%d / %d centers AT/NEAR the tested boundary but within range (KEPT).\n', sum(isAtBoundaryElAll), nAboveFloor);
fprintf('%d / %d have a degenerate (span-exceeding) sigma (excluded on this axis).\n', sum(isDegenElAll), nAboveFloor);
fprintf('%d / %d pass full elevation trust (isTrustedEl).\n', sum(isTrustedElAll), nAboveFloor);

fprintf('\n--- Combined (both axes trusted) ---\n');
fprintf('%d / %d pass the trusted set (isTrustedAz AND isTrustedEl).\n', sum(isTrustedAll), nAboveFloor);
fprintf('%d / %d of those are also isResponsive -- this is the set to use for "real, well-fit RF" analyses.\n', ...
    sum(isTrustedAll & [gaussFitResults.isResponsive]), sum(isTrustedAll));

fracPinAzAll = [gaussFitResults.fracSigmaPinnedAz];
fracPinElAll = [gaussFitResults.fracSigmaPinnedEl];
fprintf('\nMedian fraction of bootstrap resamples with sigma pinned at its upper bound:\n');
fprintf('  Azimuth:   %.1f%%\n', 100*median(fracPinAzAll, 'omitnan'));
fprintf('  Elevation: %.1f%%\n', 100*median(fracPinElAll, 'omitnan'));
fprintf(['NOTE: a high pinned fraction means the bootstrap keeps failing to see the response turn back down --\n' ...
    'that axis'' true RF width likely exceeds your sampled range for a meaningful fraction of boutons,\n' ...
    'regardless of R^2. Treat sigma (and possibly center) on that axis with real caution for those boutons.\n']);

%% ===================== combined R^2-floor x (per-axis trust) sensitivity =====================
fprintf('\n=== Combined sensitivity: R^2 floor x [isTrustedAz AND isTrustedEl] (n = %d fitted ROIs) ===\n', numBoutons);
threshList = [0.1 0.2 0.3 0.4 0.5 0.6 0.7];
for th = threshList
    passR2 = allR2 >= th;
    trustedForThresh = false(numBoutons, 1);
    for gi = 1:numel(gaussFitResults)
        kMatch = find(fitBoutonIdx == gaussFitResults(gi).iROI, 1);
        if ~isempty(kMatch) && gaussFitResults(gi).R2 >= th
            trustedForThresh(kMatch) = gaussFitResults(gi).isTrustedAz && gaussFitResults(gi).isTrustedEl;
        end
    end
    nR2Only = sum(passR2);
    nTrusted = sum(passR2 & trustedForThresh);
    fprintf('  R^2 >= %.1f : %d pass R^2 alone | %d pass ALL checks (%.1f%% of R^2-passing)\n', ...
        th, nR2Only, nTrusted, 100*nTrusted/max(nR2Only,1));
end

nPassed = numel(gaussFitResults);
fprintf('\n%d / %d ROIs passed the R^2 floor (%.2f).\n', nPassed, numBoutons, r2Floor);
fprintf('Use the combined-sensitivity table above to choose your actual final R^2 x stability criteria.\n');

%% ===================== coverage summary (uses stored per-bouton coverage, no re-run needed) =====================
pctWithinAzAll = [gaussFitResults.pctWithinStepAz];
pctWithinElAll = [gaussFitResults.pctWithinStepEl];
fprintf('\n=== Coverage summary (n = %d fitted ROIs, R^2 floor only) ===\n', nPassed);
fprintf('Median %% of resamples within one grid-step: azimuth = %.1f%%, elevation = %.1f%%\n', ...
    100*median(pctWithinAzAll, 'omitnan'), 100*median(pctWithinElAll, 'omitnan'));
fprintf('Coverage threshold in use: %.0f%% (matches published convention); %d / %d pass this + non-degenerate + not-beyond-range on BOTH axes.\n', ...
    100*coverageThreshold, sum(isTrustedAll), nPassed);

%%
r2Thresh = 0.7;   % set whatever threshold you want to test
R2vals = [gaussFitResults.R2];
idx = R2vals >= r2Thresh;

fprintf('n = %d fits with R^2 >= %.1f\n', sum(idx), r2Thresh);
fprintf('Mean sigma azimuth:   %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaX]), std([gaussFitResults(idx).sigmaX]));
fprintf('Mean sigma elevation: %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaY]), std([gaussFitResults(idx).sigmaY]));

%%
r2Thresh = 0.2;   
R2vals = [gaussFitResults.R2];
idx = R2vals >= r2Thresh & [gaussFitResults.isTrusted];

fprintf('n = %d fits with R^2 >= %.1f\n', sum(idx), r2Thresh);
fprintf('Mean sigma azimuth:   %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaX]), std([gaussFitResults(idx).sigmaX]));
fprintf('Mean sigma elevation: %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaY]), std([gaussFitResults(idx).sigmaY]));

%%
r2Thresh = 0.7;   % set whatever threshold you want to test
R2vals = [gaussFitResults.R2];
idx = R2vals >= r2Thresh;

fprintf('n = %d fits with R^2 >= %.1f\n', sum(idx), r2Thresh);
fprintf('Mean sigma azimuth:   %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaX]), std([gaussFitResults(idx).sigmaX]));
fprintf('Mean sigma elevation: %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaY]), std([gaussFitResults(idx).sigmaY]));

%%
r2Thresh = 0.3;   
R2vals = [gaussFitResults.R2];
idx = R2vals >= r2Thresh & [gaussFitResults.isTrusted];

fprintf('n = %d fits with R^2 >= %.1f\n', sum(idx), r2Thresh);
fprintf('Mean sigma azimuth:   %.2f +/- %.2f deg\n', median([gaussFitResults(idx).sigmaX]), std([gaussFitResults(idx).sigmaX]));
fprintf('Mean sigma elevation: %.2f +/- %.2f deg\n', median([gaussFitResults(idx).sigmaY]), std([gaussFitResults(idx).sigmaY]));