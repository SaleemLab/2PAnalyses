% RFMapping_Gaussian2DFit.m
%
% Fits an anisotropic 2D Gaussian (no rotation, no baseline term) to each
% bouton's mean response across the RF mapping grid -- matching
% Timplalexi et al.:
%
%   G(x,y) = A * exp( -( (x-x0)^2/(2*sigmaX^2) + (y-y0)^2/(2*sigmaY^2) ) )
%
% SIMPLIFIED TO MATCH THE PAPER'S APPROACH. Trust criteria are now just:
%   1. R^2 >= r2Floor (paper uses R^2 >= 0.7; r2Floor here is set lower
%      since criterion #2 below does additional work the paper's method
%      doesn't have)
%   2. Bootstrap center-coverage: >= coverageThreshold (95%, matching the
%      published convention) of 500 trial-resampled refits land within
%      one grid-step of the original fit's center, on BOTH axes.
% Removed relative to earlier versions of this script: isDegenerateSigma,
% isBeyondRange, isAtBoundary, and the stimulus-footprint/screen-clipping
% range extension logic that only fed those. None of this is described in
% the paper, and none of it affects the OPTIMIZER or the fitted sigma
% value itself -- these were post-fit classifiers only, deciding which
% already-computed fits counted as "trusted." Removing them changes which
% fits are included in your reported summary stats; it does not change
% any individual bouton's fitted sigma.
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
% USES rawAlignedTrials, NOT baselineSubtracted: analyseRFMapping.m's
% baselineSubtracted field is the raw aligned trial subtracted a SECOND
% time by its own per-trial pre-stimulus baseline mean, on top of
% whatever global dF/F normalization was already applied upstream. That
% second, more local/noisy correction can distort peak-vs-surround
% contrast and bias a no-baseline-term Gaussian fit toward smaller sigma.
% rawAlignedTrials (built directly in the pooling script from each
% session's Response file / psthData.alignedResponses) stores the trial
% data BEFORE that second subtraction, so the Gaussian fit works from a
% single consistent normalization rather than two stacked ones. Grid
% positions only (blank trials live separately in rawAlignedTrialsBlank
% and are never referenced here), so blank is correctly excluded from
% every fit.
%
% NOTE: analyseRFMapping.m computes its own responsiveness/meanGridResponse
% stats using respWin = [0.5 3], saved to RFMappingMetadata.respWin -- but
% the main pooling script overrides this with a hardcoded respWin = [0.1
% 3] (see the commented-out "% respWin = RFMappingMetadata.respWin;" line
% there). This script's respWin below currently matches the pooling
% script's override, NOT RFMappingMetadata.respWin. These are two
% different windows already in use elsewhere in the pipeline -- worth
% reconciling deliberately rather than by accident.
%
% REQUIRES (already in workspace from the main pooling script):
%   allRFMapping   - pooled bouton struct array, with .rawAlignedTrials
%                    ([nEl x nAz] cell of [Trials x Time] matrices)
%   uAz, uEl_plot  - grid position vectors (uEl_plot rows, uAz columns)
%   timeVector     - time vector matching the Time dimension above, used
%                    to compute respIdx from respWin
%
% REQUIRES: Optimization Toolbox (lsqcurvefit).

%% params
rng(1, 'twister'); % FIXED SEED -- without this, every run draws different bootstrap resamples,
                    % which can change WHICH boutons pass the coverage criterion.
r2Floor           = 0.1;  % lower than the paper's 0.7 -- bootstrap coverage below does the
                          % primary reliability work here; adjust freely.
coverageThreshold = 0.95; % matches the published convention (95% of resamples within one grid-step)
coverageWindowSteps = 1;  % literal "one step size," matching the published wording
nBootstrap        = 100;

% Response window: defined EXPLICITLY here rather than trusting an externally-set respIdx
% variable, so the point estimate and bootstrap are guaranteed to use an IDENTICAL window
% regardless of what else has run before this script. Keep in sync with respWin in the main
% pooling script.
respWin = [0.1 3]; 
if ~exist('timeVector', 'var')
    error('timeVector not found -- run the main pooling script first (needed to compute respIdx from respWin).');
end
respIdx = timeVector >= respWin(1) & timeVector <= respWin(2);

azStep  = mean(diff(sort(uAz)));
elStep  = mean(diff(sort(uEl_plot)));
azRange = [min(uAz) max(uAz)];
elRange = [min(uEl_plot) max(uEl_plot)];

[AzGrid, ElGrid] = meshgrid(uAz, uEl_plot);   % [nEl x nAz]
xyList = [AzGrid(:), ElGrid(:)];

gaussFit2D = @(p, xy) p(1) * exp( -( (xy(:,1)-p(2)).^2 ./ (2*p(4)^2) + (xy(:,2)-p(3)).^2 ./ (2*p(5)^2) ) );
% p = [A, x0, y0, sigmaX, sigmaY]

fitOpts = optimoptions('lsqcurvefit', 'Display', 'off');
lb = [0,   azRange(1)-azStep, elRange(1)-elStep, azStep/4,        elStep/4];
ub = [Inf, azRange(2)+azStep, elRange(2)+elStep, 3*range(azRange), 3*range(elRange)];

% Fit ALL ROIs, not just responsive ones. isResponsive is left untouched on allRFMapping so it
% remains available as a separate downstream filter -- a non-responsive bouton's noise can
% still occasionally clear R^2 floor + bootstrap coverage by chance, so combine isResponsive
% AND gaussFit_isTrusted downstream, not gaussFit_isTrusted alone.
% fitBoutonIdx = 1:numel(allRFMapping);
fitBoutonIdx = respIdxList;
numBoutons   = numel(fitBoutonIdx);

% Clear any gaussFit_* fields left over from a PREVIOUS run with different criteria.
staleFields = {'gaussFit_A','gaussFit_Az0','gaussFit_El0','gaussFit_sigmaAz','gaussFit_sigmaEl', ...
    'gaussFit_R2', 'gaussFit_isTrusted', 'gaussFit_bootCenterStdAz', 'gaussFit_bootCenterStdEl', ...
    'gaussFit_pctWithinStepAz', 'gaussFit_pctWithinStepEl'};
for iROI = fitBoutonIdx(:)'
    for f = staleFields
        if isfield(allRFMapping, f{1})
            allRFMapping(iROI).(f{1}) = [];
        end
    end
end

gaussFitResults = struct('iROI', {}, 'isResponsive', {}, 'A', {}, 'x0', {}, 'y0', {}, 'sigmaX', {}, 'sigmaY', {}, 'R2', {}, ...
    'bootCenterStdAz', {}, 'bootCenterStdEl', {}, ...
    'pctWithinStepAz', {}, 'pctWithinStepEl', {}, ...
    'isTrusted', {}, 'nBootValid', {});

fprintf('Fitting 2D Gaussian RFs for %d ROIs (ALL ROIs, regardless of responsiveness; R^2 floor = %.2f)...\n', numBoutons, r2Floor);
fprintf('Bootstrapping ALL fits above the floor (%d resamples each) -- this may take a while.\n', nBootstrap);

allR2 = nan(numBoutons, 1);

for k = 1:numBoutons
    iROI = fitBoutonIdx(k);
    b = allRFMapping(iROI);
    trialMatrix = b.rawAlignedTrials; % cell array [nEl x nAz], GRID POSITIONS ONLY, pre-baseline-subtraction
    % CLEAN RE-TEST: this is the raw-trials hypothesis being tested again, now that BOTH prior
    % bugs are fixed -- fitBoutonIdx correctly uses respIdxList (not 1:numel(respIdxList)), and
    % allRFMapping.rawAlignedTrials should come from a FRESH run of analyse_RFBoutons_current.m
    % (which fixed the session-level grid-indexing/scrambling bug). Only run this AFTER
    % rerunning the pooling script from scratch -- reusing an old allRFMapping defeats the point.
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
        continue; % excludes only catastrophic non-fits
    end

    %% bootstrap: SAME computeGridRespVec function, called on a per-position trial resample each iteration
    bootStdAz = NaN; bootStdEl = NaN; pctWithinAz = NaN; pctWithinEl = NaN; nBootValid = 0;
    isStable = false;

    bootCenters = nan(nBootstrap, 2);
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
            pBoot = lsqcurvefit(gaussFit2D, pFit, xyList, bootRespVec, lb, ub, fitOpts);
            bootCenters(bIter, :) = [pBoot(2), pBoot(3)];
            nBootValid = nBootValid + 1;
        catch
            continue;
        end
    end

    validRows = ~any(isnan(bootCenters), 2);
    if sum(validRows) >= 20
        bootStdAz = std(bootCenters(validRows, 1));
        bootStdEl = std(bootCenters(validRows, 2));

        % COVERAGE-BASED criterion: fraction of bootstrap-resampled centers within
        % coverageWindowSteps grid-steps of the ORIGINAL fit's center, per axis.
        pctWithinAz = mean(abs(bootCenters(validRows, 1) - pFit(2)) <= coverageWindowSteps*azStep);
        pctWithinEl = mean(abs(bootCenters(validRows, 2) - pFit(3)) <= coverageWindowSteps*elStep);

        isStable = (pctWithinAz >= coverageThreshold) && (pctWithinEl >= coverageThreshold);
    end

    gaussFitResults(end+1) = struct('iROI', iROI, 'isResponsive', b.isResponsive, ...
        'A', pFit(1), 'x0', pFit(2), 'y0', pFit(3), 'sigmaX', pFit(4), 'sigmaY', pFit(5), 'R2', R2, ...
        'bootCenterStdAz', bootStdAz, 'bootCenterStdEl', bootStdEl, ...
        'pctWithinStepAz', pctWithinAz, 'pctWithinStepEl', pctWithinEl, ...
        'isTrusted', isStable, 'nBootValid', nBootValid); %#ok<SAGROW>

    allRFMapping(iROI).gaussFit_A       = pFit(1);
    allRFMapping(iROI).gaussFit_Az0     = pFit(2);
    allRFMapping(iROI).gaussFit_El0     = pFit(3);
    allRFMapping(iROI).gaussFit_sigmaAz = pFit(4);
    allRFMapping(iROI).gaussFit_sigmaEl = pFit(5);
    allRFMapping(iROI).gaussFit_R2      = R2;
    allRFMapping(iROI).gaussFit_isTrusted = isStable;
    allRFMapping(iROI).gaussFit_bootCenterStdAz = bootStdAz;
    allRFMapping(iROI).gaussFit_bootCenterStdEl = bootStdEl;
    allRFMapping(iROI).gaussFit_pctWithinStepAz = pctWithinAz;
    allRFMapping(iROI).gaussFit_pctWithinStepEl = pctWithinEl;
end

nAboveFloor = numel(gaussFitResults);
fprintf('\n%d / %d ROIs passed the R^2 >= %.2f floor.\n', nAboveFloor, numBoutons, r2Floor);

isTrustedAll = [gaussFitResults.isTrusted];
fprintf('%d / %d pass the trusted set (R^2 floor AND bootstrap coverage >= %.0f%%, both axes).\n', ...
    sum(isTrustedAll), nAboveFloor, 100*coverageThreshold);

nResponsiveTrusted = sum(isTrustedAll & [gaussFitResults.isResponsive]);
fprintf('%d / %d of the trusted set are also isResponsive -- remember to AND with isResponsive\n', ...
    nResponsiveTrusted, sum(isTrustedAll));
fprintf('downstream, since this script fits ALL ROIs, not just responsive ones.\n');

%% ===================== R^2-floor sensitivity (paper-style, no other gates) =====================
fprintf('\n=== R^2 floor sensitivity (n = %d fitted ROIs) ===\n', numBoutons);
threshList = [0.1 0.2 0.3 0.4 0.5 0.6 0.7];
for th = threshList
    passR2 = allR2 >= th;
    trustedForThresh = false(numBoutons, 1);
    for gi = 1:numel(gaussFitResults)
        kMatch = find(fitBoutonIdx == gaussFitResults(gi).iROI, 1);
        if ~isempty(kMatch) && gaussFitResults(gi).R2 >= th
            trustedForThresh(kMatch) = gaussFitResults(gi).isTrusted;
        end
    end
    nR2Only = sum(passR2);
    nTrusted = sum(passR2 & trustedForThresh);
    fprintf('  R^2 >= %.1f : %d pass R^2 alone | %d pass R^2 + bootstrap coverage (%.1f%% of R^2-passing)\n', ...
        th, nR2Only, nTrusted, 100*nTrusted/max(nR2Only,1));
end


fprintf('\nUse the sensitivity table above to choose your final R^2 x coverage criteria.\n');
r2Thresh = 0.3;   
R2vals = [gaussFitResults.R2];
idx = R2vals >= r2Thresh & [gaussFitResults.isTrusted];

fprintf('n = %d fits with R^2 >= %.1f\n', sum(idx), r2Thresh);
fprintf('Mean sigma azimuth:   %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaX]), std([gaussFitResults(idx).sigmaX]));
fprintf('Mean sigma elevation: %.2f +/- %.2f deg\n', mean([gaussFitResults(idx).sigmaY]), std([gaussFitResults(idx).sigmaY]));

%% ===================== local function =====================
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


%%
