% Loads DirTuning response.mat per session, pools boutons into a struct
% array for the ANOVA/responsiveness script.
%
% UPDATED to match the current getTrialResponsePSTH.m: alignedResponses
% is NOT baseline-subtracted at all upstream anymore. This script performs the single
% baseline subtraction here, using baseWin = [-1, 0], on the trace
% trace.
%
% Data shape: response.psthData is [8x1], one entry per direction.
%   response.psthData(d).alignedResponses : [nBoutons x nTimepoints x nTrials] (interpolated, NOT baseline-subtracted)
%   response.psthData(d).timeVector       : [nTimepoints x 1]
%
% Timing variants:
%   postStimTime == 4 -> 2s-on/2s-off -> respWin = [0.5 3]
%   postStimTime == 3 -> 1s-on/2s-off -> respWin = [0.5 2]
%   baseWin = [-1 0] for both.

% soma
% DirTuningTable = filterMasterTable_usingNameSessionPairs('MouseID', ...
%     {'M25131', 'M25126'}, 'Exclude', 0, 'HasStimulus', {'DirTuning'});

DirTuningTable = filterMasterTable_usingNameSessionPairs('MouseID', ...
    {'M25132','M25133', 'M26003'}, 'Exclude', 0, 'HasStimulus', {'DirTuning'});
allMice    = DirTuningTable.MouseID;
uniqueMice = unique(allMice, 'stable');

% window definitions per timing variant
baseWin      = [-0.75 0];    % same for both variants
respWin_2sOn = [0.2 3];   % postStimTime == 4 (2s-on/2s-off)
respWin_1sOn = [0.2 3];   % postStimTime == 3 (1s-on/2s-off)

allDirTuning     = struct('baselineSubtracted', {}, 'fullTraceSub', {}, 'timeVec', {}, ...
                           'stimOnDuration', {}, 'meanDirResponse', {}, ...
                           'trialMeanResp', {}, 'trialRawResp', {}, 'trialBaselineVals', {}, ...
                           'stimValues', {}, 'stimVariant', {}, 'sessionLabel', {});
sessionLabels    = {};
haveWarnedShape  = false;

%% loop over mice and pool data
for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));
    fprintf('MOUSE: %s | %d sessions\n', thisMouse, length(mouseSessIdx));
    %% loop over sessions for this mouse
    for iSess = 1:length(mouseSessIdx)
        tableRow        = DirTuningTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        thisDay         = tableRow.DayOfExperience;
        thisArea        = tableRow.TargetArea;
        fprintf('\n--- Session: %s | Day %d | Area: %s ---\n', ...
            thisSessionName, thisDay, char(thisArea));

        %% load sessionFileInfo
        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath)
            warning('sfi missing for %s — skipping.', thisSessionName);
            continue;
        end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        DirTuneIdx = find(contains(stimNames, 'DirTuning'));
        if isempty(DirTuneIdx)
            warning('No DirTuning files for %s — skipping.', thisSessionName);
            continue;
        end

        %% load the response file for this session
        try
            load(sessionFileInfo.stimFiles(DirTuneIdx(1)).Response, 'response');
        catch ME
            warning('  Could not load response for %s: %s', thisSessionName, ME.message);
            continue;
        end
        if ~isfield(response, 'psthData') || isempty(response.psthData)
            warning('  psthData empty/missing for %s — skipping.', thisSessionName);
            continue;
        end

        %% pick respWin based on this session's timing variant
        switch response.postStimTime
            case 4
                respWin = respWin_2sOn; stimVariant = 4;
            case 3
                respWin = respWin_1sOn; stimVariant = 3;
            otherwise
                warning(['  Unexpected postStimTime (%.2f) for %s — skipping.'], ...
                    response.postStimTime, thisSessionName);
                continue;
        end

        nDir      = numel(response.psthData);
        timeVec   = response.psthData(1).timeVector(:)';
        baseIdx   = timeVec >= baseWin(1) & timeVec <= baseWin(2);
        respIdx   = timeVec >= respWin(1) & timeVec <= respWin(2);
        stimVals  = arrayfun(@(s) s.stimValue, response.psthData);

        nBoutons = size(response.psthData(1).alignedResponses, 1);
        if ~haveWarnedShape
            fprintf('  Detected shape: %d boutons x %d timepoints x %d trials (dir 1)\n', ...
                size(response.psthData(1).alignedResponses,1), ...
                size(response.psthData(1).alignedResponses,2), ...
                size(response.psthData(1).alignedResponses,3));
            haveWarnedShape = true;
        end

        %% per-bouton baseline subtraction + saving
        for iBouton = 1:nBoutons
            baselineSubtracted = cell(nDir, 1);
            fullTraceSub        = cell(nDir, 1); % NEW: full pre+post trace, baseline-subtracted (not just respIdx slice)
            meanDirResponse     = nan(nDir, 1);
            trialMeanResp       = cell(nDir, 1); % per-trial resp-window mean, BASELINE-SUBTRACTED (for tuning shape)
            trialRawResp        = cell(nDir, 1); % per-trial resp-window mean, RAW (for responsiveness threshold)
            trialBaselineVals   = cell(nDir, 1); % per-trial baseline-window mean, RAW (the blank-equivalent)

            for thisDir = 1:nDir
                traceMat = squeeze(response.psthData(thisDir).alignedResponses(iBouton, :, :));
                if isvector(traceMat)
                    traceMat = traceMat(:);
                end
                traceMat = double(traceMat)'; % [nTrials x nTimepoints]

                % take the mean during pre onset activity
                perTrialBaseline = mean(traceMat(:, baseIdx), 2, 'omitnan');
                % subtract
                traceMatSub = traceMat - perTrialBaseline;

                baselineSubtracted{thisDir} = traceMatSub(:, respIdx);
                fullTraceSub{thisDir}       = traceMatSub; % NEW: entire timecourse, not sliced to respIdx

                meanDirResponse(thisDir) = mean(mean(traceMatSub(:, respIdx), 2, 'omitnan'), 'omitnan');
                % get mean across trials this direction
                trialMeanResp{thisDir}     = mean(traceMatSub(:, respIdx), 2, 'omitnan'); % [nTrials x 1], baseline-subtracted
                trialRawResp{thisDir}      = mean(traceMat(:, respIdx), 2, 'omitnan');    % [nTrials x 1], RAW, same scale as baseline
                trialBaselineVals{thisDir} = perTrialBaseline;                            % [nTrials x 1], RAW baseline
            end

            % stim-on duration for shading in plots: per your timing
            % convention, postStimTime==4 -> 2s-on/2s-off, postStimTime==3
            % -> 1s-on/2s-off (off period always 2s)
            switch stimVariant
                case 4, stimOnDuration = 2;
                case 3, stimOnDuration = 1;
                otherwise, stimOnDuration = NaN;
            end

            allDirTuning(end+1) = struct( ...
                'baselineSubtracted', {baselineSubtracted}, ...
                'fullTraceSub',       {fullTraceSub}, ...
                'timeVec',            timeVec, ...
                'stimOnDuration',     stimOnDuration, ...
                'meanDirResponse',    meanDirResponse, ...
                'trialMeanResp',      {trialMeanResp}, ...
                'trialRawResp',       {trialRawResp}, ...
                'trialBaselineVals',  {trialBaselineVals}, ...
                'stimValues',         stimVals(:), ...
                'stimVariant',        stimVariant, ...
                'sessionLabel',       sprintf('%s_%s', thisMouse, thisSessionName));
        end

        sessionLabels = [sessionLabels; repmat({sprintf('%s_%s', thisMouse, thisSessionName)}, ...
            nBoutons, 1)];
        fprintf('  Added %d boutons (running total: %d).\n', nBoutons, numel(allDirTuning));
    end
end

fprintf('\nDone. %d boutons pooled across %d sessions.\n', ...
    numel(allDirTuning), numel(unique(sessionLabels)));


%% For each bouton test for responsiveness using anova

% One-way ANOVA of per-trial response-window means across the 8
% directions -> anovaP.

% "Blank"-equivalent responsiveness check: the pre-stimulus baseline
% period (baseWin = [-0.5 0]) is measured on every single trial, so it
% substitutes for a dedicated blank stimulus. Pooled across ALL
% trials and ALL directions for this bouton -> baseMean, baseSD.
% Bouton must have mean response to its preferred direction >
% baseMean + 2*baseSD.
%
% isResponsive(b) = true  <->  anovaP(b) < 0.05  AND  meanPreferred(b) > baseMean(b) + 2*baseSD(b)
%
% Requires allDirTuning to already be in the workspace


alpha = 0.05;
numSD   = 2;

nBoutonsTotal    = numel(allDirTuning);
anovaP           = nan(nBoutonsTotal, 1);
baseMean         = nan(nBoutonsTotal, 1);
baseSD           = nan(nBoutonsTotal, 1);
meanPreferred    = nan(nBoutonsTotal, 1);    % baseline-subtracted, for reference/plotting
meanPreferredRaw = nan(nBoutonsTotal, 1);    % raw, used in the threshold check
prefDirIdx       = nan(nBoutonsTotal, 1);
isResponsive     = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s    = allDirTuning(b);
    nDir = numel(s.trialMeanResp);

    % one-way ANOVA on baseline-subtracted per-trial values (tuning shape, relative across directions) ---
    y   = [];
    grp = [];
    for thisDir = 1:nDir
        vals = s.trialMeanResp{thisDir}(:);
        y    = [y;   vals];
        grp  = [grp; repmat(thisDir, numel(vals), 1)];
    end

    if numel(y) < nDir * 2 || all(isnan(y))
        % not enough trials / all-NaN bouton -> skip, leave as not responsive
        continue;
    end

    validRows = ~isnan(y);
    anovaP(b) = anova1(y(validRows), grp(validRows), 'off');

    % pooled "blank"-equivalent baseline across all trials, all directions (RAW scale)
    baselinePool  = vertcat(s.trialBaselineVals{:});
    baselinePool  = baselinePool(~isnan(baselinePool));
    baseMean(b)   = mean(baselinePool);
    baseSD(b)     = std(baselinePool);

    % preferred direction identified from baseline-subtracted means (relative comparison is fine here)
    [meanPreferred(b), prefDirIdx(b)] = max(s.meanDirResponse);

    % mean RAW response for that same preferred direction (matches the RAW baseline scale)
    meanPreferredRaw(b) = mean(s.trialRawResp{prefDirIdx(b)}, 'omitnan');

    % combined responsiveness criterion (RAW vs RAW, same units)
    isResponsive(b) = (anovaP(b) < alpha) && ...
                     (meanPreferredRaw(b) > baseMean(b) + numSD * baseSD(b));
end

fprintf('\n%d / %d boutons classified as visually responsive (p<%.2f, >%d SD above baseline).\n', ...
    sum(isResponsive), nBoutonsTotal, alpha, numSD);

% attach results back onto allDirTuning; Dont think i need to save this
% anywhere it is fairly quick
for b = 1:nBoutonsTotal
    allDirTuning(b).anovaP           = anovaP(b);
    allDirTuning(b).baseMean         = baseMean(b);
    allDirTuning(b).baseSD           = baseSD(b);
    allDirTuning(b).meanPreferred    = meanPreferred(b);
    allDirTuning(b).meanPreferredRaw = meanPreferredRaw(b);
    allDirTuning(b).prefDirIdx       = prefDirIdx(b);
    allDirTuning(b).isResponsive     = isResponsive(b);
end

%% debugging..
% boutonIndicesToCheck = [165 167 171 173];
%
% if isempty(boutonIndicesToCheck)
%     error('Set boutonIndicesToCheck to the bouton indices you want to inspect.');
% end
%
% for bi = 1:numel(boutonIndicesToCheck)
%     b = boutonIndicesToCheck(bi);
%     s = allDirTuning(b);
%     nDir = numel(s.trialMeanResp);
%
%     fprintf('\n========== Bouton index %d ==========\n', b);
%
%     %
%     fprintf('%-8s %8s %10s %10s\n', 'Dir#', 'nTrials', 'MeanResp', 'MeanBase');
%     for thisDir = 1:nDir
%         nTrials = numel(s.trialMeanResp{thisDir});
%         mResp   = mean(s.trialMeanResp{thisDir}, 'omitnan');
%         mBase   = mean(s.trialBaselineVals{thisDir}, 'omitnan');
%         fprintf('%-8d %8d %10.5f %10.5f\n', thisDir, nTrials, mResp, mBase);
%     end
%
%     %
%     y = []; grp = [];
%     for thisDir = 1:nDir
%         vals = s.trialMeanResp{thisDir}(:);
%         y = [y; vals];
%         grp = [grp; repmat(thisDir, numel(vals), 1)];
%     end
%     validRows = ~isnan(y);
%     anovaP = anova1(y(validRows), grp(validRows), 'off');
%     fprintf('\nANOVA p-value: %.4g  (threshold: p < %.2f -> %s)\n', ...
%         anovaP, ALPHA, string(anovaP < ALPHA));
%
%     % raw baseline pool (simple, no split-half)
%     basePool = vertcat(s.trialBaselineVals{:});
%     basePool = basePool(~isnan(basePool));
%     baseMean = mean(basePool);
%     baseSD   = std(basePool);
%     fprintf('Raw baseline pool: n=%d, mean=%.5f, SD=%.5f\n', numel(basePool), baseMean, baseSD);
%
%     % preferred direction (from baseline-subtracted means, for direction selection only)
%     [~, prefDirIdx] = max(s.meanDirResponse);
%     meanPreferredRaw = mean(s.trialRawResp{prefDirIdx}, 'omitnan');
%     threshold = baseMean + NSD * baseSD;
%     fprintf('Preferred direction: %d, meanPreferredRaw = %.5f\n', prefDirIdx, meanPreferredRaw);
%     fprintf('Threshold (baseMean + %d*baseSD) = %.5f\n', NSD, threshold);
%     fprintf('Magnitude check: %.5f > %.5f -> %s\n', meanPreferredRaw, threshold, ...
%         string(meanPreferredRaw > threshold));
%
%     isResponsive = (anovaP < ALPHA) && (meanPreferredRaw > threshold);
%     fprintf('\n isResponsive = %s \n', string(isResponsive));
%
%     % --- quick plot: this bouton's per-direction RAW mean traces if available ---
%     if isfield(s, 'baselineSubtracted')
%         figure('Name', sprintf('Bouton %d diagnostic', b));
%         hold on;
%         cmap = parula(nDir);
%         for thisDir = 1:nDir
%             respTrace = mean(s.baselineSubtracted{thisDir}, 1, 'omitnan');
%             plot(respTrace, 'Color', cmap(thisDir,:), 'LineWidth', 1.2);
%         end
%         yline(threshold - baseMean, 'r--', 'LineWidth', 1.5, 'DisplayName', 'threshold (baseline-subtracted equivalent)');
%         title(sprintf('Bouton %d | anovaP=%.3g | prefRaw=%.4f | thresh=%.4f | %s', ...
%             b, anovaP, meanPreferredRaw, threshold, string(isResponsive)));
%         xlabel('Frames within response window'); ylabel('\DeltaF/F (baseline-subtracted)');
%     end
% end

%%  ANOVA-protected t-test / Wilcoxon vs pooled baseline
% Replaces the SD heuristic (above) with 
% test comparing the preferred direction's single trials directly
% against the pooled pre-stimulus baseline trials .

% ANOVA (already computed above, anovaP) is kept as a first pass:
% we only trust the direct preferred-vs-baseline comparison if the
%  ANOVA across all 8 directions was already significant. This
% guards against the selection-bias problem of picking whichever
% direction happens to look best and testing ONLY that one, uncorrected.

% Welch's t-test does not assume equal variance between groups.
% Wilcoxon rank-sum is nonparametric -- no normality assumption, more
% robust to skewed/heavy-tailed dF/F distributions.

ttestP   = nan(nBoutonsTotal, 1);
ranksumP = nan(nBoutonsTotal, 1);
isResponsive_ttest   = false(nBoutonsTotal, 1);
isResponsive_ranksum = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s = allDirTuning(b);
    if isnan(anovaP(b)) || isnan(prefDirIdx(b))
        continue;
    end

    prefTrials = s.trialRawResp{prefDirIdx(b)}(:);
    prefTrials = prefTrials(~isnan(prefTrials));

    baselineTrials = vertcat(s.trialBaselineVals{:});
%     baselineTrials = s.trialBaselineVals{prefDirIdx(b)}; % include only
%     matched baseline? 

    baselineTrials = baselineTrials(~isnan(baselineTrials));

    if numel(prefTrials) < 2 || numel(baselineTrials) < 2
        continue; % need at least 2 trials per group for a meaningful test
    end

    [~, ttestP(b)] = ttest2(prefTrials, baselineTrials, 'Vartype', 'unequal');
    ranksumP(b) = ranksum(prefTrials, baselineTrials);

    isResponsive_ttest(b)   = (anovaP(b) < alpha) && (ttestP(b) < alpha) && ...
                              (meanPreferredRaw(b) > baseMean(b));
    isResponsive_ranksum(b) = (anovaP(b) < alpha) && (ranksumP(b) < alpha) && ...
                              (meanPreferredRaw(b) > baseMean(b));
end

fprintf('\n%d / %d boutons responsive (ANOVA-protected Welch''s t-test, p<%.2f).\n', ...
    sum(isResponsive_ttest), nBoutonsTotal, alpha);
fprintf('%d / %d boutons responsive (ANOVA-protected Wilcoxon rank-sum, p<%.2f).\n', ...
    sum(isResponsive_ranksum), nBoutonsTotal, alpha);

for b = 1:nBoutonsTotal
    allDirTuning(b).ttestP               = ttestP(b);
    allDirTuning(b).ranksumP             = ranksumP(b);
    allDirTuning(b).isResponsive_ttest   = isResponsive_ttest(b);
    allDirTuning(b).isResponsive_ranksum = isResponsive_ranksum(b);
end

%% Cross-validated tuning R^2 (calc_kfold_R2), same pattern as DotFields_TuningCurveAnalysis_compareStates_2PData.m -- no running/stationary state-splitting here, just 8 directions.
% Asks how reliable is the shape of the tuning curve.. 
% Mirrors:
%   units(thisROI).allSpikes{s, istate} = traceAccumulator;
%   minTrial = min(min(cellfun(@(x) size(x,2), units(1).allSpikes)));
%   units(thisROI).allSpikesDownsample = cellfun(@(x) x(:,1:minTrial), units(thisROI).allSpikes, 'UniformOutput', false);
%   gca = units(thisROI).allSpikesDownsample(:,1)';
%   [units(thisROI).statR2, units(thisROI).statR2_pval] = calc_kfold_R2(gca, r2opts.kval, r2opts.nPerms, ...
%       r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);
%
% NOTE: DotFields script uses RAW (non-baseline-subtracted) per-trial
% values as input to calc_kfold_R2. This section defaults to
% trialMeanResp (BASELINE-SUBTRACTED) for consistency with the rest of
% this script and because direction-dependent baseline level (seen in
% bouton 167) could bias a cross-validated R^2 computed on raw values.
% Swap trialMeanResp -> trialRawResp below (one line, marked) to match
% the DotFields script exactly.
%


r2opts.kval       = 3;
r2opts.nPerms     = 10;
r2opts.randFlag   = 1;
r2opts.validMeans = 1;
r2opts.nShuffle   = 100;

cvR2   = nan(nBoutonsTotal, 1);
cvPval = nan(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s    = allDirTuning(b);
    nDir = numel(s.trialMeanResp);

    %1 x nDir cell array, each cell = 1 x nTrials row vector ---
    gca = cell(1, nDir);
    for thisDir = 1:nDir
        vals = s.trialRawResp{thisDir}(:)';   %% SWAP HERE for s.trialMeanResp{thisDir}(:)' to match DotFields script exactly
        vals = vals(~isnan(vals));
        gca{thisDir} = vals;
    end

    %  downsample to equal trial count per direction (trials across conditions can be different because of bad-frames realted exclusion)
    trialCounts = cellfun(@numel, gca);
    if any(trialCounts == 0)
        continue; % skip boutons with a direction that has zero valid trials
    end
    minTrial = min(trialCounts);
    gcaDownsampled = cellfun(@(x) x(1:minTrial), gca, 'UniformOutput', false);

    % call edds function to do cross-val explained variace on the shape of
    % the tuning curve
    [cvR2(b), cvPval(b)] = calc_kfold_R2(gcaDownsampled, r2opts.kval, r2opts.nPerms, ...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

    if mod(b, 100) == 0
        fprintf('Processed %d / %d boutons (cross-val R^2)...\n', b, nBoutonsTotal);
    end
end

isTunedCV = cvPval < alpha;
fprintf('\n%d / %d boutons show significant cross-validated tuning (p < %.2f).\n', ...
    sum(isTunedCV, 'omitnan'), nBoutonsTotal, alpha);

for b = 1:nBoutonsTotal
    allDirTuning(b).cvR2      = cvR2(b);
    allDirTuning(b).cvPval    = cvPval(b);
    allDirTuning(b).isTunedCV = isTunedCV(b);
end

%% Population-level comparison: SD-heuristic vs t-test vs Wilcoxon vs cross-validated R^2

validCompare = ~isnan(cvPval) & ~isnan(anovaP); % boutons where all methods produced a result

n_bothPass    = sum(isResponsive(validCompare)  & isTunedCV(validCompare));
n_onlyANOVA   = sum(isResponsive(validCompare)  & ~isTunedCV(validCompare));
n_onlyCV      = sum(~isResponsive(validCompare) & isTunedCV(validCompare));
n_neitherPass = sum(~isResponsive(validCompare) & ~isTunedCV(validCompare));
nCompared     = sum(validCompare);

fprintf('\n=== Population comparison: ANOVA+threshold vs cross-validated R^2 (n=%d boutons compared) ===\n', nCompared);
fprintf('%-30s %6d  (%.1f%%)\n', 'Both methods agree: responsive',     n_bothPass,    100*n_bothPass/nCompared);
fprintf('%-30s %6d  (%.1f%%)\n', 'ANOVA+thresh only',                  n_onlyANOVA,   100*n_onlyANOVA/nCompared);
fprintf('%-30s %6d  (%.1f%%)\n', 'Cross-val R^2 only',                 n_onlyCV,      100*n_onlyCV/nCompared);
fprintf('%-30s %6d  (%.1f%%)\n', 'Both agree: NOT responsive',         n_neitherPass, 100*n_neitherPass/nCompared);

fprintf('\nTotal passing ANOVA+threshold: %d / %d (%.1f%%)\n', ...
    sum(isResponsive(validCompare)), nCompared, 100*sum(isResponsive(validCompare))/nCompared);
fprintf('Total passing ANOVA-protected t-test: %d / %d (%.1f%%)\n', ...
    sum(isResponsive_ttest(validCompare)), nCompared, 100*sum(isResponsive_ttest(validCompare))/nCompared);
fprintf('Total passing ANOVA-protected Wilcoxon: %d / %d (%.1f%%)\n', ...
    sum(isResponsive_ranksum(validCompare)), nCompared, 100*sum(isResponsive_ranksum(validCompare))/nCompared);
fprintf('Total passing cross-val R^2:   %d / %d (%.1f%%)\n', ...
    sum(isTunedCV(validCompare)), nCompared, 100*sum(isTunedCV(validCompare))/nCompared);

% Scatter: -log10(anovaP) vs cvR2, colored by cvPval significance
figure('Position', [100 100 900 400]);

subplot(1,2,1);
scatter(-log10(anovaP(validCompare)), cvR2(validCompare), 15, cvPval(validCompare), 'filled');
hold on;
xline(-log10(alpha), 'k--', 'LineWidth', 1);
xlabel('-log_{10}(ANOVA p-value)'); ylabel('Cross-validated R^2');
cb = colorbar; ylabel(cb, 'cross-val p-value');
title('ANOVA significance vs cross-val R^2 (color = cv p-value)');

subplot(1,2,2);
categories = {'SD-heuristic', 'ANOVA+t-test', 'ANOVA+Wilcoxon', 'CV-R^2'};
counts = [sum(isResponsive(validCompare)), sum(isResponsive_ttest(validCompare)), ...
          sum(isResponsive_ranksum(validCompare)), sum(isTunedCV(validCompare))];
bar(counts);
xticklabels(categories);
ylabel('Number of boutons');
title(sprintf('Method comparison (n=%d)', nCompared));

sgtitle('DirTuning: SD-heuristic vs t-test vs Wilcoxon vs cross-validated R^2, all stats saved per-ROI on allDirTuning');


%%
allDirTuning = computeDirTuningOSI(allDirTuning);
allDirTuning = computeDirTuningDSI(allDirTuning);

%%
% plotExampleDirTuningBouton(allDirTuning, 1871)
% plotExampleDirTuningBouton(allDirTuning, 789)
% plotExampleDirTuningBouton(allDirTuning, 1838)


%% 
%DirTuning_PlotAllSelectedBoutons


%% is it necessary to shuffle?
allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning); % now also tests OSI_simple/DSI_simple significance
DirTuning_SummaryFigures