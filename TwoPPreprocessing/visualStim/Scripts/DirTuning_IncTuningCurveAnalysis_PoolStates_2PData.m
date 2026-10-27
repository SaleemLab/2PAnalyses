% DirTuning_IncTuningCurveAnalysis_poolstates_2PData.m
%
% Compare-states version of the DirTuning pooling script: classifies
% EVERY trial (per direction) as stationary or running, following the
% wheel-based classification pattern in
% DotFields_IncTuningCurveAnalysis_compareStatesV2_2PData.m, then pools
% stationary + running trials TOGETHER per direction and computes ONE
% combined cross-validated R^2 per bouton -- to test whether having more
% trials per direction (regardless of behavioral state) reveals more
% tuned boutons than either state alone would.
%
% Sticks to RAW (non-baseline-subtracted) values throughout, per your
% instruction -- meanDirResponse/trialRawResp/fullTrace are raw. If you
% want to run OSI/DSI on this later, remember OSI/DSI's ratio formula is
% NOT invariant to a raw baseline offset (see earlier discussion) --
% something to account for before comparing OSI/DSI computed on this
% raw stationary/running struct against the baseline-subtracted
% running-only struct.
%
% ASSUMPTIONS ABOUT DATA SHAPE (same as the running-only script --
% please verify against your actual response.mat):
%   response.wheelData(d) -- one entry PER DIRECTION, same order as
%                             psthData(d); .alignedResponses [nWheelSamples x nTrials],
%                             .timeVector [nWheelSamples x 1]; trial ti
%                             in wheelData(d) == trial ti in psthData(d).
%
% Requires calc_kfold_R2.m on your MATLAB path.

%%
runSpeedThresh  = 3;    % mean speed above this = candidate "running"
statSpeedThresh = 0.5;  % mean speed below this = candidate "stationary"
propThresh      = 0.75;

alpha  = 0.05;
numSD  = 1;             % ADDED: needed for ANOVA+threshold responsiveness criterion
r2_thresh  = 0.1;
r2p_thresh = 0.05;

r2opts.kval       = 3;
r2opts.nPerms     = 10;
r2opts.randFlag   = 1;
r2opts.validMeans = 1;
r2opts.nShuffle   = 100;
seed = 1;

nTrialsTarget = 5; % minimum trials/direction required (after pooling stat+run) to attempt cvR2 for a bouton -- ADJUST based on the trial-count distribution printed below

baseWin = [-0.75 0];   % ADDED: pre-stimulus baseline window

%%
DirTuningTable = filterMasterTable_usingNameSessionPairs('MouseID', ...
    {'M25132','M25133', 'M26003'}, 'Exclude', 0, 'HasStimulus', {'DirTuning'});
allMice    = DirTuningTable.MouseID;
uniqueMice = unique(allMice, 'stable');

respWin_2sOn = [0.1 3];   % postStimTime == 4 (2s-on/2s-off)
respWin_1sOn = [0.1 3];   % postStimTime == 3 (1s-on/2s-off)

allDirTuning     = struct('fullTrace_stat', {}, 'fullTrace_run', {}, 'timeVec', {}, ...
                           'stimOnDuration', {}, ...
                           'meanDirResponse_stat', {}, 'meanDirResponse_run', {}, ...
                           'trialRawResp_stat', {}, 'trialRawResp_run', {}, ...
                           'trialBaselineVals_stat', {}, 'trialBaselineVals_run', {}, ... % ADDED
                           'nTrialsStatPerDir', {}, 'nTrialsRunPerDir', {}, ...
                           'stimValues', {}, 'stimVariant', {}, 'sessionLabel', {});
sessionLabels    = {};
haveWarnedShape  = false;

% running tally of trial counts per direction across all sessions, for
% the reporting block after pooling
sessionStatCounts = [];
sessionRunCounts  = [];

%% loop over mice and pool data (STATIONARY vs RUNNING, per trial)
for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));
    fprintf('MOUSE: %s | %d sessions\n', thisMouse, length(mouseSessIdx));

    for iSess = 1:length(mouseSessIdx)
        tableRow        = DirTuningTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        fprintf('\n--- Session: %s ---\n', thisSessionName);

        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath)
            warning('sfi missing for %s -- skipping.', thisSessionName);
            continue;
        end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        DirTuneIdx = find(contains(stimNames, 'DirTuning'));
        if isempty(DirTuneIdx)
            warning('No DirTuning files for %s -- skipping.', thisSessionName);
            continue;
        end

        try
            load(sessionFileInfo.stimFiles(DirTuneIdx(1)).Response, 'response');
        catch ME
            warning('  Could not load response for %s: %s', thisSessionName, ME.message);
            continue;
        end
        if ~isfield(response, 'psthData') || isempty(response.psthData)
            warning('  psthData empty/missing for %s -- skipping.', thisSessionName);
            continue;
        end
        if ~isfield(response, 'wheelData') || isempty(response.wheelData)
            warning('  wheelData missing for %s -- cannot classify trials -- skipping.', thisSessionName);
            continue;
        end

        switch response.postStimTime
            case 4
                respWin = respWin_2sOn; stimVariant = 4; stimOnDuration = 2;
            case 3
                respWin = respWin_1sOn; stimVariant = 3; stimOnDuration = 1;
            otherwise
                warning('  Unexpected postStimTime (%.2f) for %s -- skipping.', ...
                    response.postStimTime, thisSessionName);
                continue;
        end

        stimFramesMask_range = [-0.2, stimOnDuration + 0.8];

        nDir      = numel(response.psthData);
        timeVec   = response.psthData(1).timeVector(:)';
        respIdx   = timeVec >= respWin(1) & timeVec <= respWin(2);
        baseIdx   = timeVec >= baseWin(1) & timeVec <= baseWin(2);   % ADDED
        stimVals  = arrayfun(@(s) s.stimValue, response.psthData);

        nBoutons = size(response.psthData(1).alignedResponses, 1);
        if ~haveWarnedShape
            fprintf('  Detected shape: %d boutons x %d timepoints x %d trials (dir 1)\n', ...
                size(response.psthData(1).alignedResponses,1), ...
                size(response.psthData(1).alignedResponses,2), ...
                size(response.psthData(1).alignedResponses,3));
            haveWarnedShape = true;
        end

        %% classify EVERY trial per direction as stationary (0) / running (1) / ambiguous (excluded)
        statTrialIdx = cell(nDir, 1);
        runTrialIdx  = cell(nDir, 1);
        nTrialsStatPerDir = nan(1, nDir);
        nTrialsRunPerDir  = nan(1, nDir);

                for thisDir = 1:nDir
            grpWheel    = response.wheelData(thisDir);
            speedMatrix = grpWheel.alignedResponses;      % [nWheelSamples x nTrials]
            tVecWheel   = grpWheel.timeVector(:)';
            stimFramesMask = tVecWheel >= stimFramesMask_range(1) & tVecWheel <= stimFramesMask_range(2);

            nTrialsThisDir = size(speedMatrix, 2);
            thisDirFlag = nan(nTrialsThisDir, 1); % NaN = ambiguous/excluded, 0 = stat, 1 = run

            for ti = 1:nTrialsThisDir
                singleTrialTrace = speedMatrix(:, ti);
                if all(isnan(singleTrialTrace)), continue; end

                meanSpeed      = nanmean(singleTrialTrace(stimFramesMask));
                propRunning    = sum(singleTrialTrace(stimFramesMask) > statSpeedThresh) / sum(stimFramesMask);
                propStationary = sum(singleTrialTrace(stimFramesMask) < runSpeedThresh)  / sum(stimFramesMask);

                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    thisDirFlag(ti) = 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    thisDirFlag(ti) = 0;
                end
                % else: ambiguous/ITI-crossing trial -- left as NaN, excluded from BOTH states
            end

            statTrialIdx{thisDir}     = find(thisDirFlag == 0);
            runTrialIdx{thisDir}      = find(thisDirFlag == 1);
            nTrialsStatPerDir(thisDir) = numel(statTrialIdx{thisDir});
            nTrialsRunPerDir(thisDir)  = numel(runTrialIdx{thisDir});
                end

                fprintf('  Stationary trials per direction: %s\n', mat2str(nTrialsStatPerDir));
                fprintf('  Running trials per direction:    %s\n', mat2str(nTrialsRunPerDir));

                if any(nTrialsStatPerDir == 0) || any(nTrialsRunPerDir == 0)
                    warning('  At least one direction has zero trials in one state for %s -- skipping session.', thisSessionName);
                    continue;
                end

                if isempty(sessionStatCounts)
                    sessionStatCounts = nTrialsStatPerDir;
                    sessionRunCounts  = nTrialsRunPerDir;
                elseif size(sessionStatCounts, 2) == nDir
                    sessionStatCounts = [sessionStatCounts; nTrialsStatPerDir]; %#ok<AGROW>
                    sessionRunCounts  = [sessionRunCounts; nTrialsRunPerDir]; %#ok<AGROW>
                else
                    warning('  nDir (%d) differs from previous sessions -- skipping this session''s counts in the summary table.', nDir);
                end

                %% per-bouton, RAW (no baseline subtraction), split by state
                for iBouton = 1:nBoutons
                    fullTrace_stat = cell(nDir, 1);
                    fullTrace_run  = cell(nDir, 1);
                    meanDirResponse_stat = nan(nDir, 1);
                    meanDirResponse_run  = nan(nDir, 1);
                    trialRawResp_stat = cell(nDir, 1);
                    trialRawResp_run  = cell(nDir, 1);
                    trialBaselineVals_stat = cell(nDir, 1);  % ADDED
                    trialBaselineVals_run  = cell(nDir, 1);  % ADDED

                    for thisDir = 1:nDir
                        traceMat = squeeze(response.psthData(thisDir).alignedResponses(iBouton, :, :));
                        if isvector(traceMat)
                            traceMat = traceMat(:);
                        end
                        traceMat = double(traceMat)'; % [nTrials x nTimepoints]

                        traceMatStat = traceMat(statTrialIdx{thisDir}, :);
                        traceMatRun  = traceMat(runTrialIdx{thisDir}, :);

                        fullTrace_stat{thisDir} = traceMatStat; % RAW, no baseline subtraction
                        fullTrace_run{thisDir}  = traceMatRun;  % RAW, no baseline subtraction

                        meanDirResponse_stat(thisDir) = mean(mean(traceMatStat(:, respIdx), 2, 'omitnan'), 'omitnan');
                        meanDirResponse_run(thisDir)  = mean(mean(traceMatRun(:, respIdx), 2, 'omitnan'), 'omitnan');

                        trialRawResp_stat{thisDir} = mean(traceMatStat(:, respIdx), 2, 'omitnan');
                        trialRawResp_run{thisDir}  = mean(traceMatRun(:, respIdx), 2, 'omitnan');

                        trialBaselineVals_stat{thisDir} = mean(traceMatStat(:, baseIdx), 2, 'omitnan');  % ADDED
                        trialBaselineVals_run{thisDir}  = mean(traceMatRun(:, baseIdx), 2, 'omitnan');   % ADDED
                    end

                    allDirTuning(end+1) = struct( ...
                        'fullTrace_stat',       {fullTrace_stat}, ...
                        'fullTrace_run',        {fullTrace_run}, ...
                        'timeVec',              timeVec, ...
                        'stimOnDuration',       stimOnDuration, ...
                        'meanDirResponse_stat', meanDirResponse_stat, ...
                        'meanDirResponse_run',  meanDirResponse_run, ...
                        'trialRawResp_stat',    {trialRawResp_stat}, ...
                        'trialRawResp_run',     {trialRawResp_run}, ...
                        'trialBaselineVals_stat', {trialBaselineVals_stat}, ...  % ADDED
                        'trialBaselineVals_run',  {trialBaselineVals_run}, ...   % ADDED
                        'nTrialsStatPerDir',    nTrialsStatPerDir, ...
                        'nTrialsRunPerDir',     nTrialsRunPerDir, ...
                        'stimValues',           stimVals(:), ...
                        'stimVariant',          stimVariant, ...
                        'sessionLabel',         sprintf('%s_%s', thisMouse, thisSessionName));
                end

                sessionLabels = [sessionLabels; repmat({sprintf('%s_%s', thisMouse, thisSessionName)}, nBoutons, 1)];
                fprintf('  Added %d boutons (running total: %d).\n', nBoutons, numel(allDirTuning));
    end
end

nBoutonsTotal = numel(allDirTuning);
fprintf('\nDone pooling. %d boutons total across %d mice (STATIONARY vs RUNNING, raw).\n', nBoutonsTotal, numel(uniqueMice));

%% report: trial counts per direction, per state
if ~isempty(sessionStatCounts)
    fprintf('\n=== Trial counts per direction, across %d sessions ===\n', size(sessionStatCounts,1));
    fprintf('%-8s %10s %10s %10s | %10s %10s %10s\n', 'Dir#', 'Stat Mean', 'Stat Min', 'Stat Max', 'Run Mean', 'Run Min', 'Run Max');
    for d = 1:size(sessionStatCounts, 2)
        fprintf('%-8d %10.1f %10d %10d | %10.1f %10d %10d\n', d, ...
            mean(sessionStatCounts(:,d)), min(sessionStatCounts(:,d)), max(sessionStatCounts(:,d)), ...
            mean(sessionRunCounts(:,d)),  min(sessionRunCounts(:,d)),  max(sessionRunCounts(:,d)));
    end
    fprintf('\nMin stationary trials in ANY direction, ANY session: %d\n', min(sessionStatCounts(:)));
    fprintf('Min running trials in ANY direction, ANY session:    %d\n', min(sessionRunCounts(:)));
else
    fprintf('\nNo sessions produced valid trial counts to summarize.\n');
end

%% ANOVA-based responsiveness (all states combined: stationary + running pooled)  % ADDED SECTION

anovaP_comb           = nan(nBoutonsTotal, 1);
baseMean_comb         = nan(nBoutonsTotal, 1);
baseSD_comb           = nan(nBoutonsTotal, 1);
meanPreferredRaw_comb = nan(nBoutonsTotal, 1);
prefDirIdx_comb       = nan(nBoutonsTotal, 1);
isResponsive_comb     = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s    = allDirTuning(b);
    nDir = numel(s.trialRawResp_stat);

    y = []; grp = [];
    for thisDir = 1:nDir
        valsStat = s.trialRawResp_stat{thisDir}(:);
        valsRun  = s.trialRawResp_run{thisDir}(:);
        vals     = [valsStat; valsRun];
        y   = [y;   vals];
        grp = [grp; repmat(thisDir, numel(vals), 1)];
    end

    if numel(y) < nDir * 2 || all(isnan(y))
        continue;
    end

    validRows = ~isnan(y);
    anovaP_comb(b) = anova1(y(validRows), grp(validRows), 'off');

    baselinePool = [vertcat(s.trialBaselineVals_stat{:}); vertcat(s.trialBaselineVals_run{:})];
    baselinePool = baselinePool(~isnan(baselinePool));
    baseMean_comb(b) = mean(baselinePool);
    baseSD_comb(b)   = std(baselinePool);

    meanDirResponse_comb = nan(nDir, 1);
    for thisDir = 1:nDir
        valsStat = s.trialRawResp_stat{thisDir}(:);
        valsRun  = s.trialRawResp_run{thisDir}(:);
        meanDirResponse_comb(thisDir) = mean([valsStat; valsRun], 'omitnan');
    end

    [~, prefDirIdx_comb(b)] = max(meanDirResponse_comb);

    prefValsStat = s.trialRawResp_stat{prefDirIdx_comb(b)}(:);
    prefValsRun  = s.trialRawResp_run{prefDirIdx_comb(b)}(:);
    meanPreferredRaw_comb(b) = mean([prefValsStat; prefValsRun], 'omitnan');

    isResponsive_comb(b) = (anovaP_comb(b) < alpha) && ...
        (meanPreferredRaw_comb(b) > baseMean_comb(b) + numSD * baseSD_comb(b));
end

fprintf('\n%d / %d boutons responsive (ANOVA[dirs]+threshold, p<%.2f, >%dSD above baseline, ALL STATES COMBINED).\n', ...
    sum(isResponsive_comb), nBoutonsTotal, alpha, numSD);

for b = 1:nBoutonsTotal
    allDirTuning(b).anovaP_comb           = anovaP_comb(b);
    allDirTuning(b).baseMean_comb         = baseMean_comb(b);
    allDirTuning(b).baseSD_comb           = baseSD_comb(b);
    allDirTuning(b).meanPreferredRaw_comb = meanPreferredRaw_comb(b);
    allDirTuning(b).prefDirIdx_comb       = prefDirIdx_comb(b);
    allDirTuning(b).isResponsive_comb     = isResponsive_comb(b);
end

%% ANOVA-protected t-test / Wilcoxon vs pooled baseline (all states combined)  % ADDED SECTION
ttestP_comb   = nan(nBoutonsTotal, 1);
ranksumP_comb = nan(nBoutonsTotal, 1);
isResponsive_ttest_comb   = false(nBoutonsTotal, 1);
isResponsive_ranksum_comb = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s = allDirTuning(b);
    if isnan(anovaP_comb(b)) || isnan(prefDirIdx_comb(b))
        continue;
    end

    prefValsStat = s.trialRawResp_stat{prefDirIdx_comb(b)}(:);
    prefValsRun  = s.trialRawResp_run{prefDirIdx_comb(b)}(:);
    prefTrials   = [prefValsStat; prefValsRun];
    prefTrials   = prefTrials(~isnan(prefTrials));

    baselineTrials = [vertcat(s.trialBaselineVals_stat{:}); vertcat(s.trialBaselineVals_run{:})];
    baselineTrials = baselineTrials(~isnan(baselineTrials));

    if numel(prefTrials) < 2 || numel(baselineTrials) < 2
        continue;
    end

    [~, ttestP_comb(b)] = ttest2(prefTrials, baselineTrials, 'Vartype', 'unequal');
    ranksumP_comb(b) = ranksum(prefTrials, baselineTrials);

    isResponsive_ttest_comb(b)   = (anovaP_comb(b) < alpha) && (ttestP_comb(b) < alpha) && (meanPreferredRaw_comb(b) > baseMean_comb(b));
    isResponsive_ranksum_comb(b) = (anovaP_comb(b) < alpha) && (ranksumP_comb(b) < alpha) && (meanPreferredRaw_comb(b) > baseMean_comb(b));
end

fprintf('%d / %d boutons responsive (ANOVA-protected Welch t-test, p<%.2f, ALL STATES COMBINED).\n', sum(isResponsive_ttest_comb), nBoutonsTotal, alpha);
fprintf('%d / %d boutons responsive (ANOVA-protected Wilcoxon, p<%.2f, ALL STATES COMBINED).\n', sum(isResponsive_ranksum_comb), nBoutonsTotal, alpha);

for b = 1:nBoutonsTotal
    allDirTuning(b).ttestP_comb             = ttestP_comb(b);
    allDirTuning(b).ranksumP_comb           = ranksumP_comb(b);
    allDirTuning(b).isResponsive_ttest_comb   = isResponsive_ttest_comb(b);
    allDirTuning(b).isResponsive_ranksum_comb = isResponsive_ranksum_comb(b);
end

%% ============================================================
%  Cross-validated R^2 (calc_kfold_R2), ALL STATES COMBINED
%  -- pools stationary + running trials together per direction (both
%  are RAW, so this is a straightforward concatenation), to test
%  whether having MORE trials per direction (regardless of behavioral
%  state) rescues/reveals more tuned boutons than either state alone.
% ============================================================
cvR2_comb   = nan(nBoutonsTotal, 1); cvPval_comb = nan(nBoutonsTotal, 1); minTrialUsed_comb = nan(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s    = allDirTuning(b);
    nDir = numel(s.trialRawResp_stat);

    gcaComb = cell(1, nDir);
    for thisDir = 1:nDir
        valsStat = s.trialRawResp_stat{thisDir}(:)';
        valsRun  = s.trialRawResp_run{thisDir}(:)';
        valsComb = [valsStat(~isnan(valsStat)), valsRun(~isnan(valsRun))];
        gcaComb{thisDir} = valsComb;
    end

    trialCountsComb = cellfun(@numel, gcaComb);
    if any(trialCountsComb == 0)
        continue; % a direction with zero trials (either state) for this bouton -- skip
    end

    minTrialComb = min(trialCountsComb);
    if minTrialComb < nTrialsTarget
        continue; % not enough trials in the sparsest direction -- skip rather than downsample too aggressively
    end

    gcaCombDown = cellfun(@(x) x(1:minTrialComb), gcaComb, 'UniformOutput', false);
    minTrialUsed_comb(b) = minTrialComb;
    [cvR2_comb(b), cvPval_comb(b)] = calc_kfold_R2(gcaCombDown, r2opts.kval, r2opts.nPerms, ...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle, seed);

    if mod(b, 100) == 0
        fprintf('Processed %d / %d boutons (cross-val R^2, all states combined)...\n', b, nBoutonsTotal);
    end
end

isTunedCVR2_comb = cvR2_comb > r2_thresh & cvPval_comb < alpha;

for b = 1:nBoutonsTotal
    allDirTuning(b).cvR2_comb         = cvR2_comb(b);
    allDirTuning(b).cvPval_comb       = cvPval_comb(b);
    allDirTuning(b).isTunedCVR2_comb  = isTunedCVR2_comb(b);
    allDirTuning(b).minTrialUsed_comb = minTrialUsed_comb(b);
end

%% ============================================================
%  Summary: combined (all states pooled) cvR2
% ============================================================
nAttemptedComb = sum(~isnan(cvR2_comb));

fprintf('\n=== Cross-validated R^2, all states combined ===\n');
fprintf('%d / %d boutons attempted (>= %d trials/direction); %d / %d tuned (R^2>%.2f, p<%.2f).\n', ...
    nAttemptedComb, nBoutonsTotal, nTrialsTarget, sum(isTunedCVR2_comb, 'omitnan'), nBoutonsTotal, r2_thresh, alpha);

validIdx_comb = find(isTunedCVR2_comb);
fprintf('\nTrials/direction actually used after downsampling -- mean: %.1f, min: %d, max: %d\n', ...
    mean(minTrialUsed_comb, 'omitnan'), min(minTrialUsed_comb(~isnan(minTrialUsed_comb))), max(minTrialUsed_comb(~isnan(minTrialUsed_comb))));
%DirTuning_PlotAllSelectedBoutons_Combined


%% Session and mouse distribution of tuned boutons (all states combined)

validIdx_comb = find(isTunedCVR2_comb);
sessionLabelsAll = {allDirTuning.sessionLabel};

tunedSessions = sessionLabelsAll(validIdx_comb);
[uniqueTunedSessions, ~, ic] = unique(tunedSessions);
countsPerSession = accumarray(ic, 1);

[countsSorted, sortIdx] = sort(countsPerSession, 'descend');
sessionsSorted = uniqueTunedSessions(sortIdx);

fprintf('=== Session distribution of isTunedCVR2_comb boutons ===\n');
fprintf('%d tuned boutons come from %d distinct sessions (out of %d sessions total in the dataset).\n\n', ...
    numel(validIdx_comb), numel(uniqueTunedSessions), numel(unique(sessionLabelsAll)));

fprintf('%-30s %10s %10s\n', 'Session', 'N tuned', 'Percent');
for i = 1:numel(sessionsSorted)
    fprintf('%-30s %10d %9.1f%%\n', sessionsSorted{i}, countsSorted(i), 100*countsSorted(i)/numel(validIdx_comb));
end

fprintf('\n%-30s %10s %10s %12s\n', 'Session', 'N tuned', 'N total', 'Tuned rate');
for i = 1:numel(sessionsSorted)
    nTotalThisSession = sum(strcmp(sessionLabelsAll, sessionsSorted{i}));
    fprintf('%-30s %10d %10d %11.1f%%\n', sessionsSorted{i}, countsSorted(i), nTotalThisSession, ...
        100*countsSorted(i)/nTotalThisSession);
end

% Mouse-level breakdown (extract mouse ID from sessionLabel, format "MouseID_SessionName")
mouseIDsAll = cellfun(@(x) strtok(x, '_'), sessionLabelsAll, 'UniformOutput', false);
tunedMice   = mouseIDsAll(validIdx_comb);
[uniqueTunedMice, ~, icMouse] = unique(tunedMice);
countsPerMouse = accumarray(icMouse, 1);
[countsMouseSorted, sortIdxMouse] = sort(countsPerMouse, 'descend');
miceSorted = uniqueTunedMice(sortIdxMouse);

fprintf('\n=== Mouse distribution of isTunedCVR2_comb boutons ===\n');
fprintf('%d tuned boutons come from %d distinct mice (out of %d mice total in the dataset).\n\n', ...
    numel(validIdx_comb), numel(uniqueTunedMice), numel(unique(mouseIDsAll)));

fprintf('%-15s %10s %10s\n', 'Mouse', 'N tuned', 'Percent');
for i = 1:numel(miceSorted)
    fprintf('%-15s %10d %9.1f%%\n', miceSorted{i}, countsMouseSorted(i), 100*countsMouseSorted(i)/numel(validIdx_comb));
end


%% Trials per direction, stationary vs running, pooled across sessions

fprintf('\n=== Trials per direction, pooled across %d sessions ===\n', size(sessionStatCounts,1));

fprintf('Stationary: mean=%.1f, SD=%.1f (across all sessions x directions)\n', ...
    mean(sessionStatCounts(:)), std(sessionStatCounts(:)));

fprintf('Running:    mean=%.1f, SD=%.1f (across all sessions x directions)\n', ...
    mean(sessionRunCounts(:)), std(sessionRunCounts(:)));

% Combined (stat+run summed per session x direction), for the pooled analysis
sessionCombCounts = sessionStatCounts + sessionRunCounts;
fprintf('Combined:   mean=%.1f, SD=%.1f (across all sessions x directions)\n', ...
    mean(sessionCombCounts(:)), std(sessionCombCounts(:)));


%% Proportion of trials running vs stationary, pooled across sessions

totalStat = sum(sessionStatCounts(:));
totalRun  = sum(sessionRunCounts(:));
totalTrials = totalStat + totalRun;

propRunning = 100 * totalRun / totalTrials;
propStationary = 100 * totalStat / totalTrials;

fprintf('\n=== Proportion of trials running vs stationary (pooled across sessions) ===\n');
fprintf('Stationary: %d / %d trials (%.1f%%)\n', totalStat, totalTrials, propStationary);
fprintf('Running:    %d / %d trials (%.1f%%)\n', totalRun, totalTrials, propRunning);

% Per-session breakdown, in case it varies a lot session to session
fprintf('\n=== Per-session breakdown ===\n');
for s = 1:size(sessionStatCounts, 1)
    sessStat = sum(sessionStatCounts(s,:));
    sessRun  = sum(sessionRunCounts(s,:));
    sessTotal = sessStat + sessRun;
    fprintf('Session %d: %.1f%% running (%d/%d trials)\n', s, 100*sessRun/sessTotal, sessRun, sessTotal);
end

%%
%% Add a combined (pooled stat+run) meanDirResponse field, needed by computeDirTuningOSI/DSI

for b = 1:numel(allDirTuning)
    s = allDirTuning(b);
    nDir = numel(s.meanDirResponse_stat);
    meanDirResponse_comb = nan(nDir, 1);

    for thisDir = 1:nDir
        valsStat = s.trialRawResp_stat{thisDir}(:);
        valsRun  = s.trialRawResp_run{thisDir}(:);
        meanDirResponse_comb(thisDir) = mean([valsStat; valsRun], 'omitnan');
    end

    allDirTuning(b).meanDirResponse = meanDirResponse_comb;
end

%% Build combined fullTraceSub field for plotting (pooled states)

for b = 1:numel(allDirTuning)
    s = allDirTuning(b);
    nDir = numel(s.fullTrace_stat);
    fullTraceSub = cell(nDir, 1);

    for thisDir = 1:nDir
        traceStat = s.fullTrace_stat{thisDir}; % [nTrials_stat x nTimepoints]
        traceRun  = s.fullTrace_run{thisDir};  % [nTrials_run x nTimepoints]
        fullTraceSub{thisDir} = [traceStat; traceRun]; % concatenate along trials
    end

    allDirTuning(b).fullTraceSub = fullTraceSub;
end

fprintf('Added combined fullTraceSub field to %d boutons.\n', numel(allDirTuning));