% DirTuning_IncTuningCurveAnalysis_StationaryOnly_2PData.m
%
% Stationary-only version of Figure5_3_DriftingGratings_basics.m /
% DirTuning_IncTuningCurveAnalysis_RunningOnly_2PData.m, mirroring the
% running-only script's field names and structure EXACTLY (so this
% struct plugs directly into computeDirTuningOSI.m, computeDirTuningDSI.m,
% and plotExampleDirTuningBouton.m without any changes needed) -- just
% classifying trials as STATIONARY instead of running, and sticking to
% RAW (non-baseline-subtracted) values throughout, per your last
% instruction (matching the compare-states script's convention, not the
% baseline-subtracted running-only convention).
%
% Trials are classified and filtered to STATIONARY ONLY at pooling time
% (per direction, per trial), following the wheel-based classification
% pattern in DotFields_IncTuningCurveAnalysis_compareStatesV2_2PData.m.
%
% ASSUMPTIONS ABOUT DATA SHAPE (same as the running-only script -- please
% verify against your actual response.mat):
%   response.wheelData(d) -- one entry PER DIRECTION, same order as
%                             psthData(d); .alignedResponses [nWheelSamples x nTrials],
%                             .timeVector [nWheelSamples x 1]; trial ti
%                             in wheelData(d) == trial ti in psthData(d).
%
% Requires calc_kfold_R2.m, computeDirTuningOSI.m, computeDirTuningDSI.m

%%
runSpeedThresh  = 3;    % used only to define the upper bound for "stationary" (see classification below)
statSpeedThresh = 0.5;  % mean speed below this = candidate "stationary"
propThresh      = 0.75;

alpha  = 0.05;
numSD  = 1;

r2_thresh  = 0.1;
r2p_thresh = 0.05;

r2opts.kval       = 3;
r2opts.nPerms     = 10;
r2opts.randFlag   = 1;
r2opts.validMeans = 1;
r2opts.nShuffle   = 100;
seed = 1;

nTrialsTarget = 5; % minimum trials/direction required to attempt cvR2 for a bouton -- check the trial-count summary printed below before finalizing this

%%
DirTuningTable = filterMasterTable_usingNameSessionPairs('MouseID', ...
    {'M25132','M25133', 'M26003'}, 'Exclude', 0, 'HasStimulus', {'DirTuning'});
allMice    = DirTuningTable.MouseID;
uniqueMice = unique(allMice, 'stable');

respWin_2sOn = [0.1 3];   % postStimTime == 4 (2s-on/2s-off)
respWin_1sOn = [0.1 3];   % postStimTime == 3 (1s-on/2s-off)

allDirTuning     = struct('fullTraceSub', {}, 'timeVec', {}, ...
                           'stimOnDuration', {}, 'meanDirResponse', {}, ...
                           'trialMeanResp', {}, 'trialRawResp', {}, 'trialBaselineVals', {}, ...
                           'nTrialsStatPerDir', {}, ...
                           'stimValues', {}, 'stimVariant', {}, 'sessionLabel', {});
sessionLabels    = {};
haveWarnedShape  = false;

sessionStatCounts = []; % [nSessions x nDir], for the trial-count summary after pooling

%% loop over mice and pool data (STATIONARY TRIALS ONLY)
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
            warning('  wheelData missing for %s -- cannot classify stationary trials -- skipping.', thisSessionName);
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
        baseIdx   = timeVec >= -0.75 & timeVec <= 0; % kept for trialBaselineVals only -- NOT subtracted from the trace (see below)
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

        %% classify stationary trials per direction (KEEP STATIONARY ONLY)
        statTrialIdx      = cell(nDir, 1);
        nTrialsStatPerDir = nan(1, nDir);

        for thisDir = 1:nDir
            grpWheel    = response.wheelData(thisDir);
            speedMatrix = grpWheel.alignedResponses;      % [nWheelSamples x nTrials]
            tVecWheel   = grpWheel.timeVector(:)';
            stimFramesMask = tVecWheel >= stimFramesMask_range(1) & tVecWheel <= stimFramesMask_range(2);

            nTrialsThisDir = size(speedMatrix, 2);
            thisDirStatFlag = false(nTrialsThisDir, 1);

            for ti = 1:nTrialsThisDir
                singleTrialTrace = speedMatrix(:, ti);
                if all(isnan(singleTrialTrace)), continue; end

                meanSpeed      = nanmean(singleTrialTrace(stimFramesMask));
                propStationary = sum(singleTrialTrace(stimFramesMask) < runSpeedThresh) / sum(stimFramesMask);

                if propStationary >= propThresh && meanSpeed < statSpeedThresh
                    thisDirStatFlag(ti) = true;
                end
            end

            statTrialIdx{thisDir}       = find(thisDirStatFlag);
            nTrialsStatPerDir(thisDir)  = numel(statTrialIdx{thisDir});
        end

        fprintf('  Stationary trials per direction: %s\n', mat2str(nTrialsStatPerDir));

        if any(nTrialsStatPerDir == 0)
            warning('  At least one direction has zero stationary trials for %s -- skipping session.', thisSessionName);
            continue;
        end

        if isempty(sessionStatCounts)
            sessionStatCounts = nTrialsStatPerDir;
        elseif size(sessionStatCounts, 2) == nDir
            sessionStatCounts = [sessionStatCounts; nTrialsStatPerDir]; %#ok<AGROW>
        else
            warning('  nDir (%d) differs from previous sessions -- skipping this session''s counts in the summary table (data itself is still pooled).', nDir);
        end

        %% per-bouton, STATIONARY TRIALS ONLY, RAW (no baseline subtraction)
        for iBouton = 1:nBoutons
            fullTraceSub      = cell(nDir, 1); % NOTE: name kept for compatibility with plotExampleDirTuningBouton.m/computeDirTuningOSI.m/DSI.m; contents are RAW (not baseline-subtracted)
            meanDirResponse   = nan(nDir, 1);
            trialMeanResp     = cell(nDir, 1); % identical to trialRawResp -- kept as a separate field for downstream compatibility
            trialRawResp      = cell(nDir, 1);
            trialBaselineVals = cell(nDir, 1);

            for thisDir = 1:nDir
                traceMat = squeeze(response.psthData(thisDir).alignedResponses(iBouton, :, :));
                if isvector(traceMat)
                    traceMat = traceMat(:);
                end
                traceMat = double(traceMat)'; % [nTrials x nTimepoints]

                traceMat = traceMat(statTrialIdx{thisDir}, :); % KEEP STATIONARY TRIALS ONLY

                % NO baseline subtraction -- RAW throughout, per your
                % instruction. perTrialBaseline still computed (kept for
                % reference/potential future use), just never subtracted.
                perTrialBaseline = mean(traceMat(:, baseIdx), 2, 'omitnan');

                fullTraceSub{thisDir}    = traceMat; % RAW
                meanDirResponse(thisDir) = mean(mean(traceMat(:, respIdx), 2, 'omitnan'), 'omitnan');

                trialMeanResp{thisDir}     = mean(traceMat(:, respIdx), 2, 'omitnan');
                trialRawResp{thisDir}      = mean(traceMat(:, respIdx), 2, 'omitnan');
                trialBaselineVals{thisDir} = perTrialBaseline;
            end

            allDirTuning(end+1) = struct( ...
                'fullTraceSub',       {fullTraceSub}, ...
                'timeVec',            timeVec, ...
                'stimOnDuration',     stimOnDuration, ...
                'meanDirResponse',    meanDirResponse, ...
                'trialMeanResp',      {trialMeanResp}, ...
                'trialRawResp',       {trialRawResp}, ...
                'trialBaselineVals',  {trialBaselineVals}, ...
                'nTrialsStatPerDir',  nTrialsStatPerDir, ...
                'stimValues',         stimVals(:), ...
                'stimVariant',        stimVariant, ...
                'sessionLabel',       sprintf('%s_%s', thisMouse, thisSessionName));
        end

        sessionLabels = [sessionLabels; repmat({sprintf('%s_%s', thisMouse, thisSessionName)}, nBoutons, 1)];
        fprintf('  Added %d boutons (running total: %d).\n', nBoutons, numel(allDirTuning));
    end
end

nBoutonsTotal = numel(allDirTuning);
fprintf('\nDone pooling. %d boutons total across %d mice (STATIONARY trials only, raw).\n', nBoutonsTotal, numel(uniqueMice));

%% report: how many trials per direction survive the stationary filter
if ~isempty(sessionStatCounts)
    fprintf('\n=== Stationary-trial counts per direction, across %d sessions ===\n', size(sessionStatCounts,1));
    fprintf('%-8s %8s %8s %8s\n', 'Dir#', 'Mean', 'Min', 'Max');
    for d = 1:size(sessionStatCounts, 2)
        fprintf('%-8d %8.1f %8d %8d\n', d, mean(sessionStatCounts(:,d)), min(sessionStatCounts(:,d)), max(sessionStatCounts(:,d)));
    end
    fprintf('\nMin stationary trials in ANY direction, ANY session: %d\n', min(sessionStatCounts(:)));
else
    fprintf('\nNo sessions produced valid stationary-trial counts to summarize.\n');
end

%% ANOVA-based responsiveness (stationary trials only)
anovaP           = nan(nBoutonsTotal, 1);
baseMean         = nan(nBoutonsTotal, 1);
baseSD           = nan(nBoutonsTotal, 1);
meanPreferredRaw = nan(nBoutonsTotal, 1);
prefDirIdx       = nan(nBoutonsTotal, 1);
isResponsive     = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s    = allDirTuning(b);
    nDir = numel(s.trialMeanResp);

    y = []; grp = [];
    for thisDir = 1:nDir
        vals = s.trialMeanResp{thisDir}(:);
        y    = [y;   vals];
        grp  = [grp; repmat(thisDir, numel(vals), 1)];
    end

    if numel(y) < nDir * 2 || all(isnan(y))
        continue;
    end

    validRows = ~isnan(y);
    anovaP(b) = anova1(y(validRows), grp(validRows), 'off');

    baselinePool = vertcat(s.trialBaselineVals{:});
    baselinePool = baselinePool(~isnan(baselinePool));
    baseMean(b)  = mean(baselinePool);
    baseSD(b)    = std(baselinePool);

    [~, prefDirIdx(b)] = max(s.meanDirResponse);
    meanPreferredRaw(b) = mean(s.trialRawResp{prefDirIdx(b)}, 'omitnan');

    isResponsive(b) = (anovaP(b) < alpha) && ...
                      (meanPreferredRaw(b) > baseMean(b) + numSD * baseSD(b));
end

fprintf('\n%d / %d boutons responsive (ANOVA[dirs]+threshold, p<%.2f, >%dSD above baseline, STATIONARY only).\n', ...
    sum(isResponsive), nBoutonsTotal, alpha, numSD);

for b = 1:nBoutonsTotal
    allDirTuning(b).anovaP_stat           = anovaP(b);
    allDirTuning(b).baseMean_stat         = baseMean(b);
    allDirTuning(b).baseSD_stat           = baseSD(b);
    allDirTuning(b).meanPreferredRaw_stat = meanPreferredRaw(b);
    allDirTuning(b).prefDirIdx_stat       = prefDirIdx(b);
    allDirTuning(b).isResponsive_stat     = isResponsive(b);
end

%% ANOVA-protected t-test / Wilcoxon vs pooled baseline (stationary trials only)
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
    baselineTrials = baselineTrials(~isnan(baselineTrials));

    if numel(prefTrials) < 2 || numel(baselineTrials) < 2
        continue;
    end

    [~, ttestP(b)] = ttest2(prefTrials, baselineTrials, 'Vartype', 'unequal');
    ranksumP(b) = ranksum(prefTrials, baselineTrials);

    isResponsive_ttest(b)   = (anovaP(b) < alpha) && (ttestP(b) < alpha) && (meanPreferredRaw(b) > baseMean(b));
    isResponsive_ranksum(b) = (anovaP(b) < alpha) && (ranksumP(b) < alpha) && (meanPreferredRaw(b) > baseMean(b));
end

fprintf('%d / %d boutons responsive (ANOVA-protected Welch t-test, p<%.2f, STATIONARY only).\n', sum(isResponsive_ttest), nBoutonsTotal, alpha);
fprintf('%d / %d boutons responsive (ANOVA-protected Wilcoxon, p<%.2f, STATIONARY only).\n', sum(isResponsive_ranksum), nBoutonsTotal, alpha);

for b = 1:nBoutonsTotal
    allDirTuning(b).ttestP_stat             = ttestP(b);
    allDirTuning(b).ranksumP_stat           = ranksumP(b);
    allDirTuning(b).isResponsive_ttest_stat   = isResponsive_ttest(b);
    allDirTuning(b).isResponsive_ranksum_stat = isResponsive_ranksum(b);
end

%% ============================================================
%  Cross-validated R^2 (calc_kfold_R2), STATIONARY TRIALS ONLY
% ============================================================
cvR2   = nan(nBoutonsTotal, 1);
cvPval = nan(nBoutonsTotal, 1);
minTrialUsed = nan(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s    = allDirTuning(b);
    nDir = numel(s.trialMeanResp);

    gca = cell(1, nDir);
    for thisDir = 1:nDir
        vals = s.trialRawResp{thisDir}(:)';
        vals = vals(~isnan(vals));
        gca{thisDir} = vals;
    end

    trialCounts = cellfun(@numel, gca);
    if any(trialCounts == 0)
        continue;
    end

    minTrial = min(trialCounts);
    if minTrial < nTrialsTarget
        continue;
    end

    gcaDownsampled = cellfun(@(x) x(1:minTrial), gca, 'UniformOutput', false);
    minTrialUsed(b) = minTrial;

    [cvR2(b), cvPval(b)] = calc_kfold_R2(gcaDownsampled, r2opts.kval, r2opts.nPerms, ...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle, seed);

    if mod(b, 100) == 0
        fprintf('Processed %d / %d boutons (cross-val R^2, stationary only)...\n', b, nBoutonsTotal);
    end
end

isTunedCVR2_stat = cvR2 > r2_thresh & cvPval < alpha;
fprintf('\n%d / %d boutons show significant cross-validated tuning (p<%.2f, R^2>%.2f, STATIONARY only).\n', ...
    sum(isTunedCVR2_stat, 'omitnan'), nBoutonsTotal, alpha, r2_thresh);

nAttempted = sum(~isnan(cvR2));
fprintf('%d / %d boutons had >= %d stationary trials/direction and were attempted.\n', ...
    nAttempted, nBoutonsTotal, nTrialsTarget);
if nAttempted > 0
    fprintf('Trials/direction actually used after downsampling -- mean: %.1f, min: %d, max: %d\n', ...
        mean(minTrialUsed, 'omitnan'), min(minTrialUsed), max(minTrialUsed));
end

for b = 1:nBoutonsTotal
    allDirTuning(b).cvR2_stat         = cvR2(b);
    allDirTuning(b).cvPval_stat       = cvPval(b);
    allDirTuning(b).isTunedCVR2_stat  = isTunedCVR2_stat(b);
    allDirTuning(b).minTrialUsed_stat = minTrialUsed(b);
end

validIdx_stat = find(cat(1, allDirTuning.cvR2_stat) > r2_thresh & cat(1, allDirTuning.cvPval_stat) < r2p_thresh);
fprintf('\n%d / %d boutons pass R^2 filter (R^2>%.2f, p<%.2f), STATIONARY only.\n', ...
    numel(validIdx_stat), nBoutonsTotal, r2_thresh, r2p_thresh);

%% (optional) OSI/DSI on stationary-only data
allDirTuning = computeDirTuningOSI(allDirTuning, 'isTunedCVR2_stat');
allDirTuning = computeDirTuningDSI(allDirTuning, 'isTunedCVR2_stat');


DirTuning_PlotAllSelectedBoutons_StationaryOnly