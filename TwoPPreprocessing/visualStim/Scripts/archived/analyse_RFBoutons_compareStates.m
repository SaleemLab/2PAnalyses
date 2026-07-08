% RFMapping_PooledActivityVsRunning.m
%
% Pools the activity-vs-running check across ALL RF mapping
% sessions/mice (same loop structure as the other pooling scripts).
% For each bouton, pools per-trial (runSpeed, activity) pairs across ALL
% grid positions + blank, computes a correlation, and classifies each
% trial as running/stationary (same thresholds as the DotFields
% scripts) so you can see, per session, how much running data is even
% available before deciding whether "not enough running" could explain
% low RF responsiveness rates.
%
% NOTE ON WINDOWS: uses the SAME respWin=[0.5,3], baseWin=[-1,0] as
% analyseRFMapping.m for activity, and evaluates running state over that
% SAME window (not a separate stimFramesMask like DotFields uses) --
% since RF mapping doesn't have a separately documented stim-on duration
% in the metadata used here. Flagging this as an assumption worth
% checking if you have the actual RF stimulus timing available.
%
% Requires response.wheelData to exist and correspond 1:1 with
% response.psthData (same check as the single-session version) --
% sessions failing this check are skipped with a warning, not silently
% guessed.

%% 
stimName  = 'RFMapping';

respWin = [0.5 3];
baseWin = [-1.0 0];

runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;
%% ===================================================

pairs=struct;
pairs.M25132 = {'20260219','20260223','20260226','20260228','20260303','20260313','20260306'};
pairs.M25133 = {'20260219','20260223','20260221'};
pairs.M26003 = {'20260316','20260322','20260324','20260325'};
filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'RFMapping'});

allMice    = filteredTable.MouseID;
uniqueMice = unique(allMice, 'stable');

allRFRunningCheck = struct('mouseID', {}, 'sessionName', {}, 'sessionLabel', {}, 'roiIdx', {}, ...
    'runSpeedTrials', {}, 'activityTrials', {}, 'runFlagTrials', {}, 'corrR', {}, 'corrP', {});

sessionSummary = struct('sessionLabel', {}, 'nTrialsTotal', {}, 'nRunning', {}, 'nStationary', {}, 'nAmbiguous', {});
haveWarnedShape = false;

for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));
    fprintf('MOUSE: %s | %d sessions\n', thisMouse, length(mouseSessIdx));

    for iSess = 1:length(mouseSessIdx)
        tableRow        = filteredTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        sessionLabel    = sprintf('%s_%s', thisMouse, thisSessionName);
        fprintf('  --- Session: %s ---\n', thisSessionName);

        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath), warning('    sfi missing -- skipping.'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        iStim = find(contains(stimNames,stimName), 1);
        if isempty(iStim), warning('    No %s file -- skipping.', stimName); continue; end

        try
            load(sessionFileInfo.stimFiles(iStim).Response, 'response');
        catch ME
            warning('    Could not load response: %s', ME.message); continue;
        end

        if ~isfield(response, 'wheelData') || isempty(response.wheelData)
            warning('    No response.wheelData for %s -- skipping (cannot check running).', sessionLabel);
            continue;
        end
        if numel(response.wheelData) ~= numel(response.psthData)
            warning('    wheelData (%d groups) does not match psthData (%d groups) for %s -- skipping.', ...
                numel(response.wheelData), numel(response.psthData), sessionLabel);
            continue;
        end

        psthData = response.psthData;
        nGroups  = numel(psthData);
        nROI     = size(psthData(1).alignedResponses, 1);

        timeVec = psthData(1).timeVector(:)';
        respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
        baseIdx = timeVec >= baseWin(1) & timeVec < baseWin(2);

        wheelTimeVec = response.wheelData(1).timeVector(:)';
        wheelRespIdx = wheelTimeVec >= respWin(1) & wheelTimeVec <= respWin(2);

        if ~haveWarnedShape
            fprintf('    Detected shape: %d ROIs x %d timepoints x %d trials (group 1)\n', ...
                size(psthData(1).alignedResponses,1), size(psthData(1).alignedResponses,2), ...
                size(psthData(1).alignedResponses,3));
            haveWarnedShape = true;
        end

        %% session-level running/stationary trial-count summary (pooled across ROIs, since running state is a trial/session property) ---
        nTrialsTotal = 0; nRunning = 0; nStationary = 0; nAmbiguous = 0;
        for g = 1:nGroups
            wheelTrials = response.wheelData(g).alignedResponses; % [nTimepoints x nTrials]
            for ti = 1:size(wheelTrials, 2)
                trace = wheelTrials(wheelRespIdx, ti);
                if all(isnan(trace)), continue; end
                meanSpeed      = nanmean(trace);
                propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
                propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdx);

                nTrialsTotal = nTrialsTotal + 1;
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    nRunning = nRunning + 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    nStationary = nStationary + 1;
                else
                    nAmbiguous = nAmbiguous + 1;
                end
            end
        end

        sessionSummary(end+1) = struct('sessionLabel', sessionLabel, 'nTrialsTotal', nTrialsTotal, ...
            'nRunning', nRunning, 'nStationary', nStationary, 'nAmbiguous', nAmbiguous); %#ok<SAGROW>

        fprintf('    Trials: %d total | %d running (%.1f%%) | %d stationary (%.1f%%) | %d ambiguous\n', ...
            nTrialsTotal, nRunning, 100*nRunning/max(nTrialsTotal,1), ...
            nStationary, 100*nStationary/max(nTrialsTotal,1), nAmbiguous);

        %%  per-bouton pooled activity vs running (all groups + blank, all trials) 
        for iROI = 1:nROI
            allActivity = []; allRunSpeed = []; allRunFlag = [];

            for g = 1:nGroups
                trialsThisGroup = squeeze(psthData(g).alignedResponses(iROI, :, :));
                if isvector(trialsThisGroup), trialsThisGroup = trialsThisGroup(:); end
                nTrialsHere = size(trialsThisGroup, 2);

                trialBaselines = mean(trialsThisGroup(baseIdx, :), 1, 'omitnan');
                correctedTrials = trialsThisGroup - trialBaselines;
                activityVals = mean(correctedTrials(respIdx, :), 1, 'omitnan')';

                wheelTrials = response.wheelData(g).alignedResponses;
                if size(wheelTrials, 2) ~= nTrialsHere
                    continue; % mismatched trial count for this group -- skip silently here (already warned at session level if pervasive)
                end
                runSpeedVals = mean(wheelTrials(wheelRespIdx, :), 1, 'omitnan')';

                runFlagVals = nan(nTrialsHere, 1);
                for ti = 1:nTrialsHere
                    trace = wheelTrials(wheelRespIdx, ti);
                    if all(isnan(trace)), continue; end
                    meanSpeed      = nanmean(trace);
                    propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
                    propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdx);
                    if propRunning >= propThresh && meanSpeed > runSpeedThresh
                        runFlagVals(ti) = 1;
                    elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                        runFlagVals(ti) = 0;
                    end
                end

                allActivity = [allActivity; activityVals];
                allRunSpeed = [allRunSpeed; runSpeedVals]; 
                allRunFlag  = [allRunFlag; runFlagVals];    
            end

            validRows = ~isnan(allActivity) & ~isnan(allRunSpeed);
            corrR = NaN; corrP = NaN;
            if sum(validRows) >= 3
                [corrR, corrP] = corr(allRunSpeed(validRows), allActivity(validRows), 'rows', 'complete');
            end

            allRFRunningCheck(end+1) = struct( ...
                'mouseID',        thisMouse, ...
                'sessionName',    thisSessionName, ...
                'sessionLabel',   sessionLabel, ...
                'roiIdx',         iROI, ...
                'runSpeedTrials', allRunSpeed, ...
                'activityTrials', allActivity, ...
                'runFlagTrials',  allRunFlag, ...
                'corrR',          corrR, ...
                'corrP',          corrP);
        end

        fprintf('    Added %d boutons (running total: %d).\n', nROI, numel(allRFRunningCheck));
    end
end

fprintf('\nDone. %d boutons pooled across %d sessions.\n', numel(allRFRunningCheck), numel(sessionSummary));

%% ===================== session-level summary table =====================
fprintf('\n=== Session summary: trial counts by running state ===\n');
fprintf('%-30s %10s %10s %10s %10s\n', 'Session', 'Total', 'Running', 'Stationary', 'Ambiguous');
for s = 1:numel(sessionSummary)
    ss = sessionSummary(s);
    fprintf('%-30s %10d %10d %10d %10d\n', ss.sessionLabel, ss.nTrialsTotal, ss.nRunning, ss.nStationary, ss.nAmbiguous);
end

%% ===================== population summary: activity-vs-running correlation =====================
validCorr = ~isnan([allRFRunningCheck.corrR]);
fprintf('\n%d / %d boutons have a valid activity-vs-running correlation.\n', sum(validCorr), numel(allRFRunningCheck));
fprintf('%d boutons show a significant positive correlation (r>0, p<0.05).\n', ...
    sum(validCorr & [allRFRunningCheck.corrR] > 0 & [allRFRunningCheck.corrP] < 0.05));

figure('Position', [100 100 500 400]);
histogram([allRFRunningCheck(validCorr).corrR], 30);
xlabel('Correlation (activity vs running speed)'); ylabel('Number of boutons');
title(sprintf('Activity-vs-running correlation distribution (n=%d)', sum(validCorr)));


%%
% RFMapping_BehaviorSplit.m
%
% Splits RF mapping trials by behavioral state (stationary/running) and
% recomputes analyseRFMapping.m's exact responsiveness logic (ANOVA
% across grid positions + blank, then preferred-position mean vs
% blankMean + 2*blankSD) SEPARATELY for each state.
%
% Starts with a HAND-PICKED subset of sessions that had a reasonable
% mix of both states (from the trial-count summary already run) --
% edit sessionPairs below to add/remove sessions.
%
% Uses the SAME fixed windows as analyseRFMapping.m (respWin=[0.5,3],
% baseWin=[-1,0]) and the SAME running/stationary thresholds as the
% DotFields scripts (propThresh=0.75, runSpeedThresh=3,
% statSpeedThresh=0.5), evaluated over respWin (see note in
% RFMapping_PooledActivityVsRunning.m about this being an assumption
% worth checking against actual RF stimulus timing).

essionPairs = struct();
sessionPairs.M25132 = {'20260228', '20260306'};  % 94/59 and 75/79 run/stat -- good mixes
sessionPairs.M26003 = {'20260316', '20260324'};  % 30/121 and 124/17 -- more skewed but usable
% Edit/add more sessions here as you like. Sessions with very few trials
% in one state (e.g. M25133's sessions, ~2-14 running trials) are
% deliberately excluded for now -- not enough data to say anything
% meaningful about that state.
 
stimName = 'RFMapping';
respWin  = [0.5 3];
baseWin  = [-1.0 0];
 
runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;
 
ALPHA = 0.05;
NSD   = 2;
%% ===================================================
 
mouseNames = fieldnames(sessionPairs);
allRFBehaviorSplit = struct('mouseID', {}, 'sessionName', {}, 'sessionLabel', {}, 'roiIdx', {}, ...
    'pValANOVA_stat', {}, 'isResponsive_stat', {}, 'centerAz_stat', {}, 'centerEl_stat', {}, ...
    'pValANOVA_run',  {}, 'isResponsive_run',  {}, 'centerAz_run',  {}, 'centerEl_run',  {}, ...
    'pValANOVA_combined', {}, 'isResponsive_combined', {}, 'centerAz_combined', {}, 'centerEl_combined', {});
 
for iMouse = 1:numel(mouseNames)
    thisMouse = mouseNames{iMouse};
    sessList  = sessionPairs.(thisMouse);
    fprintf('MOUSE: %s | %d hand-picked sessions\n', thisMouse, numel(sessList));
 
    for iSess = 1:numel(sessList)
        thisSessionName = sessList{iSess};
        sessionLabel = sprintf('%s_%s', thisMouse, thisSessionName);
        fprintf('  --- Session: %s ---\n', thisSessionName);
 
        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath), warning('    sfi missing -- skipping.'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        iStim = find(contains(stimNames, stimName), 1);
        if isempty(iStim), warning('    No %s file -- skipping.', stimName); continue; end
 
        try
            load(sessionFileInfo.stimFiles(iStim).Response, 'response');
        catch ME
            warning('    Could not load response: %s', ME.message); continue;
        end
 
        if ~isfield(response, 'wheelData') || isempty(response.wheelData)
            warning('    No wheelData for %s -- skipping.', sessionLabel); continue;
        end
        if numel(response.wheelData) ~= numel(response.psthData)
            warning('    wheelData/psthData group mismatch for %s -- skipping.', sessionLabel); continue;
        end
 
        psthData = response.psthData;
        stimVs   = vertcat(psthData.stimValue);
        nROI     = size(psthData(1).alignedResponses, 1);
 
        blankIdx  = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
        gridMask  = stimVs(:,1) ~= 200;
        gridPSTHIdx = find(gridMask);
        gridStim  = stimVs(gridMask, :);
 
        uAz      = sort(unique(gridStim(:,1)), 'ascend');
        uEl_plot = sort(unique(gridStim(:,2)), 'descend');
 
        timeVec = psthData(1).timeVector(:)';
        respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
        baseIdx = timeVec >= baseWin(1) & timeVec < baseWin(2);
 
        wheelTimeVec = response.wheelData(1).timeVector(:)';
        wheelRespIdx = wheelTimeVec >= respWin(1) & wheelTimeVec <= respWin(2);
 
        %% --- classify EVERY trial (grid + blank) as running/stationary, per group ---
        runFlagByGroup = cell(numel(psthData), 1);
        for g = 1:numel(psthData)
            wheelTrials = response.wheelData(g).alignedResponses;
            nTrialsHere = size(wheelTrials, 2);
            rf = nan(nTrialsHere, 1);
            for ti = 1:nTrialsHere
                trace = wheelTrials(wheelRespIdx, ti);
                if all(isnan(trace)), continue; end
                meanSpeed      = nanmean(trace);
                propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
                propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdx);
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    rf(ti) = 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    rf(ti) = 0;
                end
            end
            runFlagByGroup{g} = rf;
        end
 
        %% --- per-bouton, per-state: replicate analyseRFMapping.m's ANOVA+threshold logic ---
        for iROI = 1:nROI
            results = struct();
            for si = 1:3 % 1 = stationary, 2 = running, 3 = ALL TRIALS COMBINED (no state filter, matches original analyseRFMapping.m)
                if si <= 2
                    stateVal = si - 1;
                end
 
                % blank trials for this pass
                bTrials = squeeze(psthData(blankIdx).alignedResponses(iROI, :, :));
                if isvector(bTrials), bTrials = bTrials(:); end
                if si <= 2
                    blankRunFlag = runFlagByGroup{blankIdx};
                    blankStateMask = (blankRunFlag == stateVal);
                    bTrialsState = bTrials(:, blankStateMask);
                else
                    bTrialsState = bTrials; % ALL trials, no filtering -- matches original script exactly
                end
 
                if isempty(bTrialsState) || size(bTrialsState, 2) < 2
                    results(si).pValANOVA = NaN; results(si).isResponsive = false;
                    results(si).centerAz = NaN; results(si).centerEl = NaN;
                    continue;
                end
 
                trialBaselinesB = mean(bTrialsState(baseIdx, :), 1, 'omitnan');
                bTrialsCorrected = bTrialsState - trialBaselinesB;
                blankTrialMeans = mean(bTrialsCorrected(respIdx, :), 1, 'omitnan')';
 
                allTrialMeans = blankTrialMeans(:);
                groupLabels   = repmat(numel(gridPSTHIdx) + 1, numel(blankTrialMeans), 1); % blank = last group
 
                meanGridResponse = nan(numel(uEl_plot), numel(uAz));
 
                for pIdx = 1:numel(gridPSTHIdx)
                    g = gridPSTHIdx(pIdx);
                    trialsAtPos = squeeze(psthData(g).alignedResponses(iROI, :, :));
                    if isvector(trialsAtPos), trialsAtPos = trialsAtPos(:); end
                    if si <= 2
                        posRunFlag = runFlagByGroup{g};
                        posStateMask = (posRunFlag == stateVal);
                        trialsAtPosState = trialsAtPos(:, posStateMask);
                    else
                        trialsAtPosState = trialsAtPos; % ALL trials, no filtering
                    end
 
                    if isempty(trialsAtPosState), continue; end
 
                    trialBaselines = mean(trialsAtPosState(baseIdx, :), 1, 'omitnan');
                    correctedTrials = trialsAtPosState - trialBaselines;
                    posTrialMeans = mean(correctedTrials(respIdx, :), 1, 'omitnan')';
 
                    allTrialMeans = [allTrialMeans; posTrialMeans(:)]; %#ok<AGROW>
                    groupLabels   = [groupLabels; repmat(pIdx, numel(posTrialMeans), 1)]; %#ok<AGROW>
 
                    rowIdx = find(uEl_plot == gridStim(pIdx, 2), 1);
                    colIdx = find(uAz == gridStim(pIdx, 1), 1);
                    meanGridResponse(rowIdx, colIdx) = mean(posTrialMeans, 'omitnan');
                end
 
                validIdx = ~isnan(allTrialMeans) & ~isnan(groupLabels);
                allTrialMeans = allTrialMeans(validIdx);
                groupLabels   = groupLabels(validIdx);
 
                if isempty(allTrialMeans) || numel(unique(groupLabels)) < 2
                    results(si).pValANOVA = NaN; results(si).isResponsive = false;
                    results(si).centerAz = NaN; results(si).centerEl = NaN;
                    continue;
                end
 
                pValANOVA = anova1(allTrialMeans, groupLabels, 'off');
                blankMean = mean(blankTrialMeans, 'omitnan');
                blankStd  = std(blankTrialMeans, 'omitnan');
                [prefVal, mI] = max(meanGridResponse(:), [], 'omitnan');
 
                isResponsive = (pValANOVA < ALPHA) && (prefVal > (blankMean + NSD * blankStd));
 
                results(si).pValANOVA = pValANOVA;
                results(si).isResponsive = isResponsive;
                if isResponsive && ~isnan(mI)
                    [rPeak, cPeak] = ind2sub(size(meanGridResponse), mI);
                    results(si).centerAz = uAz(cPeak);
                    results(si).centerEl = uEl_plot(rPeak);
                else
                    results(si).centerAz = NaN; results(si).centerEl = NaN;
                end
            end
 
            allRFBehaviorSplit(end+1) = struct( ...
                'mouseID', thisMouse, 'sessionName', thisSessionName, 'sessionLabel', sessionLabel, 'roiIdx', iROI, ...
                'pValANOVA_stat', results(1).pValANOVA, 'isResponsive_stat', results(1).isResponsive, ...
                'centerAz_stat', results(1).centerAz, 'centerEl_stat', results(1).centerEl, ...
                'pValANOVA_run', results(2).pValANOVA, 'isResponsive_run', results(2).isResponsive, ...
                'centerAz_run', results(2).centerAz, 'centerEl_run', results(2).centerEl, ...
                'pValANOVA_combined', results(3).pValANOVA, 'isResponsive_combined', results(3).isResponsive, ...
                'centerAz_combined', results(3).centerAz, 'centerEl_combined', results(3).centerEl); %#ok<SAGROW>
        end
 
        fprintf('    Added %d boutons (running total: %d).\n', nROI, numel(allRFBehaviorSplit));
    end
end
 
fprintf('\nDone. %d boutons pooled across %d hand-picked sessions.\n', ...
    numel(allRFBehaviorSplit), sum(structfun(@numel, sessionPairs)));
 
nRespStat = sum([allRFBehaviorSplit.isResponsive_stat]);
nRespRun  = sum([allRFBehaviorSplit.isResponsive_run]);
nRespCombined = sum([allRFBehaviorSplit.isResponsive_combined]);
fprintf('\nResponsive (stationary): %d / %d (%.1f%%)\n', nRespStat, numel(allRFBehaviorSplit), 100*nRespStat/numel(allRFBehaviorSplit));
fprintf('Responsive (running):    %d / %d (%.1f%%)\n', nRespRun, numel(allRFBehaviorSplit), 100*nRespRun/numel(allRFBehaviorSplit));
fprintf('Responsive (ALL TRIALS COMBINED, no state split): %d / %d (%.1f%%)\n', ...
    nRespCombined, numel(allRFBehaviorSplit), 100*nRespCombined/numel(allRFBehaviorSplit));
 
nRespEither = sum([allRFBehaviorSplit.isResponsive_stat] | [allRFBehaviorSplit.isResponsive_run]);
nRespOnlyRun  = sum(~[allRFBehaviorSplit.isResponsive_stat] & [allRFBehaviorSplit.isResponsive_run]);
nRespOnlyStat = sum([allRFBehaviorSplit.isResponsive_stat] & ~[allRFBehaviorSplit.isResponsive_run]);
fprintf('Responsive in EITHER state (split): %d (%.1f%%)\n', nRespEither, 100*nRespEither/numel(allRFBehaviorSplit));
fprintf('Responsive ONLY during running:    %d\n', nRespOnlyRun);
fprintf('Responsive ONLY during stationary: %d\n', nRespOnlyStat);
 
%%  combined vs split comparison 

RFMapping_SelectSplitOnlyExamples
RFMapping_CheckPositionCoverageByState