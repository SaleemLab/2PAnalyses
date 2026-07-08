% RFMapping_RunningOnly_SufficientCoverage.m
%
% Instead of comparing running vs stationary (which we found suffers
% from uneven/thin trial counts creating artifacts), this uses ONLY
% running trials, but first requires a session to have SUFFICIENT
% running trials at EVERY grid position (+ blank) before including any
% of its boutons. "Sufficient" is set empirically from the actual data
% distribution (e.g. the 33rd percentile of running-trial-counts across
% all session x position combinations), not an arbitrary fixed number.
%
% STEP 1: characterize the distribution of running-trial-counts per
%         position, across ALL sessions -- plot it, pick a percentile-
%         based threshold.
% STEP 2: a session qualifies only if EVERY position (+ blank) has at
%         least that many running trials. Sessions that don't qualify
%         are excluded ENTIRELY (not just the thin positions) -- since
%         analyseRFMapping.m's ANOVA needs all positions present for a
%         valid, comparable test.
% STEP 3: for qualifying sessions, recompute isResponsive using ONLY
%         running trials (identical logic to analyseRFMapping.m), and
%         report the responsive fraction among this "sufficient
%         running" pool.

%%
pairs = struct;
pairs.M25132 = {'20260219','20260223','20260226','20260228','20260303','20260313','20260306'};
pairs.M25133 = {'20260219','20260223','20260221'};
pairs.M26003 = {'20260316','20260322','20260324','20260325'};
pairs.M25132 = {'20260204A','20260204B','20260204C', ...
                '20260212A','20260212B','20260212C','20260212D','20260212E','20260212F', ...
                '20260214A','20260214B','20260214C','20260214D', ...
                '20260216A','20260216B','20260216C'};
pairs.M25133 = {'20260205A','20260205B','20260205C', ...
                '20260212A','20260212B','20260212C','20260212D','20260212E', ...
                '20260216A','20260216B','20260216C','20260216D'};
pairs.M26003 = {'20260307A','20260307B', ...
                '20260313A','20260313B','20260313C'};

stimName = 'RFMapping';
respWin  = [0.5 3];
baseWin  = [-1.0 0];

runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;

ALPHA = 0.05;
NSD   = 2;

percentileForThreshold = 33.33; %
%% 

filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'RFMapping'});
allMice    = filteredTable.MouseID;
uniqueMice = unique(allMice, 'stable');

%% ===================== STEP 1: characterize running-trial-count distribution per position =====================
sessionInfoCache = struct('sessionLabel', {}, 'thisMouse', {}, 'thisSessionName', {}, ...
    'runCountsByGroup', {}, 'groupIdx', {}, 'nROI', {});

allRunCounts = []; % pooled across all session x position combinations, for the distribution

for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));

    for iSess = 1:length(mouseSessIdx)
        tableRow        = filteredTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        sessionLabel    = sprintf('%s_%s', thisMouse, thisSessionName);

        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath), warning('sfi missing for %s -- skipping.', sessionLabel); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        iStim = find(contains(stimNames, stimName), 1);
        if isempty(iStim), warning('No %s file for %s -- skipping.', stimName, sessionLabel); continue; end

        try
            load(sessionFileInfo.stimFiles(iStim).Response, 'response');
        catch ME
            warning('Could not load response for %s: %s', sessionLabel, ME.message); continue;
        end
        if ~isfield(response, 'wheelData') || numel(response.wheelData) ~= numel(response.psthData)
            warning('wheelData missing/mismatched for %s -- skipping.', sessionLabel); continue;
        end

        psthData = response.psthData;
        nROI = size(psthData(1).alignedResponses, 1);
        wheelTimeVec = response.wheelData(1).timeVector(:)';
        wheelRespIdx = wheelTimeVec >= respWin(1) & wheelTimeVec <= respWin(2);

        nGroups = numel(psthData);
        runCountsByGroup = nan(nGroups, 1);
        for g = 1:nGroups
            wheelTrials = response.wheelData(g).alignedResponses;
            nRun = 0;
            for ti = 1:size(wheelTrials, 2)
                trace = wheelTrials(wheelRespIdx, ti);
                if all(isnan(trace)), continue; end
                meanSpeed   = nanmean(trace);
                propRunning = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    nRun = nRun + 1;
                end
            end
            runCountsByGroup(g) = nRun;
        end

        sessionInfoCache(end+1) = struct('sessionLabel', sessionLabel, 'thisMouse', thisMouse, ...
            'thisSessionName', thisSessionName, 'runCountsByGroup', runCountsByGroup, ...
            'groupIdx', (1:nGroups)', 'nROI', nROI); 

        allRunCounts = [allRunCounts; runCountsByGroup(:)]; 

        fprintf('%s: %d groups, running trial counts range %d-%d\n', ...
            sessionLabel, nGroups, min(runCountsByGroup), max(runCountsByGroup));
    end
end

%% ===================== plot the distribution + pick threshold =====================
minRunTrialThreshold = round(prctile(allRunCounts, percentileForThreshold));
fprintf('\n=== Running-trial-count distribution (pooled across %d session x position combinations) ===\n', numel(allRunCounts));
fprintf('Median: %.1f | Mean: %.1f | %gth percentile: %.1f --> using threshold = %d\n', ...
    median(allRunCounts), mean(allRunCounts), percentileForThreshold, ...
    prctile(allRunCounts, percentileForThreshold), minRunTrialThreshold);

figure('Position', [100 100 500 400]);
histogram(allRunCounts, 0:1:max(allRunCounts)+1);
hold on;
xline(minRunTrialThreshold, 'r--', 'LineWidth', 2, 'Label', sprintf('%gth pctile = %d', percentileForThreshold, minRunTrialThreshold));
xlabel('Running trials per position (pooled across all sessions)'); ylabel('Count');
title('Distribution used to set the minimum-trial threshold');

%% ===================== STEP 2: which sessions qualify (ALL positions >= threshold)? =====================
minRunTrialThreshold = 5;
qualifyingSessions = {};
for s = 1:numel(sessionInfoCache)
    if all(sessionInfoCache(s).runCountsByGroup >= minRunTrialThreshold)
        qualifyingSessions{end+1} = sessionInfoCache(s).sessionLabel;
    end
end

fprintf('\n%d / %d sessions qualify (ALL positions have >= %d running trials):\n', ...
    numel(qualifyingSessions), numel(sessionInfoCache), minRunTrialThreshold);
for s = 1:numel(qualifyingSessions)
    fprintf('  %s\n', qualifyingSessions{s});
end

if isempty(qualifyingSessions)
    error('No sessions qualify at this threshold -- consider a lower percentile or excluding blank from the requirement.');
end

%% ===================== STEP 3: recompute isResponsive using RUNNING TRIALS ONLY, qualifying sessions =====================
allRFRunningOnly = struct('mouseID', {}, 'sessionName', {}, 'sessionLabel', {}, 'roiIdx', {}, ...
    'pValANOVA', {}, 'isResponsive', {}, 'centerAz', {}, 'centerEl', {});

for s = 1:numel(sessionInfoCache)
    sc = sessionInfoCache(s);
    if ~ismember(sc.sessionLabel, qualifyingSessions), continue; end

    infoPath = findSessionFileInfoFilePath(sc.thisMouse, sc.thisSessionName);
    loadedInfo      = load(infoPath, 'sessionFileInfo');
    sessionFileInfo = loadedInfo.sessionFileInfo;
    stimNames       = {sessionFileInfo.stimFiles.name};
    iStim = find(contains(stimNames, stimName), 1);
    load(sessionFileInfo.stimFiles(iStim).Response, 'response');

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

    % classify running trials only (stationary/ambiguous just excluded)
    runFlagByGroup = cell(numel(psthData), 1);
    for g = 1:numel(psthData)
        wheelTrials = response.wheelData(g).alignedResponses;
        nTrialsHere = size(wheelTrials, 2);
        rf = false(nTrialsHere, 1);
        for ti = 1:nTrialsHere
            trace = wheelTrials(wheelRespIdx, ti);
            if all(isnan(trace)), continue; end
            meanSpeed   = nanmean(trace);
            propRunning = sum(trace > statSpeedThresh) / sum(wheelRespIdx);
            rf(ti) = (propRunning >= propThresh) && (meanSpeed > runSpeedThresh);
        end
        runFlagByGroup{g} = rf;
    end

    for iROI = 1:nROI
        bTrials = squeeze(psthData(blankIdx).alignedResponses(iROI, :, :));
        if isvector(bTrials), bTrials = bTrials(:); end
        bTrials = bTrials(:, runFlagByGroup{blankIdx});

        tB = mean(bTrials(baseIdx, :), 1, 'omitnan');
        bCorr = bTrials - tB;
        blankTrialMeans = mean(bCorr(respIdx, :), 1, 'omitnan')';

        allTrialMeans = blankTrialMeans(:);
        groupLabels   = repmat(numel(gridPSTHIdx)+1, numel(blankTrialMeans), 1);
        meanGridResponse = nan(numel(uEl_plot), numel(uAz));

        for pIdx = 1:numel(gridPSTHIdx)
            g = gridPSTHIdx(pIdx);
            trialsAtPos = squeeze(psthData(g).alignedResponses(iROI, :, :));
            if isvector(trialsAtPos), trialsAtPos = trialsAtPos(:); end
            trialsAtPos = trialsAtPos(:, runFlagByGroup{g});
            if isempty(trialsAtPos), continue; end

            tB2 = mean(trialsAtPos(baseIdx, :), 1, 'omitnan');
            corr2 = trialsAtPos - tB2;
            posMeans = mean(corr2(respIdx, :), 1, 'omitnan')';
            allTrialMeans = [allTrialMeans; posMeans(:)]; %#ok<AGROW>
            groupLabels   = [groupLabels; repmat(pIdx, numel(posMeans), 1)]; %#ok<AGROW>

            rowIdx = find(uEl_plot == gridStim(pIdx,2), 1);
            colIdx = find(uAz == gridStim(pIdx,1), 1);
            meanGridResponse(rowIdx, colIdx) = mean(posMeans, 'omitnan');
        end

        validIdx = ~isnan(allTrialMeans) & ~isnan(groupLabels);
        isResponsive = false; pValANOVA = NaN; centerAz = NaN; centerEl = NaN;
        if numel(unique(groupLabels(validIdx))) >= 2
            pValANOVA = anova1(allTrialMeans(validIdx), groupLabels(validIdx), 'off');
            bMean = mean(blankTrialMeans, 'omitnan');
            bStd  = std(blankTrialMeans, 'omitnan');
            [prefVal, mI] = max(meanGridResponse(:), [], 'omitnan');
            isResponsive = (pValANOVA < ALPHA) && (prefVal > (bMean + NSD*bStd));
            if isResponsive && ~isnan(mI)
                [rPeak, cPeak] = ind2sub(size(meanGridResponse), mI);
                centerAz = uAz(cPeak); centerEl = uEl_plot(rPeak);
            end
        end

        allRFRunningOnly(end+1) = struct('mouseID', sc.thisMouse, 'sessionName', sc.thisSessionName, ...
            'sessionLabel', sc.sessionLabel, 'roiIdx', iROI, 'pValANOVA', pValANOVA, ...
            'isResponsive', isResponsive, 'centerAz', centerAz, 'centerEl', centerEl); %#ok<SAGROW>
    end

    fprintf('%s: added %d boutons.\n', sc.sessionLabel, nROI);
end

nResp = sum([allRFRunningOnly.isResponsive]);
fprintf('\n=== FINAL: running-only, sufficient-coverage pool ===\n');
fprintf('%d qualifying sessions, %d total boutons.\n', numel(qualifyingSessions), numel(allRFRunningOnly));
fprintf('Responsive (running trials only): %d / %d (%.1f%%)\n', ...
    nResp, numel(allRFRunningOnly), 100*nResp/numel(allRFRunningOnly));