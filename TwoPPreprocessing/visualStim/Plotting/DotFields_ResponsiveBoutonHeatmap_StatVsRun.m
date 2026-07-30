% DotFields_ResponsiveBoutonHeatmap_StatVsRun.m
%
% Builds a "many individual boutons" heatmap panel, analogous to Fig 3.5a
% (running-in-darkness) and Horrocks et al. Fig 1d: each ROW is one
% bouton's trial-averaged PSTH (at that bouton's own preferred speed),
% NOT one trial. Stationary and running heatmaps are shown side by side
% using the SAME row order (same boutons, same sort), so shape
% differences between states are directly visible per-bouton.
%
% Population included: boutons responsive in BOTH stationary and running
% (SD-heuristic, isResponsive_stat & isResponsive_run), matching the
% "reliable in both states" logic used in Horrocks et al. Fig 1d.
%
% Requires: allDotUnits already built by
% DotFields_TuningCurveAnalysis_compareStatesV2_2PData.m (needs
% mouseID, sessionName, roiIdx, isResponsive_stat/run, prefSpeedIdx_stat/run,
% gaussChar_stat/run fields).
%
% Sorting options (set sortMethod below):
%   'peakLatency'   -- sort rows by time-of-peak in the STATIONARY trace
%   'classification' -- sort rows by gaussChar_stat (High/Band/Low/Trough),
%                       then by peak latency within each class

%%
sortMethod = 'classification';  % 'peakLatency' or 'classification'
respWin_range = [0.1 3];        % same final window as main pipeline

stimFramesMask_range = [0 2.0];
runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;

classNames  = {'High-pass','Band-pass','Low-pass','Trough'};
classColors = {[0.7 0.75 0.2], [0.45 0.5 0.15], [0.4 0.7 0.85], [0.15 0.35 0.55]}; % loosely matches Fig 3.5 palette

%% identify boutons responsive in BOTH states
bothIdx = find([allDotUnits.isResponsive_stat] & [allDotUnits.isResponsive_run]);
fprintf('%d boutons responsive in BOTH stationary and running (SD-heuristic).\n', numel(bothIdx));

nBoth = numel(bothIdx);
statTraces = cell(nBoth, 1);
runTraces  = cell(nBoth, 1);
peakLatencyStat = nan(nBoth, 1);
gaussCharStat   = nan(nBoth, 1);
timeVecRef = [];

%% group by session so raw data is loaded once per session, not once per bouton
sessionLabels_both = {allDotUnits(bothIdx).sessionLabel};
uniqueSessions = unique(sessionLabels_both, 'stable');

for iSess = 1:numel(uniqueSessions)
    thisSessionLabel = uniqueSessions{iSess};
    boutonPosInBoth = find(strcmp(sessionLabels_both, thisSessionLabel));
    if isempty(boutonPosInBoth), continue; end

    % get mouse/session names from the first matching unit
    exampleUnit = allDotUnits(bothIdx(boutonPosInBoth(1)));
    thisMouse = exampleUnit.mouseID;
    thisSessionName = exampleUnit.sessionName;

    fprintf('Loading %s (%d boutons)...\n', thisSessionLabel, numel(boutonPosInBoth));

    %% load + classify trials (same logic as main pipeline)
    infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
    if ~isfile(infoPath), warning('sfi missing for %s -- skipping.', thisSessionLabel); continue; end
    loadedInfo      = load(infoPath, 'sessionFileInfo');
    sessionFileInfo = loadedInfo.sessionFileInfo;
    stimNames       = {sessionFileInfo.stimFiles.name};
    dotIdx = find(contains(stimNames, 'DotMotion_SpeedTuning'), 1);
    if isempty(dotIdx), warning('No DotMotion_SpeedTuning file for %s -- skipping.', thisSessionLabel); continue; end

    load(sessionFileInfo.stimFiles(dotIdx).Response, 'response');
    load(sessionFileInfo.stimFiles(dotIdx).BonsaiData, 'bonsaiData'); %#ok<NASGU>

    nGroups = numel(response.wheelData);
    trialsSpeed2D = struct('VelX1', {}, 'numDots1', {}, 'runFlag', {}, 'origGroup', {}, 'origTrialInGroup', {});
    trialCounter = 1;
    for g = 1:nGroups
        grpWheel  = response.wheelData(g);
        speedMatrix = grpWheel.alignedResponses;
        tVecWheel   = grpWheel.timeVector;
        stimFramesMask = (tVecWheel >= stimFramesMask_range(1) & tVecWheel <= stimFramesMask_range(2));
        for ti = 1:size(speedMatrix, 2)
            singleTrialTrace = speedMatrix(:, ti);
            if all(isnan(singleTrialTrace)), continue; end
            meanSpeed      = nanmean(singleTrialTrace(stimFramesMask));
            propRunning    = sum(singleTrialTrace(stimFramesMask) > statSpeedThresh) / sum(stimFramesMask);
            propStationary = sum(singleTrialTrace(stimFramesMask) < runSpeedThresh)  / sum(stimFramesMask);
            runFlag = NaN;
            if propRunning >= propThresh && meanSpeed > runSpeedThresh
                runFlag = 1;
            elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                runFlag = 0;
            end
            if isnan(runFlag), continue; end
            trialsSpeed2D(trialCounter).VelX1 = grpWheel.stimValue;
            trialsSpeed2D(trialCounter).numDots1 = (grpWheel.stimValue ~= 1) * 573;
            trialsSpeed2D(trialCounter).runFlag          = runFlag;
            trialsSpeed2D(trialCounter).origGroup        = g;
            trialsSpeed2D(trialCounter).origTrialInGroup  = ti;
            trialCounter = trialCounter + 1;
        end
    end
    tsd = trialsSpeed2D;
    temp_tsd = tsd([tsd.numDots1] == 573);
    uniqueVelocities = unique(abs([temp_tsd.VelX1]));
    timeVec = response.psthData(1).timeVector(:)';
    if isempty(timeVecRef), timeVecRef = timeVec; end

    %% for each bouton in this session, pull PSTH at its OWN preferred speed, per state
    for k = 1:numel(boutonPosInBoth)
        rowIdx = boutonPosInBoth(k); % position within bothIdx / statTraces / runTraces
        unit = allDotUnits(bothIdx(rowIdx));
        thisROI = unit.roiIdx;

        prefSpeedIdx_stat = unit.prefSpeedIdx_stat;
        prefSpeedIdx_run  = unit.prefSpeedIdx_run;

        % stationary trace at its own preferred speed
        matchStat = find(abs([temp_tsd.VelX1]) == uniqueVelocities(prefSpeedIdx_stat) & [temp_tsd.runFlag] == 0);
        traceMatStat = nan(numel(matchStat), numel(timeVec));
        for mt = 1:numel(matchStat)
            og = temp_tsd(matchStat(mt)).origGroup; ot = temp_tsd(matchStat(mt)).origTrialInGroup;
            traceMatStat(mt,:) = squeeze(response.psthData(og).alignedResponses(thisROI, :, ot));
        end
        statTraces{rowIdx} = mean(traceMatStat, 1, 'omitnan');

        % running trace at its own preferred speed
        matchRun = find(abs([temp_tsd.VelX1]) == uniqueVelocities(prefSpeedIdx_run) & [temp_tsd.runFlag] == 1);
        traceMatRun = nan(numel(matchRun), numel(timeVec));
        for mt = 1:numel(matchRun)
            og = temp_tsd(matchRun(mt)).origGroup; ot = temp_tsd(matchRun(mt)).origTrialInGroup;
            traceMatRun(mt,:) = squeeze(response.psthData(og).alignedResponses(thisROI, :, ot));
        end
        runTraces{rowIdx} = mean(traceMatRun, 1, 'omitnan');

        gaussCharStat(rowIdx) = unit.gaussChar_stat;
        respIdxMask = timeVec >= respWin_range(1) & timeVec <= respWin_range(2);
        [~, pkI] = max(statTraces{rowIdx}(respIdxMask));
        respTimeVec = timeVec(respIdxMask);
        peakLatencyStat(rowIdx) = respTimeVec(pkI);
    end
end

%% drop any boutons that failed to load
validRows = ~cellfun(@isempty, statTraces) & ~cellfun(@isempty, runTraces);
statTraces = statTraces(validRows);
runTraces  = runTraces(validRows);
peakLatencyStat = peakLatencyStat(validRows);
gaussCharStat = gaussCharStat(validRows);

statMat = cat(1, statTraces{:});
runMat  = cat(1, runTraces{:});

%% normalize each row (min-max), same convention as Fig 3.5a / Horrocks Fig 1d
normRow = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);
statMatNorm = cell2mat(arrayfun(@(i) normRow(statMat(i,:)), (1:size(statMat,1))', 'UniformOutput', false));
runMatNorm  = cell2mat(arrayfun(@(i) normRow(runMat(i,:)),  (1:size(runMat,1))',  'UniformOutput', false));

%% sort rows
switch sortMethod
    case 'peakLatency'
        [~, sortOrder] = sort(peakLatencyStat);
    case 'classification'
        % sort by class (1-4), then by peak latency within class
        sortTable = table(gaussCharStat, peakLatencyStat, (1:numel(gaussCharStat))', ...
            'VariableNames', {'class','latency','origIdx'});
        sortTable = sortrows(sortTable, {'class','latency'});
        sortOrder = sortTable.origIdx;
    otherwise
        error('sortMethod must be ''peakLatency'' or ''classification''.');
end

statMatSorted = statMatNorm(sortOrder, :);
runMatSorted  = runMatNorm(sortOrder, :);
classSorted   = gaussCharStat(sortOrder);

%% Fig: side-by-side heatmaps, shared row order
figure('Color','w','Position',[100 50 1000 800]);

axStat = subplot(1,2,1);
imagesc(timeVecRef, 1:size(statMatSorted,1), statMatSorted, [0 1]);
colormap(axStat, flipud(gray));
xlabel('Time (s)'); ylabel(sprintf('%d boutons', size(statMatSorted,1)));
title('Stationary');
hold on; xline(0,'r:','LineWidth',1); xline(2,'r:','LineWidth',1);

axRun = subplot(1,2,2);
imagesc(timeVecRef, 1:size(runMatSorted,1), runMatSorted, [0 1]);
colormap(axRun, flipud(gray));
xlabel('Time (s)');
title('Running');
hold on; xline(0,'r:','LineWidth',1); xline(2,'r:','LineWidth',1);
cb = colorbar; cb.Label.String = 'norm. \DeltaF/F';

sgtitle(sprintf('Responsive boutons (both states, n=%d), sorted by %s', ...
    size(statMatSorted,1), sortMethod), 'FontWeight','bold');

%% if sorted by classification, add class boundary lines + labels
if strcmp(sortMethod, 'classification')
    boundaries = find(diff(classSorted) ~= 0);
    for b = boundaries'
        yline(axStat, b+0.5, 'b-', 'LineWidth', 1);
        yline(axRun,  b+0.5, 'b-', 'LineWidth', 1);
    end
    % class label ticks at the midpoint of each block
    classStarts = [1; boundaries+1];
    classEnds   = [boundaries; numel(classSorted)];
    for c = 1:numel(classStarts)
        midY = mean([classStarts(c), classEnds(c)]);
        thisClass = classSorted(classStarts(c));
        if ~isnan(thisClass) && thisClass >= 1 && thisClass <= 4
            text(axStat, -1.8, midY, classNames{thisClass}, 'FontSize', 8, ...
                'HorizontalAlignment','right', 'Color', classColors{thisClass}, 'FontWeight','bold');
        end
    end
end
