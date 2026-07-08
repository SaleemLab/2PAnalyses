function plotRFExampleBoutonStateSplit(thisMouse, thisSessionName, thisROI)
% PLOTRFEXAMPLEBOUTONSTATESPLIT  Plots one bouton's RF heatmap+traces
% (using the existing plotRFHeatmapWithTraces function) separately for
% stationary trials, running trials, and all trials combined -- side by
% side, for direct visual comparison.
%
%   plotRFExampleBoutonStateSplit(thisMouse, thisSessionName, thisROI)
%
% Reloads raw session data and reclassifies trials by running/stationary
% state using the SAME thresholds as RFMapping_BehaviorSplit.m. For each
% of the three conditions (stationary / running / combined), builds a
% boutonData struct with the fields plotRFHeatmapWithTraces expects
% (meanGridResponse, meanTemporalResponse, meanBlankResponse,
% peakAmplitude, centerAz, centerEl, isResponsive), then calls that
% existing function to draw each panel -- no heatmap-plotting logic is
% reimplemented here.
%
% Requires plotRFHeatmapWithTraces.m on your MATLAB path.

%% ===================== CONFIG (keep in sync with RFMapping_BehaviorSplit.m) =====================
stimName = 'RFMapping';
respWin  = [0.5 3];
baseWin  = [-1.0 0];

runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;

ALPHA = 0.05;
NSD   = 2;

%% ===================== load raw session data =====================
infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
if ~isfile(infoPath), error('sessionFileInfo not found for %s %s', thisMouse, thisSessionName); end
loadedInfo      = load(infoPath, 'sessionFileInfo');
sessionFileInfo = loadedInfo.sessionFileInfo;
stimNames       = {sessionFileInfo.stimFiles.name};
iStim = find(contains(stimNames, stimName), 1);
if isempty(iStim), error('No %s file found for %s %s', stimName, thisMouse, thisSessionName); end

load(sessionFileInfo.stimFiles(iStim).Response, 'response');

if ~isfield(response, 'wheelData') || numel(response.wheelData) ~= numel(response.psthData)
    error('wheelData missing or mismatched with psthData for %s %s.', thisMouse, thisSessionName);
end

psthData = response.psthData;
stimVs   = vertcat(psthData.stimValue);

blankIdx  = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
gridMask  = stimVs(:,1) ~= 200;
gridPSTHIdx = find(gridMask);
gridStim  = stimVs(gridMask, :);

uAz      = sort(unique(gridStim(:,1)), 'ascend');
uEl_plot = sort(unique(gridStim(:,2)), 'descend');
nAz = numel(uAz); nEl = numel(uEl_plot);

timeVec = psthData(1).timeVector(:);
respIdx = timeVec >= respWin(1) & timeVec <= respWin(2);
baseIdx = timeVec >= baseWin(1) & timeVec < baseWin(2);
nTpts = numel(timeVec);

wheelTimeVec = response.wheelData(1).timeVector(:)';
wheelRespIdx = wheelTimeVec >= respWin(1) & wheelTimeVec <= respWin(2);

%% ===================== classify every trial (grid + blank) by state =====================
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

%% ===================== build boutonData struct for each condition, then plot =====================
condNames  = {'Stationary', 'Running', 'Combined (all trials)'};
condFilter = {0, 1, []}; % stateVal to filter to, or [] for no filter (combined)

figure('Color', 'w', 'Position', [50 50 1500 500]);
sgtitle(sprintf('%s | %s | ROI %d', thisMouse, thisSessionName, thisROI), 'Interpreter', 'none', 'FontWeight', 'bold');

for ci = 1:3
    stateVal = condFilter{ci};

    %% --- blank trials for this condition ---
    bTrials = squeeze(psthData(blankIdx).alignedResponses(thisROI, :, :));
    if isvector(bTrials), bTrials = bTrials(:); end
    if ~isempty(stateVal)
        blankMask = (runFlagByGroup{blankIdx} == stateVal);
        bTrials = bTrials(:, blankMask);
    end
    trialBaselinesB = mean(bTrials(baseIdx, :), 1, 'omitnan');
    bTrialsCorrected = bTrials - trialBaselinesB; % nTpts x nTrials
    meanBlankResponse = mean(bTrialsCorrected, 2, 'omitnan'); % nTpts x 1
    blankTrialMeans = mean(bTrialsCorrected(respIdx, :), 1, 'omitnan')';

    %% --- grid trials for this condition ---
    meanGridResponse = nan(nEl, nAz);
    meanTemporalResponse = nan(nTpts, nEl, nAz);
    maxResponseScale = 1e-6;

    allTrialMeans = blankTrialMeans(:);
    groupLabels   = repmat(numel(gridPSTHIdx) + 1, numel(blankTrialMeans), 1);

    for pIdx = 1:numel(gridPSTHIdx)
        g = gridPSTHIdx(pIdx);
        trialsAtPos = squeeze(psthData(g).alignedResponses(thisROI, :, :));
        if isvector(trialsAtPos), trialsAtPos = trialsAtPos(:); end
        if ~isempty(stateVal)
            posMask = (runFlagByGroup{g} == stateVal);
            trialsAtPos = trialsAtPos(:, posMask);
        end
        if isempty(trialsAtPos), continue; end

        trialBaselines = mean(trialsAtPos(baseIdx, :), 1, 'omitnan');
        correctedTrials = trialsAtPos - trialBaselines; % nTpts x nTrials

        avgAtPos = mean(correctedTrials, 2, 'omitnan'); % nTpts x 1
        maxResponseScale = max(maxResponseScale, max(avgAtPos));

        rowIdx = find(uEl_plot == gridStim(pIdx, 2), 1);
        colIdx = find(uAz == gridStim(pIdx, 1), 1);

        meanTemporalResponse(:, rowIdx, colIdx) = avgAtPos;
        meanGridResponse(rowIdx, colIdx) = mean(avgAtPos(respIdx), 'omitnan');

        posTrialMeans = mean(correctedTrials(respIdx, :), 1, 'omitnan')';
        allTrialMeans = [allTrialMeans; posTrialMeans(:)]; %#ok<AGROW>
        groupLabels   = [groupLabels; repmat(pIdx, numel(posTrialMeans), 1)]; %#ok<AGROW>
    end

    %% --- responsiveness (same logic as analyseRFMapping.m / RFMapping_BehaviorSplit.m) ---
    validIdx = ~isnan(allTrialMeans) & ~isnan(groupLabels);
    isResponsive = false; centerAz = NaN; centerEl = NaN;
    if sum(validIdx) > 0 && numel(unique(groupLabels(validIdx))) >= 2
        pValANOVA = anova1(allTrialMeans(validIdx), groupLabels(validIdx), 'off');
        blankMean = mean(blankTrialMeans, 'omitnan');
        blankStd  = std(blankTrialMeans, 'omitnan');
        [prefVal, mI] = max(meanGridResponse(:), [], 'omitnan');
        isResponsive = (pValANOVA < ALPHA) && (prefVal > (blankMean + NSD * blankStd));
        if isResponsive && ~isnan(mI)
            [rPeak, cPeak] = ind2sub(size(meanGridResponse), mI);
            centerAz = uAz(cPeak); centerEl = uEl_plot(rPeak);
        end
    end

    %% --- assemble boutonData struct matching plotRFHeatmapWithTraces's expected fields ---
    boutonData = struct( ...
        'meanGridResponse', meanGridResponse, ...
        'meanTemporalResponse', meanTemporalResponse, ...
        'meanBlankResponse', meanBlankResponse, ...
        'peakAmplitude', maxResponseScale, ...
        'centerAz', centerAz, ...
        'centerEl', centerEl, ...
        'isResponsive', isResponsive);

    %% --- plot using the existing function ---
    axPanel = subplot(1, 3, ci);
    plotRFHeatmapWithTraces(axPanel, boutonData, uAz, uEl_plot, timeVec, 'Smooth', true);
    title(axPanel, sprintf('%s (isResponsive=%d)', condNames{ci}, isResponsive), ...
        'FontSize', 10, 'FontWeight', 'bold');
end

end
