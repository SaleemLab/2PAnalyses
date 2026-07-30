% RFMapping_RunningStationaryScatterPerBouton.m
%
% For each pooled bouton, computes its mean evoked response separately on
% RUNNING trials and STATIONARY trials (ambiguous trials discarded), at
% that bouton's own preferred grid position, then scatters
% running-response vs stationary-response -- one point per bouton -- to
% ask whether averaging over locomotor state (as the main pipeline
% currently does) is masking a real state-dependent shift in visual
% responses.
%
% Running/stationary/ambiguous classification uses the same wheel-speed
% window & thresholds as RFMapping_OverallRunningStationaryFraction.m --
% this is deliberate, so the two scripts are self-consistent.
%
% REQUIRES (already in the workspace from the main pooling script +
% analyseRFBoutons_basic.m):
%   allRFMapping        - pooled bouton struct array
%   sessionLabels        - cell array, one label per bouton, '{mouse}_{session}'
%   RFMappingMetadata     - shared reference grid/time vector
%   uAz, uEl_plot, timeVector, respIdx  - already computed in the main script


%% params 
runSpeedThresh       = 3;        
statSpeedThresh      = 0.5;      
propThresh           = 0.75;     
stimFramesMask_range = [-0.2 2.8];  
                                  
                                  
minTrialsPerState    = 3;        
if ~exist('respIdx', 'var')
    error('respIdx (bouton response-window mask on timeVector) not found -- run the main pooling script first.');
end
numBoutons  = numel(allRFMapping);
runResp     = nan(numBoutons, 1);   
statResp    = nan(numBoutons, 1);   
nRunTrials  = nan(numBoutons, 1);
nStatTrials = nan(numBoutons, 1);
uniqueSessions = unique(sessionLabels, 'stable');
fprintf('Computing per-bouton running/stationary responses across %d sessions...\n', numel(uniqueSessions));
for s = 1:numel(uniqueSessions)
    thisLabel = uniqueSessions{s};
    boutonIdxThisSess = find(strcmp(sessionLabels, thisLabel));
    if isempty(boutonIdxThisSess)
        continue;
    end
    us              = strsplit(thisLabel, '_');
    thisMouse       = us{1};
    thisSessionName = strjoin(us(2:end), '_'); 
    infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
    if ~isfile(infoPath)
        warning('sfi missing for %s -- skipping.', thisLabel);
        continue;
    end
    loadedInfo      = load(infoPath, 'sessionFileInfo');
    sessionFileInfo = loadedInfo.sessionFileInfo;
    stimNamesHere   = {sessionFileInfo.stimFiles.name};
    iStim           = find(contains(stimNamesHere, 'RFMapping'), 1);
    if isempty(iStim)
        warning('No RFMapping file for %s -- skipping.', thisLabel);
        continue;
    end
    try
        load(sessionFileInfo.stimFiles(iStim).Response, 'response');
    catch ME
        warning('Could not load response for %s: %s', thisLabel, ME.message);
        continue;
    end
    if ~isfield(response, 'wheelData') || numel(response.wheelData) ~= numel(response.psthData)
        warning('wheelData missing/mismatched for %s -- skipping.', thisLabel);
        continue;
    end
    psthDataHere = response.psthData;
    stimVsHere   = vertcat(psthDataHere.stimValue);
    gridIdxHere  = find(stimVsHere(:,1) ~= 200);   
    wheelTimeVecHere = response.wheelData(1).timeVector(:)';
    wheelRespIdxHere = wheelTimeVecHere >= stimFramesMask_range(1) & wheelTimeVecHere <= stimFramesMask_range(2);
    trialState = cell(size(response.wheelData)); 
    for g = 1:numel(response.wheelData)
        wheelTrials = response.wheelData(g).alignedResponses;
        nTrialsHere = size(wheelTrials, 2);
        stateVec    = nan(nTrialsHere, 1);
        for ti = 1:nTrialsHere
            trace = wheelTrials(wheelRespIdxHere, ti);
            if all(isnan(trace)), continue; end
            meanSpeed      = nanmean(trace);
            propRunning    = sum(trace > statSpeedThresh) / sum(wheelRespIdxHere);
            propStationary = sum(trace < runSpeedThresh)  / sum(wheelRespIdxHere);
            if propRunning >= propThresh && meanSpeed > runSpeedThresh
                stateVec(ti) = 1;
            elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                stateVec(ti) = 0;
            end
        end
        trialState{g} = stateVec;
    end
    for k = 1:numel(boutonIdxThisSess)
        iROI = boutonIdxThisSess(k);
        b    = allRFMapping(iROI);
        
        trialMatrix = b.baselineSubtracted; 
        meanGridResponse = b.meanGridResponse;
        if isempty(trialMatrix) || isempty(meanGridResponse)
            continue;
        end
        
        [~, maxIdx] = max(double(meanGridResponse(:)));
        [rPeak, cPeak] = ind2sub(size(meanGridResponse), maxIdx);
        
        correctedTrials = trialMatrix{rPeak, cPeak};
        if isempty(correctedTrials)
            continue;
        end
        correctedTrials = double(correctedTrials);
        posTrialMeans   = mean(correctedTrials(:, respIdx), 2, 'omitnan'); 
        
        targetAz = uAz(cPeak);
        targetEl = uEl_plot(rPeak);
        
        gThis = [];
        for wG = 1:numel(response.wheelData)
            wVal = response.wheelData(wG).stimValue;
            if numel(wVal) >= 2 && wVal(1) == targetAz && wVal(2) == targetEl
                gThis = wG;
                break;
            end
        end
        if isempty(gThis) || gThis > numel(trialState)
            continue;
        end
        stateVec = trialState{gThis};
        if numel(stateVec) ~= numel(posTrialMeans)
            continue;
        end
        runMask  = stateVec == 1;
        statMask = stateVec == 0;
        if sum(runMask) >= minTrialsPerState && sum(statMask) >= minTrialsPerState
            runResp(iROI)     = mean(posTrialMeans(runMask), 'omitnan');
            statResp(iROI)    = mean(posTrialMeans(statMask), 'omitnan');
            nRunTrials(iROI)  = sum(runMask);
            nStatTrials(iROI) = sum(statMask);
        end
    end
end
validBoutons = ~isnan(runResp) & ~isnan(statResp);
fprintf('\n%d / %d boutons have >= %d running AND >= %d stationary trials.\n', ...
    sum(validBoutons), numBoutons, minTrialsPerState, minTrialsPerState);
statVals    = statResp(validBoutons);
runVals     = runResp(validBoutons);
isRespValid = logical([allRFMapping(validBoutons).isResponsive]);

statResp_respOnly = statVals(isRespValid);
runResp_respOnly  = runVals(isRespValid);
nAbove = sum(runResp_respOnly > statResp_respOnly);
nBelow = sum(runResp_respOnly < statResp_respOnly);
nTotalRespValid = numel(runResp_respOnly);

pctAbove = (nAbove / nTotalRespValid) * 100;
pctBelow = (nBelow / nTotalRespValid) * 100;

maxVal = max([statVals; runVals], [], 'omitnan');
minVal = min([statVals; runVals], [], 'omitnan');
padding = (maxVal - minVal) * 0.05;
axMin = minVal - padding;
axMax = maxVal + padding;   % FIX: was "maxVal - padding", which clipped the topmost points

figScatter = figure('Color', 'w', 'Position', [100 100 520 520], ...
    'Name', 'Running vs stationary response per bouton');
hold on;
plot([axMin axMax], [axMin axMax], 'k--', 'LineWidth', 1);
scatter(statVals(~isRespValid), runVals(~isRespValid), 25, [0.7 0.7 0.7], 'filled', ...
    'MarkerFaceAlpha', 0.4, 'DisplayName', sprintf('Non-responsive (\\geq%d trials)', minTrialsPerState));
hResp = scatter(statVals(isRespValid),  runVals(isRespValid),  35, [0.8 0.1 0.1], 'filled', ...
    'MarkerFaceAlpha', 0.8, 'DisplayName', sprintf('Responsive (\\geq%d trials)', minTrialsPerState));
if ~isempty(hResp)
    hResp.MarkerEdgeColor = [0.5 0.0 0.0];
    hResp.LineWidth = 0.5;
end
text(axMin + (axMax-axMin)*0.05, axMax - (axMax-axMin)*0.05, ...
    {sprintf('Enhanced (Above): %d/%d (%.1f%%)', nAbove, nTotalRespValid, pctAbove), ...
     sprintf('Suppressed (Below): %d/%d (%.1f%%)', nBelow, nTotalRespValid, pctBelow)}, ...
    'FontName', 'Arial', 'FontSize', 9, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'top');
xlabel('Stationary-trial response (\DeltaF/F)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('Running-trial response (\DeltaF/F)', 'FontName', 'Arial', 'FontSize', 10);
legend('Location', 'southeast', 'Box', 'off', 'FontName', 'Arial', 'FontSize', 9);
xlim([axMin axMax]); ylim([axMin axMax]);
axis square;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
title({sprintf('Locomotor State Shift at Preferred Grid Position'), ...
       sprintf('(n = %d total filtered boutons)', sum(validBoutons))}, ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');

%% ===================== signed-rank tests: full population AND responsive-only =====================
% NOTE: these two populations are NOT the same thing and should be reported
% separately -- the full-filtered population (n ~= all validBoutons) mixes
% responsive and non-responsive boutons together, which can wash out an
% effect that is real and strong within the responsive (visually-tuned)
% subset. The 63/103-style enhanced/suppressed count above was already
% computed on the responsive-only subset, so the magnitude test below is
% matched to it.

if sum(validBoutons) >= 5
    [p_signrank_all, ~, statsOut_all] = signrank(runVals, statVals);
    medianDiff_all = median(runVals - statVals, 'omitnan');
    fprintf('\n[FULL population] Wilcoxon signed-rank (running vs stationary), n = %d: p = %.4f | median(run - stat) = %.4f\n', ...
        sum(validBoutons), p_signrank_all, medianDiff_all);
else
    fprintf('\n[FULL population] Fewer than 5 boutons pass the trial-count criteria -- skipping signed-rank test.\n');
end

if nTotalRespValid >= 5
    [p_signrank_resp, ~, statsOut_resp] = signrank(runResp_respOnly, statResp_respOnly);
    medianDiff_resp = median(runResp_respOnly - statResp_respOnly, 'omitnan');
    fprintf('[RESPONSIVE-ONLY] Wilcoxon signed-rank (running vs stationary), n = %d: p = %.4f | median(run - stat) = %.4f\n', ...
        nTotalRespValid, p_signrank_resp, medianDiff_resp);

    % direction-count test (independent of magnitude): are enhanced/suppressed
    % counts themselves different from chance? complements the signrank test,
    % which tests magnitude+direction jointly.
    p_binom_resp = 2 * min(binocdf(nAbove, nTotalRespValid, 0.5), 1 - binocdf(nAbove - 1, nTotalRespValid, 0.5));
    fprintf('[RESPONSIVE-ONLY] Binomial test (enhanced vs suppressed counts, n = %d): %d above / %d below, p = %.4f\n', ...
        nTotalRespValid, nAbove, nBelow, p_binom_resp);
else
    fprintf('[RESPONSIVE-ONLY] Fewer than 5 responsive boutons pass the trial-count criteria -- skipping tests.\n');
end

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\running_vs_stationary_scatter';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figScatter, 'Visible', 'off');
saveFigureFormats(figScatter, fullfile(outputDir, 'running_vs_stationary_response_perBouton'));