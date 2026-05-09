function [sessionFileInfo, RFMapping, RFMappingMetadata, allCenters] = analyseRFMapping(sessionFileInfo, stimName)
% analyseRFMapping: Classifies visual responsiveness based on ANOVA and Signal-to-Noise.
% This version uses the MEAN response with a relaxed 1.5 SD threshold for debugging.

% 1. Load response data
iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
load(sessionFileInfo.stimFiles(iStim).Response, 'response');
psthData = response.psthData;
stimVs = vertcat(psthData.stimValue);
nROI = size(psthData(1).alignedResponses, 1);

% Filter grid and blank response [200 0]
blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
gridMask = stimVs(:,1) ~= 200;
gridPSTH = psthData(gridMask);
gridStim = stimVs(gridMask, :);
blankPSTH = psthData(blankIdx);

% Grid Arrangement
uAz = sort(unique(gridStim(:,1)), 'ascend');  
uEl_plot = sort(unique(gridStim(:,2)), 'descend'); 
nAz = length(uAz); nEl = length(uEl_plot);
timeVector = psthData(1).timeVector(:);
nTpts = length(timeVector);

% --- FIXED PARAMETERS ---
respWin = [0.5 3]; 
baseWin = [-1.0 0];  
respIdx = timeVector >= respWin(1) & timeVector <= respWin(2);
baseIdx = timeVector >= baseWin(1) & timeVector < baseWin(2);

% Initialize structure 
RFMapping = struct('meanGridResponse', cell(nROI, 1), ... 
                   'meanTemporalResponse', cell(nROI, 1), ...
                   'meanBlankResponse', cell(nROI, 1), ...
                   'baselineSubtracted', cell(nROI, 1), ...      
                   'baselineSubtractedBlank', cell(nROI, 1), ... 
                   'peakAmplitude', cell(nROI, 1), ...
                   'centerAz', cell(nROI, 1), ...
                   'centerEl', cell(nROI, 1), ...
                   'pValANOVA', cell(nROI, 1), ...
                   'isResponsive', cell(nROI, 1)); 
allCenters = nan(nROI, 2);

%% Loop for each ROI 
for iROI = 1:nROI
    meanGridResponse = nan(nEl, nAz); 
    peakGridResponse = nan(nEl, nAz); 
    temporalStack = nan(nTpts, nEl, nAz);
    trialMatrix = cell(nEl, nAz); 
    maxResponseScale = 1e-6;
    
    % Accumulators for ANOVA [Observations x 1]
    allTrialMeans = []; 
    groupLabels = [];
    
    % --- Process Blank trials ---
    bTrials = blankPSTH.alignedResponses(iROI, :, :);
    trialBaselinesB = mean(bTrials(1, baseIdx, :), 2, 'omitnan'); 
    bTrialsCorrected = bTrials - trialBaselinesB; 
    meanBlankResponse = mean(bTrialsCorrected, 3, 'omitnan');
    
    % Collect Blank trial-means for ANOVA and Noise Floor
    blankTrialMeans = squeeze(mean(bTrialsCorrected(1, respIdx, :), 2, 'omitnan'));
    allTrialMeans = [allTrialMeans; blankTrialMeans(:)];
    groupLabels = [groupLabels; repmat(17, numel(blankTrialMeans), 1)]; % Group 17 is Blank 
    
    % --- Process Grid trials ---
    for thisPos = 1:numel(gridPSTH)
        allTrialsAtPos = gridPSTH(thisPos).alignedResponses(iROI, :, :);
        trialBaselines = mean(allTrialsAtPos(1, baseIdx, :), 2, 'omitnan');
        correctedTrials = allTrialsAtPos - trialBaselines; 
        
        avgAtPos = mean(correctedTrials, 3, 'omitnan');
        avgAtPos = avgAtPos(:);
        maxResponseScale = max(maxResponseScale, max(avgAtPos));
        
        rowIdx = find(uEl_plot == gridStim(thisPos, 2), 1); 
        colIdx = find(uAz == gridStim(thisPos, 1), 1);
        
        temporalStack(:, rowIdx, colIdx) = avgAtPos;
        
        % Store Mean (used for primary gate in this version)
        meanGridResponse(rowIdx, colIdx) = mean(avgAtPos(respIdx), 'omitnan');
        % Store Peak (used for coordinate centering)
        peakGridResponse(rowIdx, colIdx) = max(avgAtPos(respIdx), [], 'omitnan');
        
        % Collect trial-means for ANOVA 
        posTrialMeans = squeeze(mean(correctedTrials(1, respIdx, :), 2, 'omitnan'));
        allTrialMeans = [allTrialMeans; posTrialMeans(:)];
        groupLabels = [groupLabels; repmat(thisPos, numel(posTrialMeans), 1)];
        
        trialMatrix{rowIdx, colIdx} = squeeze(correctedTrials)'; 
    end
    
    % Statistical tests 
    
    %  ANOVA across all locations + blank (p < 0.05) 
    pValANOVA = anova1(allTrialMeans, groupLabels, 'off');
    % Calculate blank trial statistics
    blankStd  = std(blankTrialMeans, 'omitnan');
    blankMean = mean(blankTrialMeans, 'omitnan');
    
    % Find the index of the preferred location (the bin with the highest mean response)
    [~, mI] = max(meanGridResponse(:));
    
    % Use that index to get the mean response at that preferred location
    prefVal = meanGridResponse(mI);
    
    % Logic: ANOVA p < 0.05 AND Mean(Pref) > (Mean Blank + 2*SD Blank)
    isResponsive = (pValANOVA < 0.05) && (prefVal > (blankMean + 2 * blankStd));
    
    % Final Storage
    RFMapping(iROI).meanGridResponse = meanGridResponse; 
    RFMapping(iROI).meanTemporalResponse = temporalStack;
    RFMapping(iROI).meanBlankResponse = meanBlankResponse(:);
    RFMapping(iROI).baselineSubtracted = trialMatrix;
    RFMapping(iROI).baselineSubtractedBlank = squeeze(bTrialsCorrected)'; 
    RFMapping(iROI).peakAmplitude = maxResponseScale;
    RFMapping(iROI).pValANOVA = pValANOVA;
    RFMapping(iROI).isResponsive = isResponsive;
    
    if isResponsive
        % Fix: Use linear index to avoid empty find results from floating point errors
        [~, mI] = max(meanGridResponse(:));
        [rPeak, cPeak] = ind2sub(size(meanGridResponse), mI);
        
        RFMapping(iROI).centerAz = uAz(cPeak);
        RFMapping(iROI).centerEl = uEl_plot(rPeak);
        allCenters(iROI, :) = [uAz(cPeak), uEl_plot(rPeak)];
    else
        allCenters(iROI, :) = [NaN, NaN];
    end
end

% Save results - Ensuring uAz and uEl are included for plotRFMapping
RFMappingMetadata = struct('uAz', uAz, 'uEl', uEl_plot, 'timeVector', timeVector, ...
    'respWin', respWin, 'baseWin', baseWin, 'zThreshold', 1);

savePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
save(savePath, 'RFMapping', 'RFMappingMetadata', '-append');

fprintf('Processing complete. %d/%d ROIs classified as visually responsive.\n', sum([RFMapping.isResponsive]), nROI);
end