function [sessionFileInfo, RFMapping, RFMappingMetadata, allCenters] = analyseRFMapping(sessionFileInfo, stimName)
% Uses the psth data from response, seperates grid trials and blank trials;
% Identifies unique Az and El position; 
% For every position grid each trial is baseline subtracted before mean for
% a given position is saved. 
% load data
iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
load(sessionFileInfo.stimFiles(iStim).Response, 'response');
psthData = response.psthData;
stimVs = vertcat(psthData.stimValue);
nROI = size(psthData(1).alignedResponses, 1);
% filter grid and blank response [200 0] are the blank trials
blankIdx = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
gridMask = stimVs(:,1) ~= 200;
gridPSTH = psthData(gridMask);
gridStim = stimVs(gridMask, :);
blankPSTH = psthData(blankIdx);
% arrangement of grids
uAz = sort(unique(gridStim(:,1)), 'ascend');  
uEl_plot = sort(unique(gridStim(:,2)), 'descend'); 
nAz = length(uAz); nEl = length(uEl_plot);
timeVector = psthData(1).timeVector(:);
nTpts = length(timeVector);
%this is currently hardcoded 
respWin = [0.5 3];  % changed from 2 t0 3: 24/03
baseWin = [-0.75 0];  %changed 24/03 
% create structure 
RFMapping = struct('meanGridResponse', cell(nROI, 1), ... 
                   'meanTemporalResponse', cell(nROI, 1), ...
                   'meanBlankResponse', cell(nROI, 1), ...
                   'correctedTrialsGrid', cell(nROI, 1), ... 
                   'correctedTrialsBlank', cell(nROI, 1), ...
                   'peakAmplitude', cell(nROI, 1), ...
                   'centerAz', cell(nROI, 1), ...
                   'centerEl', cell(nROI, 1));
allCenters = nan(nROI, 2);
%% Loop for each roi 
for iROI = 1:nROI
    meanGridResponse = nan(nEl, nAz); 
    temporalStack = nan(nTpts, nEl, nAz); 
    correctedGridTrials = cell(nEl, nAz); 
    maxResponseScale = 1e-6; % Gemini's suggestion to prevent normalisation from crashing.. 
    
    %  Process Blank trials [trial-by-trial baseline]
    % bTrials shape: [1, nTpts, nTrials]
    bTrials = blankPSTH.alignedResponses(iROI, :, :);
    baseIdxB = timeVector >= baseWin(1) & timeVector < baseWin(2);
    
    % Calculate baseline for EACH trial and subtract it from that trial
    trialBaselinesB = mean(bTrials(1, baseIdxB, :), 2, 'omitnan'); % Mean of baseline window per trial
    bTrialsCorrected = bTrials - trialBaselinesB; % MATLAB implicitly expands trialBaselinesB
    
    % Now average the corrected trials
    meanBlankResponse = mean(bTrialsCorrected, 3, 'omitnan');
    meanBlankResponse = meanBlankResponse(:); 
    % Process Grid trials
    for thisPos = 1:numel(gridPSTH)
        % allTrialsAtPos shape: [1, nTpts, nTrials]
        allTrialsAtPos = gridPSTH(thisPos).alignedResponses(iROI, :, :);
        
        % Trial-by-trial baseline subtraction
        baseIdx = timeVector >= baseWin(1) & timeVector < baseWin(2);
        trialBaselines = mean(allTrialsAtPos(1, baseIdx, :), 2, 'omitnan');
        correctedTrials = allTrialsAtPos - trialBaselines; % Subtract each trial's own baseline
        
        % Now average the corrected trials to get the mean PSTH for this position
        avgAtPos = mean(correctedTrials, 3, 'omitnan');
        avgAtPos = avgAtPos(:); 
        
        % track peak for scaling
        maxResponseScale = max(maxResponseScale, max(avgAtPos));
        
        % map to indices
        rowIdx = find(uEl_plot == gridStim(thisPos, 2), 1); 
        colIdx = find(uAz == gridStim(thisPos, 1), 1);
        
        % time series and grid mean
        temporalStack(:, rowIdx, colIdx) = avgAtPos;
        meanGridResponse(rowIdx, colIdx) = mean(avgAtPos(timeVector >= respWin(1) & timeVector <= respWin(2)), 'omitnan');
        
        % Save individual baseline-subtracted trials
        correctedGridTrials{rowIdx, colIdx} = squeeze(correctedTrials);
    end
    
    % average response magnitude at each (Elevation, Azimuth) location.
    RFMapping(iROI).meanGridResponse = meanGridResponse;
    % Full trial-averaged time-courses for every grid position.
    RFMapping(iROI).meanTemporalResponse = temporalStack;
    % Average activity during blank trials.
    RFMapping(iROI).meanBlankResponse = meanBlankResponse;
    
    % Store trial-by-trial data
    RFMapping(iROI).correctedTrialsGrid = correctedGridTrials;
    RFMapping(iROI).correctedTrialsBlank = squeeze(bTrialsCorrected);
    
    % The highest average response recorded across all positions
    RFMapping(iROI).peakAmplitude = maxResponseScale;
    
    % find peaks in all grids
    [~, mI] = max(meanGridResponse(:));
    [rPeak, cPeak] = ind2sub(size(meanGridResponse), mI);
    RFMapping(iROI).centerAz = uAz(cPeak);
    RFMapping(iROI).centerEl = uEl_plot(rPeak);
    allCenters(iROI, :) = [uAz(cPeak), uEl_plot(rPeak)];
end
% save the metadata 
RFMappingMetadata = struct('stimName', stimName, 'timeVector', timeVector, ...
    'uAz', uAz, 'uEl', uEl_plot, 'respWin', respWin, 'baseWin', baseWin);
savePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
save(savePath, 'RFMapping', 'RFMappingMetadata', '-append');
fprintf('RF Mapping for %d ROIs saved to: %s\n', nROI, savePath);
end