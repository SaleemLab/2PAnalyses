function crossValExpVar = getCrossValidatedExplainedVariance(sessionFileInfo, VRStimName, response, nFolds, excludeFlaggedLaps, applySmoothing, plotFlag, isOpenLoop)
% Calculates n-fold cross-validated explained variance of mean tuning curves for each ROI,
% generates a shuffled null distribution to test significance,
% and appends the results to sessionROIData.
% TODO: change to run across all signals; or at least dff neu and spikes 
%
% Input:
%   sessionFileInfo     - Structure containing directory/file paths
%   response            - Structure containing activity data tensors
%   signalToUse         - String key for the signal data type (Optional/Deprecated - now loops through all)
%   nFolds              - Number of validation folds (default: 5)
%   excludeFlaggedLaps  - Boolean to filter out manually flagged artifact laps
%   plotFlag            - Boolean to display summary figures
%
% Output:
%   crossValExpVar      - Structure containing true metrics, null distribution, and p-values
%% Handle  inputs
% Load variables if response is not provided 
if nargin < 3 || isempty(response)
    stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
    if isempty(stimIdx)
        error('Specify valid vr stim name or the response struct \n.');
    end
    response = load(sessionFileInfo.stimFiles(stimIdx).Response, 'lapPositionActivity', 'trialIndicesByCondition', 'stimName', 'flaggedLaps');
end 
% Define other defualt variables
if nargin < 4 || isempty(nFolds); nFolds = 5; end
if nargin < 5 || isempty(excludeFlaggedLaps); excludeFlaggedLaps = true; end
if nargin < 5 || isempty(applySmoothing); applySmoothing = true; end
if nargin < 7 || isempty(plotFlag); plotFlag = true; end
if nargin <8  || isempty(isOpenLoop); isOpenLoop = false; end 
% Define number of shuffles for the null distribution
nShuffles = 1000; 
% Threshold for explained variance; this is used for plotting; threshold can be applied later on  
%threshold = 0.1;
    
% Identify all signals available in the response file
signalNames = fieldnames(response.lapPositionActivity);
disp('Extracting and processing cross-validated explained variance for all signals...');

% Initialize the main output structure
crossValExpVar = struct();

%% Loop through every signal 
for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};
    fprintf('\n=== Processing signal: %s ===\n', currentSignalName);

    %% Get and filter baseline data
    % Extract the raw activity matrix [ROI x Laps x Position]
    % signals =  fieldnames(response.lapPositionActivity);
    rawActivity = response.lapPositionActivity.(currentSignalName);
    nTotalLaps = size(rawActivity, 2);
    unflaggedMask = true(1, nTotalLaps);
    baseTrialsMask = false(1, nTotalLaps);
    
    % Identify manually flagged artifact laps if requested
    if excludeFlaggedLaps && isfield(response, 'flaggedLaps')
        disp('Excluding flagged laps before computing this metric...');
        unflaggedMask(response.flaggedLaps) = false;
    end 
    % Identify baseline trials
    if isfield(response, 'trialIndicesByCondition') && isfield(response.trialIndicesByCondition, 'Baseline')
        disp('Isolating Baseline condition trials...');
        baseTrialsMask(response.trialIndicesByCondition.Baseline) = true;
    else
        error('Field "response.trialIndicesByCondition.Baseline" is required to find baseline trials.');
    end
    % Intersect masks to get final clean, baseline trials
    finalValidLapsMask = unflaggedMask & baseTrialsMask;
    lapPositionActivity = rawActivity(:, finalValidLapsMask, :);
    % Establish size metrics based on filtered subset
    [nROIs, numLaps, numPositions] = size(lapPositionActivity);

    % apply smoothning
    if applySmoothing
        w = gausswin(15); w = w / sum(w);
        for iCell = 1:nROIs
            for iLap = 1:numLaps
                trace = squeeze(lapPositionActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace);
                trace(nanMask) = 0;
                smoothed = filtfilt(w, 1, trace);
                smoothed(nanMask) = NaN;
                lapPositionActivity(iCell, iLap, :) = smoothed;
            end
        end
    end

    fprintf('Analyzing %d baseline/clean laps out of %d total recorded laps.\n', numLaps, nTotalLaps);
    
    %% Compute cross validated variance explained 
    fprintf('Computing cross-validated explained variance across %d folds...\n', nFolds);
    cvExpVar = NaN(nROIs, nFolds);
    % null distribution matrix: [nROIs x nShuffles]
    cvExpVarNull = NaN(nROIs, nShuffles);
    for iCell = 1:nROIs
        cellActivity = squeeze(lapPositionActivity(iCell, :, :)); % [Laps x Position]
        
        if ~all(isnan(cellActivity), 'all')
            % Compute the TRUE cross-validated explained variance
            cvExpVar(iCell, :) = getAllExpVar(cellActivity, nFolds); 
            
            % Compute the SHUFFLED null distribution for this ROI
            for thisShuff = 1:nShuffles
                shuffledActivity = NaN(numLaps, numPositions);
                s = RandStream('mt19937ar','Seed',thisShuff);
                
                % Circularly shift position columns independently for each lap
                for iLap = 1:numLaps
                    lapData = cellActivity(iLap, :);
                    if ~all(isnan(lapData))
                        % Shift by a random amount between 1 and the number of position bins
                        randomShift = randi(s, numPositions);
                        shuffledActivity(iLap, :) = circshift(lapData, randomShift);
                    else
                        shuffledActivity(iLap, :) = lapData; % keep NaNs as-is
                    end
                end
                
                % Run cross-validation on the scrambled data
                shuffledFolds = getAllExpVar(shuffledActivity, nFolds);
                
                % Store the fold average for this shuffle item
                cvExpVarNull(iCell, thisShuff) = mean(shuffledFolds, 'omitnan');
            end
        end
    end
    % Calculate performance metrics across folds
    medianExpVar = median(cvExpVar, 2, 'omitnan');
    meanExpVar   = mean(cvExpVar, 2, 'omitnan');
    % Count how many shuffles equaled or beat the true mean
    countBetterNulls = sum(cvExpVarNull >= meanExpVar, 2, 'omitnan');
    % Instead of assuming the data fits a standard textbook curve (like a normal distribution or a t-test),
    %  it counts how the real biological data compares directly to the 1,000 randomised shuffles.
    pValues = (countBetterNulls + 1) / (nShuffles + 1); % countBetterNulls + 1 (+1 prevents it from being zero; nshuffled +1 This is the total size comparison pool (the 1,000 shuffles + 1 real dataset = 1,001 total possibilities).
    % Calculate threshold percentiles from the null 
    null95thPercentile = prctile(cvExpVarNull, 95, 2); % alpha = 0.05 
    null99thPercentile = prctile(cvExpVarNull, 99, 2); % alpha = 0.01 
    
    %% Save results to nested signal structure
    crossValExpVar.(currentSignalName).cvExpVar = cvExpVar;
    crossValExpVar.(currentSignalName).medianExpVar = medianExpVar;
    crossValExpVar.(currentSignalName).meanExpVar = meanExpVar;
    crossValExpVar.(currentSignalName).cvExpVarNull = cvExpVarNull; 
    crossValExpVar.(currentSignalName).null95thPercentile = null95thPercentile;
    crossValExpVar.(currentSignalName).null99thPercentile = null99thPercentile;
    crossValExpVar.(currentSignalName).pValues = pValues;
    crossValExpVar.(currentSignalName).nFolds = nFolds;
    crossValExpVar.(currentSignalName).nShuffles = nShuffles;

    if isOpenLoop

        OpenLoop.crossValExpVar.(currentSignalName).cvExpVar = cvExpVar;
        OpenLoop.crossValExpVar.(currentSignalName).medianExpVar = medianExpVar;
        OpenLoop.crossValExpVar.(currentSignalName).meanExpVar = meanExpVar;
        OpenLoop.crossValExpVar.(currentSignalName).cvExpVarNull = cvExpVarNull;
        OpenLoop.crossValExpVar.(currentSignalName).null95thPercentile = null95thPercentile;
        OpenLoop.crossValExpVar.(currentSignalName).null99thPercentile = null99thPercentile;
        OpenLoop.crossValExpVar.(currentSignalName).pValues = pValues;
        OpenLoop.crossValExpVar.(currentSignalName).nFolds = nFolds;
        OpenLoop.crossValExpVar.(currentSignalName).nShuffles = nShuffles;

    end
    % crossValExpVar.(currentSignalName).signalUsed = currentSignalName;

    %% Plot
    if plotFlag
        % A cell is significantly stable if its real variance is greater than 
        % its individual 99th percentile null threshold (p < 0.05)
        significantCellsIdx = find(pValues < 0.01); %& crossValExpVar.meanExpVar > threshold); 
        
        fprintf('[%s] Found %d statistically significant, stable cells out of %d total ROIs.\n', ...
            currentSignalName, length(significantCellsIdx), length(meanExpVar));
            
        % Overall distribution 
        figure('Name', ['Cross-Val Summary: ' currentSignalName], 'Color', 'w', 'Position', [100, 100, 1100, 400]);
        subplot(1, 4, 1); 
        
        histogram(meanExpVar, 'BinWidth', 0.01, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.6);
       
        xlabel('Explained Variance (R^2)');
        ylabel('Probability');
        title(sprintf('Population Distribution (%s)\n(%d Stable ROIs (p < 0.01)', currentSignalName, length(significantCellsIdx)), 'Interpreter', 'none');
        legend('Real Data (Mean)', 'Location', 'NorthEast');
        set(gca, 'TickDir', 'out', 'box', 'off');
        legend boxoff;
        
        % Example ROIs
        [~, sortedIdx] = sort(pValues);
        
        % First, middle and last
        if length(sortedIdx) >= 3
            exampleROIs = [sortedIdx(1), sortedIdx(round(length(sortedIdx)/2)), sortedIdx(end)];
            roiLabels = {'high', 'Middle', 'low'};
        else
            exampleROIs = 1:min(3, nROIs);
            roiLabels = repmat({'ROI'}, 1, 3);
        end
        
        for idx = 1:length(exampleROIs)
            roiIdx = exampleROIs(idx);
            subplot(1, 4, 1 + idx);
            
            histogram(cvExpVarNull(roiIdx, :), 'BinWidth', 0.01, ...
                'Normalization', 'probability', ...
                'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'w', 'FaceAlpha', 0.5);
            hold on;
            
            % Vertical line for the 99th percentile confidence limit
            xline(null95thPercentile(roiIdx), 'r--', 'LineWidth', 1.5);
            % Vertical line for the 99th percentile confidence limit
            xline(null99thPercentile(roiIdx), 'g--', 'LineWidth', 1.5);
            
            % Vertical line for the actual true mean performance
            xline(meanExpVar(roiIdx), 'b-', 'LineWidth', 2);
            
            xlabel('R^2');
            title(sprintf('%s (ROI %d)\np = %0.3f', roiLabels{idx}, roiIdx, pValues(roiIdx)), 'FontSize', 9);
            set(gca, 'TickDir', 'out', 'box', 'off');
            
            if idx == 1
                ylabel('Probability');
                legend({'Null Dist', '95% CI Limit', 'True Mean'}, 'Location', 'NorthWest', 'FontSize', 8);
                legend boxoff;
            end
        end
        drawnow;
    end
end

%% Save results
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    disp(['Saving Cross-Val ExpVar and Null results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, "crossValExpVar", '-append')
elseif isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2 && isOpenLoop
    disp(['Saving Cross-Val ExpVar and Null results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, "crossValExpVar", '-append')
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append cross-val data.', sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save cross-val data.');
end

disp('Cross-validation calculation complete for all signals.');
end