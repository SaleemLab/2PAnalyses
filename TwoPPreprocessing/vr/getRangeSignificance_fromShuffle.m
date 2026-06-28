function [isSignificantByRangeShuffling_local, realRangePercentileRank_local] = getRangeSignificance_fromShuffle(sessionFileInfo, response)
% Calculates the range (Max - Min) of each cell's spatial tuning curve,
% compares it to the distribution of ranges from the shuffled data,
% and determines which cells have a range greater than the 95th percentile.
%
% Aman and Sonali November 2025

%% Input Checks and Preparation
if nargin < 2
    error('Must provide the response struct.');
end

if ~isfield(response, 'lapPositionActivity_ShuffleMatrix')
    error('Shuffled data (lapPositionActivity_ShuffleMatrix) not found in response struct. Run getLapPositionActivity_withShuffle first.');
end

% Default signal list
signalNames = {'dFF', 'dFFNeuropilCorrected'};

% Dynamic integration of spks if available in both real and shuffled data fields
if isfield(response.lapPositionActivity, 'spks') && isfield(response.lapPositionActivity_ShuffleMatrix, 'spks')
    signalNames{end+1} = 'spks';
end

% Index of significant cells (1 = significant, 0 = not significant)
isSignificantByRangeShuffling_local = struct(); % Local container for saving
% Percentile of real peak relative to shuffle peaks
realRangePercentileRank_local = struct(); % Local container for saving
% Target percentile for significance
targetPercentile = 95;

%% Main Processing Loop
disp('Starting range significance calculation...');
for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};
    
    % Safety guard if a signal is unexpectedly missing from a container
    if ~isfield(response.lapPositionActivity, currentSignalName) || ...
       ~isfield(response.lapPositionActivity_ShuffleMatrix, currentSignalName)
        warning('Field %s not found across all tracking matrices. Skipping.', currentSignalName);
        continue;
    end
    
    disp(['Processing signal: ', currentSignalName]);
    
    % Real activity (Neurons x Laps x Bins)
    realActivity = response.lapPositionActivity.(currentSignalName); 
    % Shuffled activity (Neurons x Bins x Shuffles)
    shuffleMatrix = response.lapPositionActivity_ShuffleMatrix.(currentSignalName);
    
    numROIs = size(realActivity, 1);
    numShuffles = size(shuffleMatrix, 3);
    
    % Calculate the mean real tuning curve
    % Average activity across all laps to get the real tuning curve (Neurons x Bins)
    meanRealTuningCurve = mean(realActivity, 2, 'omitnan'); 
    meanRealTuningCurve = squeeze(meanRealTuningCurve); % (numROIs x numBins)
    
    % Handle single-ROI squeeze edge case (preserves horizontal bin structure)
    if numROIs == 1 && size(meanRealTuningCurve, 1) ~= 1
        meanRealTuningCurve = meanRealTuningCurve';
    end
    
    % Find real range
    realMaxValues = max(meanRealTuningCurve, [], 2, 'omitnan'); % (numROIs x 1)
    realMinValues = min(meanRealTuningCurve, [], 2, 'omitnan'); % (numROIs x 1)
    realRangeValues = realMaxValues - realMinValues; % (numROIs x 1)
    
    % Find shuffle range 
    shuffleMaxValues = squeeze(max(shuffleMatrix, [], 2, 'omitnan')); % (numROIs x numShuffles)
    shuffleMinValues = squeeze(min(shuffleMatrix, [], 2, 'omitnan')); % (numROIs x numShuffles)
    
    % Handle single-ROI squeeze edge case for shuffle extractions
    if numROIs == 1
        if size(shuffleMaxValues, 1) ~= 1, shuffleMaxValues = shuffleMaxValues'; end
        if size(shuffleMinValues, 1) ~= 1, shuffleMinValues = shuffleMinValues'; end
    end
    
    % Calculate the shuffle range distribution (Max - Min) for every shuffle
    shuffleRangeDistribution = shuffleMaxValues - shuffleMinValues; % (numROIs x numShuffles)
    
    % Sig and percentiles 
    significanceIndex = zeros(numROIs, 1);
    percentileValues = nan(numROIs, 1);
    
    % Loop over each ROI to compare its real range against the shuffle range distribution
    for thisROI = 1:numROIs
        % Distribution of range values for the current neuron across all shuffles
        currentShuffleDistribution = shuffleRangeDistribution(thisROI, :);
        % The actual range value for the current neuron
        currentRealRange = realRangeValues(thisROI);
        
        % Calculate the target significance threshold (95th percentile)
        rangeThreshold_95th = prctile(currentShuffleDistribution, targetPercentile);
        
        % Check for significance: Is the Real Range > 95th Percentile Threshold?
        if currentRealRange > rangeThreshold_95th
            significanceIndex(thisROI) = 1; % Mark as significant
        end
        
        % Calculate the percentile of the real range relative to the shuffle distribution
        numRangesBelowReal = sum(currentShuffleDistribution <= currentRealRange);
        percentile = (numRangesBelowReal / numShuffles) * 100;
        percentileValues(thisROI) = percentile;
    end
    
    isSignificantByRangeShuffling_local.(currentSignalName) = significanceIndex;
    realRangePercentileRank_local.(currentSignalName) = percentileValues;
    
    fprintf('  -> Found %d significant ROIs (Range > %dth Percentile) for %s.\n', ...
        sum(significanceIndex), targetPercentile, currentSignalName);
end
disp('Range significance calculation complete.');

%% Saving Data to sessionROIData 
nullDist_RangeTuningMetric.isSignificantByRange = isSignificantByRangeShuffling_local;
nullDist_RangeTuningMetric.realRangePercentileRank = realRangePercentileRank_local;
nullDist_RangeTuningMetric.targetPercentile = targetPercentile; 

% Check if the file path exists and save the variables
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    disp(['Saving range significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
        "nullDist_RangeTuningMetric", ...
         '-append')
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append range significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save range significance data.');
end
end
% %% Testing random permutations 
% function [nullDist_RangeTuningMetric] = getRangeSignificance_fromShuffle(sessionFileInfo, response)
% % Calculates the range (Max - Min) of each cell's spatial tuning curve,
% % compares it to the distribution of ranges from the shuffled data,
% % and determines which cells have a range greater than the 95th percentile.
% %
% % Updated to process both standard shuffle using circular shifts and RandPerm shuffle metrics.
% %
% % Aman and Sonali November 2025 (Updated May 2026)
% 
% %% Input Checks and Preparation
% if nargin < 2
%     error('Must provide the response struct.');
% end
% 
% signalNames = {'dFF', 'dFFNeuropilCorrected'};
% targetPercentile = 95;
% 
% % Identify which shuffle matrices exist in the response struct
% shuffleTypesToRun = {};
% if isfield(response, 'lapPositionActivity_ShuffleMatrix')
%     shuffleTypesToRun{end+1} = 'lapPositionActivity_ShuffleMatrix';
% end
% if isfield(response, 'lapPositionActivity_RandPermMatrix')
%     shuffleTypesToRun{end+1} = 'lapPositionActivity_RandPermMatrix';
% end
% 
% if isempty(shuffleTypesToRun)
%     error('No shuffled matrices (ShuffleMatrix or RandPermMatrix) found in response struct.');
% end
% 
% % Initialize the unified tracking structure
% nullDist_RangeTuningMetric = struct();
% nullDist_RangeTuningMetric.targetPercentile = targetPercentile;
% 
% %% Main Processing Loop
% disp('Starting range significance calculation...');
% 
% for iType = 1:length(shuffleTypesToRun)
%     currentShuffleField = shuffleTypesToRun{iType};
% 
%     % Create clean identifiers for storing results
% 
%     if strcmp(currentShuffleField, 'lapPositionActivity_ShuffleMatrix')
%         shufLabel = 'CircShift';
%     else
%         shufLabel = 'RandPerm';
%     end
% 
%     disp(['=== Running Analysis for Shuffle Type: ', shufLabel, ' ===']);
% 
%     for iSignal = 1:length(signalNames)
%         currentSignalName = signalNames{iSignal};
% 
%         % Check if this specific signal exists in the current shuffle cube
%         if ~isfield(response.(currentShuffleField), currentSignalName)
%             warning('Signal %s not found in %s. Skipping.', currentSignalName, currentShuffleField);
%             continue;
%         end
% 
%         disp(['Processing signal: ', currentSignalName]);
% 
%         % Real activity (Neurons x Laps x Bins)
%         realActivity = response.lapPositionActivity.(currentSignalName); 
% 
%         % Shuffled activity (Neurons x Bins x Shuffles)
%         shuffleMatrix = response.(currentShuffleField).(currentSignalName);
% 
%         numROIs = size(realActivity,1);
%         numShuffles = size(shuffleMatrix,3);
% 
%         %% Vectorized Range Calculation
%         % 1. Real Range
%         meanRealTuningCurve = squeeze(mean(realActivity, 2, 'omitnan')); % (numROIs x numBins)
%         realMaxValues = max(meanRealTuningCurve, [], 2, 'omitnan'); 
%         realMinValues = min(meanRealTuningCurve, [], 2, 'omitnan'); 
%         realRangeValues = realMaxValues - realMinValues; % (numROIs x 1)
% 
%         % 2. Shuffle Range Distribution
%         shuffleMaxValues = squeeze(max(shuffleMatrix, [], 2, 'omitnan')); % (numROIs x numShuffles)
%         shuffleMinValues = squeeze(min(shuffleMatrix, [], 2, 'omitnan')); % (numROIs x numShuffles)
%         shuffleRangeDistribution = shuffleMaxValues - shuffleMinValues; % (numROIs x numShuffles)
% 
%         %% Vectorized Significance and Percentile Calculations
%         % Calculate the 95th percentile threshold across all shuffles at once
%         rangeThreshold_95th = prctile(shuffleRangeDistribution, targetPercentile, 2); % (numROIs x 1)
% 
%         % Check significance via logical matrix comparison (True range > 95th threshold)
%         significanceIndex = double(realRangeValues > rangeThreshold_95th); % (numROIs x 1)
% 
%         % Count how many shuffle values fall below or equal the real range value for each ROI
%         numRangesBelowReal = sum(shuffleRangeDistribution <= realRangeValues, 2); % (numROIs x 1)
%         percentileValues = (numRangesBelowReal / numShuffles) * 100; % (numROIs x 1)
% 
%         %% Assign to the local output structure dynamically
%         nullDist_RangeTuningMetric.(shufLabel).isSignificantByRange.(currentSignalName) = significanceIndex;
%         nullDist_RangeTuningMetric.(shufLabel).realRangePercentileRank.(currentSignalName) = percentileValues;
% 
%         fprintf('  -> [%s] Found %d significant ROIs (Range > %dth Percentile) for %s.\n', ...
%             shufLabel, sum(significanceIndex), targetPercentile, currentSignalName);
%     end
% end
% 
% disp('Range significance calculation complete.');
% 
% %% Saving Data to sessionROIData 
% if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
%     disp(['Saving range significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
%     save(sessionFileInfo.otherSessFilePaths.sessionROIData, "nullDist_RangeTuningMetric", '-append');
% elseif isfield(sessionFileInfo, 'otherSessFilePaths')
%     warning('sessionROIData file not found at: %s. Cannot append range significance data.', ...
%         sessionFileInfo.otherSessFilePaths.sessionROIData);
% else
%     warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save range significance data.');
% end
% 
% end