function [isSignificantByPeakShuffling_local, realPeakPercentileRank_local] = getPeakSignificance_fromShuffle(sessionFileInfo, response)
% Calculates the peak activity for each cell's spatial tuning curve,
% compares it to the distribution of peak activities from the shuffled data,
% and determines which cells have a peak activity greater than the 95th percentile.
%
% Aman and Sonali November 2025
%% Input Checks and Preparation
if nargin < 2
    error('Must provide the response struct.');
end
signalNames = {'dFF', 'dFFNeuropilCorrected'};
%signalNames = fieldnames(response.lapPositionActivity);
if ~isfield(response, 'lapPositionActivity_ShuffleMatrix')
    error('Shuffled data (lapPositionActivity_ShuffleMatrix) not found in response struct. Run getLapPositionActivity_withShuffle first.');
end
% Index of significant cells (1 = significant, 0 = not significant)
isSignificantByPeakShuffling_local = struct(); % Local container for saving
% Percentile of real peak relative to shuffle peaks
realPeakPercentileRank_local = struct(); % Local container for saving
% Target percentile for significance
targetPercentile = 95;
%% Main Processing Loop
disp('Starting peak significance calculation...');
for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};

    disp(['Processing signal: ', currentSignalName]);

    % Real activity (Neurons x Laps x Bins)
    realActivity = response.lapPositionActivity.(currentSignalName); 
    % Shuffled activity (Neurons x Bins x Shuffles)
    shuffleMatrix = response.lapPositionActivity_ShuffleMatrix.(currentSignalName);

    [numROIs, ~, ~] = size(realActivity);
    [~, ~, numShuffles] = size(shuffleMatrix);

    %Calculate the mean real tuning curve
    % Average activity across all laps to get the real tuning curve (Neurons x Bins)
    % This is the reference for the "real peak."
    meanRealTuningCurve = mean(realActivity, 2, 'omitnan'); % (numROIs x 1 x numBins; this is to keep sonali sane)
    meanRealTuningCurve = squeeze(meanRealTuningCurve); % (numROIs x numBins)

    % Find the real peak value for each neuron 
    % The peak is the maximum activity across all position bins
    realPeakValues = max(meanRealTuningCurve, [], 2, 'omitnan'); % (numROIs x 1)

    % Find the peak value for each shuffle, for each nuron 
    % Find the maximum bin activity (Peak) within each shuffle
    %  
    shufflePeakValues_temp = max(shuffleMatrix, [], 2, 'omitnan'); %(Neurons x 1 x Shuffles)
    shufflePeakValues = squeeze(shufflePeakValues_temp); % (numROIs x numShuffles)

    % Sig and percentiles 

    significanceIndex = zeros(numROIs, 1);
    percentileValues = nan(numROIs, 1);


    % Loop over each ROI to compare its real peak against the shuffle peak distribution
    for thisROI = 1:numROIs

        % Distribution of peak values for the current neuron across all 1000 shuffles
        currentShuffleDistribution = shufflePeakValues(thisROI, :);

        % The actual peak value for the current neuron
        currentRealPeak = realPeakValues(thisROI);

        % Calculate the target significance threshold (95th percentile)
        % This is the peak value that only 5% of the shuffled peaks exceed.
        peakThreshold_95th = prctile(currentShuffleDistribution, targetPercentile);

        % Check for significance: Is the Real Peak > 95th Percentile Threshold?
        if currentRealPeak > peakThreshold_95th
            significanceIndex(thisROI) = 1; % Mark as significant
        end

        % Calculate the percentile of the real peak relative to the shuffle distribution
        % Formula: (Number of shuffle peaks <= Real Peak) / Total Shuffles * 100

        % Count how many shuffled peaks are less than or equal to the real peak
        numPeaksBelowReal = sum(currentShuffleDistribution <= currentRealPeak);

        % Calculate the percentile
        percentile = (numPeaksBelowReal / numShuffles) * 100;

        percentileValues(thisROI) = percentile;
    end

    % Store the results in the local structs instead of the response struct
    isSignificantByPeakShuffling_local.(currentSignalName) = significanceIndex;
    realPeakPercentileRank_local.(currentSignalName) = percentileValues;

    fprintf('  -> Found %d significant ROIs (Peak > %dth Percentile) for %s.\n', ...
        sum(significanceIndex), targetPercentile, currentSignalName);

end
disp('Peak significance calculation complete.');

%% Saving Data to sessionROIData 

% Assign the local structs to the final variables used in the save command
nullDist_PeakTuningMetric.isSignificantByPeakShuffling = isSignificantByPeakShuffling_local;
nullDist_PeakTuningMetric.realPeakPercentileRank = realPeakPercentileRank_local;
nullDist_PeakTuningMetric.targetPercentile = targetPercentile; 

% Check if the file path exists and save the variables
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2

    disp(['Saving peak significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);

    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
       "nullDist_PeakTuningMetric", ...
         '-append')

elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
end

end 

% %% Testing with random permutations 
% function [nullDist_PeakTuningMetric] = getPeakSignificance_fromShuffle(sessionFileInfo, response)
% % Calculates the peak activity for each cell's spatial tuning curve,
% % compares it to the distribution of peak activities from the shuffled data,
% % and determines which cells have a peak activity greater than the 95th percentile.
% %
% % Updated to process both standard shuffle and RandPerm shuffle metrics.
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
% nullDist_PeakTuningMetric = struct();
% nullDist_PeakTuningMetric.targetPercentile = targetPercentile;
% 
% %% Main Processing Loop
% disp('Starting peak significance calculation...');
% 
% for iType = 1:length(shuffleTypesToRun)
%     currentShuffleField = shuffleTypesToRun{iType};
% 
%     % Create clean identifiers for storing results
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
%         numROIs = size(realActivity, 1);
%         numShuffles = size(shuffleMatrix, 3);
% 
%         %% Vectorized Peak Calculation
%         % 1. Real Peaks (Squeeze to numROIs x numBins to keep Sonali sane!)
%         meanRealTuningCurve = squeeze(mean(realActivity, 2, 'omitnan')); 
%         realPeakValues = max(meanRealTuningCurve, [], 2, 'omitnan'); % (numROIs x 1)
% 
%         % 2. Shuffle Peaks Distribution
%         shufflePeakValues = squeeze(max(shuffleMatrix, [], 2, 'omitnan')); % (numROIs x numShuffles)
% 
%         %% Vectorized Significance and Percentile Calculations
%         % Calculate the 95th percentile threshold across all shuffles at once
%         peakThreshold_95th = prctile(shufflePeakValues, targetPercentile, 2); % (numROIs x 1)
% 
%         % Check significance via logical matrix comparison (True peak > 95th threshold)
%         significanceIndex = double(realPeakValues > peakThreshold_95th); % (numROIs x 1)
% 
%         % Count how many shuffle peaks fall below or equal the real peak value for each ROI
%         numPeaksBelowReal = sum(shufflePeakValues <= realPeakValues, 2); % (numROIs x 1)
%         percentileValues = (numPeaksBelowReal / numShuffles) * 100; % (numROIs x 1)
% 
%         %% Assign to the local output structure dynamically
%         nullDist_PeakTuningMetric.(shufLabel).isSignificantByPeakShuffling.(currentSignalName) = significanceIndex;
%         nullDist_PeakTuningMetric.(shufLabel).realPeakPercentileRank.(currentSignalName) = percentileValues;
% 
%         fprintf('  -> [%s] Found %d significant ROIs (Peak > %dth Percentile) for %s.\n', ...
%             shufLabel, sum(significanceIndex), targetPercentile, currentSignalName);
%     end
% end
% 
% disp('Peak significance calculation complete.');
% 
% %% Saving Data to sessionROIData 
% if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
%     disp(['Saving peak significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
%     save(sessionFileInfo.otherSessFilePaths.sessionROIData, "nullDist_PeakTuningMetric", '-append');
% elseif isfield(sessionFileInfo, 'otherSessFilePaths')
%     warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
%         sessionFileInfo.otherSessFilePaths.sessionROIData);
% else
%     warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
% end
% 
% end