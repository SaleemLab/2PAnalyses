function [isSignificantByPeakShuffling, realPeakPercentileRank] = getPeakSignificance_fromShuffle(sessionFileInfo, response)
% Calculates the peak activity for each cell's spatial tuning curve,
% compares it to the distribution of peak activities from the shuffled data,
% and determines which cells have a peak activity greater than the 95th percentile.
%
% Aman and Sonali November 2025
%% Input Checks and Preparation
if nargin < 2
    error('Must provide the response struct.');
end
signalNames = fieldnames(response.lapPositionActivity);
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
isSignificantByPeakShuffling = isSignificantByPeakShuffling_local;
realPeakPercentileRank = realPeakPercentileRank_local;

% Check if the file path exists and save the variables
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    
    disp(['Saving peak significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
        "isSignificantByPeakShuffling", ...
        "realPeakPercentileRank", ...
         '-append')
         
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
end

end 