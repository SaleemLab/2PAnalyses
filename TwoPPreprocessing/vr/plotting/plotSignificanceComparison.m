function plotSignificanceComparison(sessionFileInfo, signalName, roiIndex, targetPercentile)
% Plots the distribution of shuffled peak activities and compares it 
% to the real peak activity for a single specified ROI.
%
% INPUTS:
% sessionFileInfo: Struct containing the path to sessionROIData
% signalName: The name of the signal (e.g., 'dFF', 'dFFNeuropilCorrected')
% roiIndex: The 1-based index of the ROI to plot (e.g., 1 to 450)
% targetPercentile: The significance threshold (default: 95)
%
% Usage Example (assuming signalName = 'dFF'):
% plotSignificanceComparison(sessionFileInfo, 'dFF', 55, 95); 

%% 1. Load Data
if nargin < 4, targetPercentile = 95; end

% Load the necessary variables from the sessionROIData file
roiDataFile = sessionFileInfo.otherSessFilePaths.sessionROIData;

if exist(roiDataFile, 'file') == 2
    % Load the calculated significance data
    loadedVars = load(roiDataFile, 'isSignificantByPeakShuffling', 'realPeakPercentileRank');
else
    error('sessionROIData file not found at: %s. Cannot plot results.', roiDataFile);
end

% You also need the original shuffle matrix to reconstruct the shuffle peaks
% It is typically saved in the Response struct, so we must load the Response file.
stimIdx = find(strcmp(signalName, fieldnames(loadedVars.isSignificantByPeakShuffling)));
if isempty(stimIdx), error('Signal name not found in loaded significance data.'); end

% We need to access the response struct to get the raw shuffle matrix data
stimIdx = find(strcmp(signalName, fieldnames(loadedVars.isSignificantByPeakShuffling)));
if isempty(stimIdx), error('Signal name not found in loaded significance data.'); end

% Assuming the response struct is stored in a separate file (e.g., Response.mat)
% You'll need the path to the Response file to load the raw shuffle matrix.
% (This path must be available in sessionFileInfo, often under stimFiles)
try
    stimFileIdx = find(strcmp(signalName, {sessionFileInfo.stimFiles.name}));
    responseFilePath = sessionFileInfo.stimFiles(stimFileIdx).Response;
    responseStruct = load(responseFilePath, 'lapPositionActivity_ShuffleMatrix');
    
    % Get the shuffle matrix: (Neurons x Bins x Shuffles)
    shuffleMatrix = responseStruct.lapPositionActivity_ShuffleMatrix.(signalName);
catch
     error('Could not load the raw lapPositionActivity_ShuffleMatrix from the expected Response file.');
end

numROIs = size(shuffleMatrix, 1);
if roiIndex > numROIs || roiIndex < 1
    error('Invalid ROI index specified. Max ROI index is %d.', numROIs);
end

%% 2. Calculate Peak Distributions

% Get the shuffle matrix for the selected ROI: (Bins x Shuffles)
roiShuffleMatrix = squeeze(shuffleMatrix(roiIndex, :, :));

% Calculate the PEAK VALUE for each of the 1000 shuffles (1 x Shuffles)
shufflePeakDistribution = max(roiShuffleMatrix, [], 1); 

% Calculate the REAL PEAK value
% (Need to recalculate the real peak as it was not saved as a stand-alone variable)
% Assuming real activity is available in a similar path (or you can pass it as an argument)
try
    responseStructLapAct = load(responseFilePath, 'lapPositionActivity');
    realActivity = responseStructLapAct.lapPositionActivity.(signalName);
    meanRealTuningCurve = squeeze(mean(realActivity(roiIndex, :, :), 2, 'omitnan'));
    realPeakValue = max(meanRealTuningCurve);
catch
     error('Could not load the raw lapPositionActivity from the expected Response file.');
end

% Get the significance threshold (95th percentile)
peakThreshold = prctile(shufflePeakDistribution, targetPercentile);

% Get the final verdict and percentile rank
isSig = loadedVars.isSignificantByPeakShuffling.(signalName)(roiIndex);
rankVal = loadedVars.realPeakPercentileRank.(signalName)(roiIndex);


%% 3. Plotting the Results
figure('Name', ['ROI ', num2str(roiIndex), ' - Peak Significance - ', signalName]);

% A. Histogram of Shuffled Peaks
h = histogram(shufflePeakDistribution, 'Normalization', 'probability');
xlabel('Peak Activity (Shuffle Distribution)');
ylabel('Probability');

% B. Significance Threshold (Red Dashed Line)
hold on;
yMax = max(h.Values);
line([peakThreshold peakThreshold], [0 yMax], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--', ...
     'DisplayName', [num2str(targetPercentile), 'th Percentile Threshold']);

% C. Real Peak Value (Green Solid Line)
line([realPeakValue realPeakValue], [0 yMax], 'Color', [0 0.5 0], 'LineWidth', 3, ...
     'DisplayName', ['Real Peak (Rank: ', num2str(round(rankVal, 1)), '%)']);

% D. Title and Legend
if isSig == 1
    titleStr = ['ROI ', num2str(roiIndex), ' is SIGNIFICANT (Real Peak > ', num2str(targetPercentile), 'th P)'];
else
    titleStr = ['ROI ', num2str(roiIndex), ' is NOT Significant'];
end
title(titleStr, 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;

end