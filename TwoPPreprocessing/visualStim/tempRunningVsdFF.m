%% Load example session 
load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25012\Analysis\20250507\M25012_20250507_processed2PData_M25012_GrayScreen_20250507_00002.mat")
load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25012\Analysis\20250507\M25012_20250507_PeripheralData_M25012_GrayScreen_20250507_00002.mat")

%% Note: computeAndExtractPeripheralData_GrayScreen.m; output can be saved in response and update sessionFileInfo
tickToCmConversion = 3.1415 * 20 / 1024;  % Wheel radius 20 cm, 1024 ticks per revolution
displacement = [0; diff(peripheralData.Wheel.Value * tickToCmConversion)];
% Handle unrealistic large changes (e.g., due to teleportation or resets)
displacement(displacement < -100) = 0;  % Negative large jumps
displacement(displacement > 100) = 0;   % Positive large jumps
% Calculate speed (in cm/s) 
wheelSpeed = displacement ./ [0; diff(peripheralData.Wheel.sampleTimes)];

%% All data interpolated to two-p frame times; labelled timeVec here 
timeVec = processedTwoPData.(processedTwoPData.resample2PTimeUsed);

%%
dFF = processedTwoPData.zScoredProcessedSignals.dFF;
stationaryIdx = wheelSpeed < 1;
movingIdx = wheelSpeed > 1;
% moving bins to represent 7% of the TOTAL data.
numTotalPoints = length(wheelSpeed);
pointsPerBin = floor(0.07 * numTotalPoints);
movingSpeeds = wheelSpeed(movingIdx);
% quantiles to find edges where each bin has equal data density
numMovingBins = floor(sum(movingIdx) / pointsPerBin);
if numMovingBins > 0
    % Calculate edges based on percentiles of the moving speed
    movingEdges = quantile(movingSpeeds, linspace(0, 1, numMovingBins + 1));
else
    movingEdges = [];
end
% Bin 1 is stationary, others are the moving quantiles
allEdges = [1, movingEdges(2:end)]; 
numFinalBins = length(allEdges) - 1;
%%
roiIdx = 310; 
roiActivity = dFF(roiIdx,:);
binMeanSpeed = zeros(numFinalBins, 1);
binMeanActivity = zeros(numFinalBins, 1);
binErrorActivity = zeros(numFinalBins, 1);

for thisBin = 1:numFinalBins
    % find indices for this specific speed bin
    idx = wheelSpeed >= allEdges(thisBin) & wheelSpeed < allEdges(thisBin+1);
    
    if any(idx)
        % speed in any bin was the mean speed during that time
        binMeanSpeed(thisBin) = mean(wheelSpeed(idx));
        
        % mean and sem
        binMeanActivity(thisBin) = mean(roiActivity(idx));
        binErrorActivity(thisBin) = std(roiActivity(idx)) / sqrt(sum(idx));
    end
end

plotIdx = 1:numFinalBins; 
plotSpeeds = binMeanSpeed(plotIdx);
plotActivity = binMeanActivity(plotIdx);
plotErrors = binErrorActivity(plotIdx);

figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 600 350]);
errorbar(plotSpeeds, plotActivity, plotErrors, 'ko-', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'CapSize', 0);

set(gca, 'XScale', 'log');

xMin = 1; 
xMax = max(plotSpeeds) * 1.1;
xlim([xMin, xMax]);

% Create ticks that "fill" the log space visually
potentialTicks = [5, 10, 20, 30, 40, 60, 100];
finalTicks = sort([1, potentialTicks(potentialTicks > 1 & potentialTicks < max(plotSpeeds)), max(plotSpeeds)]);

xticks(finalTicks);
xticklabels(string(round(finalTicks, 1)));

% Styling to match the paper (Clean axes, no top/right box)
xlabel('Running Speed (cm/s)');
ylabel('Mean (z-scored) DF/F ');
title(sprintf('ROI %d: Speed Tuning', roiIdx), 'FontSize', 13);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 13);

saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\RunningSpeed'; %
if ~exist(saveDir, 'dir')
    mkdir(saveDir); % Create the folder if it doesn't exist
end

fileName = sprintf('ROI_%d_SpeedTuning.png', roiIdx);
fullPath = fullfile(saveDir, fileName);

% Save at 300 DPI for high quality
exportgraphics(gcf, fullPath, 'Resolution', 300);