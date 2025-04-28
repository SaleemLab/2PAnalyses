function plotSpatialRFforROI(psthData, roiIdx, applySmoothing)
% plotSingleRFMapFromPSTH - Plot spatial RF heatmap for one ROI using psthData
%
% Inputs:
%   - psthData: output from getTrialResponsePSTHsV4
%   - roiIdx: ROI index to plot
%   - applySmoothing: true/false to apply Gaussian smoothing

if nargin < 3
    applySmoothing = false;
end

% Extract unique azimuth/elevation values from stimValue
stimValues = vertcat(psthData.stimValue);  % Nx2
azimuth = unique(stimValues(:, 1));
elevation = unique(stimValues(:, 2));

% Preallocate RF map
RFmap = nan(length(elevation), length(azimuth));

% Fill RFmap
for i = 1:length(psthData)
    stimPos = psthData(i).stimValue;  % [az, el]
    stimAz = stimPos(1);
    stimEl = stimPos(2);

    azIdx = find(azimuth == stimAz);
    elIdx = find(elevation == stimEl);

    if isempty(azIdx) || isempty(elIdx)
        warning('Stimulus position [%d %d] not found in grid.', stimAz, stimEl);
        continue;
    end

    roiTrace = psthData(i).meanResponse(roiIdx, :);  % 1 x T
 
    timeVec = psthData(i).timeVector;
    timeMask = timeVec >= 0 & timeVec <= 0.5;
    meanResp = mean(roiTrace(timeMask), 'omitnan');

    RFmap(elIdx, azIdx) = meanResp;
end

% Optional smoothing
if applySmoothing
    RFmap = imgaussfilt(RFmap, 1);  % 1-pixel Gaussian
end

% Plot
figure;
imagesc(azimuth, elevation, RFmap);
set(gca, 'YDir', 'normal', ...
         'TickDir', 'out', ...
         'FontSize', 14, ...
         'Box', 'off', ...
         'Color', 'none');
xline(0, 'k--', 'LineWidth', 1.2);
yline(0, 'k--', 'LineWidth', 1.2);
xlabel('Azimuth (\circ)');
ylabel('Elevation (\circ)');
title(sprintf('RF map - ROI %d', roiIdx));
colormap(redWhiteBlue);
colorbar;

% Output raw map to console
disp(['RFmap (raw) for ROI ' num2str(roiIdx) ':']);
disp(RFmap);
end
