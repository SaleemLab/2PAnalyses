function plotSingleROIHeatmap_VRandRF(sessionFileInfo, VRCorrStimName, RFMappingStimName, roiIdx, applySmoothing)

if nargin < 5
    applySmoothing = false;
end

%% Load VR data
iVR = find(strcmp(VRCorrStimName, {sessionFileInfo.stimFiles.name}), 1);
VRResponseStruct = load(sessionFileInfo.stimFiles(iVR).Response, 'response');
VRResponse = VRResponseStruct.response;

%% Load RF Mapping data
iRF = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
RFResponseStruct = load(sessionFileInfo.stimFiles(iRF).Response, 'response');
RFResponse = RFResponseStruct.response;
psthData = RFResponse.psthData;

%% --- Prepare VR lap activity
lapActivityFull = VRResponse.lapPositionActivity;
roiActivity = squeeze(lapActivityFull(roiIdx, :, :));

if applySmoothing
    w = gausswin(9); w = w / sum(w);
    for iLap = 1:size(roiActivity, 1)
        trace = roiActivity(iLap, :);
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace);
        trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace);
        smoothed(nanMask) = NaN;
        roiActivity(iLap, :) = smoothed;
    end
end

normActivity = normalize(roiActivity, 2, 'range');
meanActivity = mean(roiActivity, 1, 'omitnan');
semActivity = std(roiActivity, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiActivity), 1));

%% --- Prepare RF map & traces
stimValues = vertcat(psthData.stimValue);
azimuth = unique(stimValues(:, 1));
elevation = unique(stimValues(:, 2));
RFmap = nan(length(elevation), length(azimuth));
meanRFTraces = cell(length(elevation), length(azimuth));
semRFTraces = cell(length(elevation), length(azimuth));

for i = 1:length(psthData)
    stimPos = psthData(i).stimValue;
    azIdx = find(azimuth == stimPos(1));
    elIdx = find(elevation == stimPos(2));
    if isempty(azIdx) || isempty(elIdx), continue; end

    roiTrace = psthData(i).alignedResponses(roiIdx, :, :);
    roiTrace = squeeze(roiTrace);
    if size(roiTrace, 1) > size(roiTrace, 2)
        roiTrace = roiTrace';  % Ensure Trials x Time
    end
    if isempty(roiTrace) || size(roiTrace, 2) ~= length(psthData(i).timeVector), continue; end

    tVec = psthData(i).timeVector;
    timeMask = tVec >= 0 & tVec <= 0.5;
    if any(timeMask)
        meanResp = mean(mean(roiTrace(:, timeMask), 2, 'omitnan'), 'omitnan');
        RFmap(elIdx, azIdx) = meanResp;
    end

    meanRFTraces{elIdx, azIdx} = mean(roiTrace, 1, 'omitnan');
    semRFTraces{elIdx, azIdx} = std(roiTrace, 0, 1, 'omitnan') ./ sqrt(size(roiTrace, 1));
end

%% --- Plot layout: 2 rows x 2 columns
fig = figure('Position', [100 100 1600 900]);

% -- 1. Mean ± SEM (VR)
subplot(2, 2, 1);
hold on;
x = 1:size(meanActivity, 2);
fill([x fliplr(x)], [meanActivity + semActivity, fliplr(meanActivity - semActivity)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(x, meanActivity, 'k', 'LineWidth', 2);
xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
xticks([0 50 70 90 110 140]);
xticklabels({'0', '50', '70', '90', '110', '140'});
xlabel('Position (cm)');
ylabel('Mean Delta F/F');
title(sprintf('%s - ROI %d (VR)', sessionFileInfo.animal_name, roiIdx));
box off; set(gca, 'FontSize', 12);

% -- 2. Lap-by-position heatmap (VR)
subplot(2, 2, 2);
imagesc(normActivity);
caxis([0 1]); colormap(gca, flipud(gray));
xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
xticks([0 50 70 90 110 140]);
xticklabels({'0', '50', '70', '90', '110', '140'});
xlabel('Position (cm)');
ylabel('Lap #');
title('Lap-by-position activity');
colorbar; ylabel(colorbar, 'deltaf/f (normalised)');
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');

% -- 3. RF Spatial Map
subplot(2, 2, 3);
imagesc(azimuth, elevation, RFmap);
set(gca, 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 12, 'Box', 'off');
xline(0, 'k--'); yline(0, 'k--');
xlabel('Azimuth (°)');
ylabel('Elevation (°)');
title(sprintf('RF Map - ROI %d', roiIdx));
colormap(gca, redWhiteBlue);
colorbar;

% -- 4. RF Temporal Profiles Grid
subplot(2, 2, 4); hold on;
tVec = psthData(1).timeVector;
for el = 1:length(elevation)
    for az = 1:length(azimuth)
        meanTrace = meanRFTraces{el, az};
        semTrace = semRFTraces{el, az};
        if isempty(meanTrace) || length(meanTrace) ~= length(tVec), continue; end

        % Offset placement
        offsetX = (az - 1) * (max(tVec) + 0.2);
        offsetY = (length(elevation) - el) * 1.5;

        % Plot SEM fill + mean line
        fill([tVec fliplr(tVec)] + offsetX, ...
             [meanTrace + semTrace, fliplr(meanTrace - semTrace)] + offsetY, ...
             [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot(tVec + offsetX, meanTrace + offsetY, 'k', 'LineWidth', 1.2);
        xline(offsetX, 'k:', 'Alpha', 0.4);
    end
end
xlim([0, offsetX + max(tVec)]);
xlabel('Time (s)');
ylabel('Offset traces by Elevation');
title('RF Temporal Profiles (Grid Layout)');
set(gca, 'FontSize', 12);

%% Save to PDF
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end

pdfName = sprintf('%s_%s_ROI%d_LapPositionActivity_and_RFMapping.pdf', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, roiIdx);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [16 9]);
set(gcf, 'PaperPosition', [0 0 16 9]);
print(gcf, fullfile(figSaveDir, pdfName), '-dpdf', '-r300');
close(gcf);

end
