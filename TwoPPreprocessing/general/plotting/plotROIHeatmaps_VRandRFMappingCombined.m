function plotROIHeatmaps_VRandRFMappingCombined(sessionFileInfo, VRCorrStimName, RFMappingStimName, applySmoothing)

if nargin < 4
    applySmoothing = false;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'PerROIHeatmaps_AllROIs');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

%% Load VR data
iVR = find(strcmp(VRCorrStimName, {sessionFileInfo.stimFiles.name}), 1);
load(sessionFileInfo.stimFiles(iVR).Response, 'response');
VRResponse = response;

%% Load RF mapping data
iRF = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
load(sessionFileInfo.stimFiles(iRF).Response, 'response');
load(sessionFileInfo.stimFiles(iRF).BonsaiData, 'bonsaiData');
RFResponse = response;

psthFile = fullfile(sessionFileInfo.Directories.save_folder, ...
    sprintf('%s_%s_PSTH_%s.mat', sessionFileInfo.animal_name, sessionFileInfo.session_name, RFMappingStimName));
load(psthFile, 'psthData');

%% VR heatmap: pre-smoothing
lapActivityFull = VRResponse.lapPositionActivity;
nROIs = size(lapActivityFull, 1);

if applySmoothing
    w = gausswin(9); w = w / sum(w);
    for iCell = 1:nROIs
        for iLap = 1:size(lapActivityFull, 2)
            trace = squeeze(lapActivityFull(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapActivityFull(iCell, iLap, :) = smoothed;
        end
    end
end

%% Azimuth / Elevation grid for RF mapping
stimValues = vertcat(psthData.stimValue);
azimuth = unique(stimValues(:, 1));
elevation = unique(stimValues(:, 2));

%% Loop through ROIs
for roiIdx = 1:nROIs
    roiActivity = squeeze(lapActivityFull(roiIdx, :, :)); % laps x bins
    if all(isnan(roiActivity), 'all'), continue; end
    normActivity = normalize(roiActivity, 2, 'range');

    % --- Build RF map
    RFmap = nan(length(elevation), length(azimuth));
    for i = 1:length(psthData)
        stimPos = psthData(i).stimValue;
        azIdx = find(azimuth == stimPos(1));
        elIdx = find(elevation == stimPos(2));
        if isempty(azIdx) || isempty(elIdx), continue; end

        roiTrace = psthData(i).meanResponse(roiIdx, :);
        timeVec = psthData(i).timeVector;
        timeMask = timeVec >= 0 & timeVec <= 0.5;
        meanResp = mean(roiTrace(timeMask), 'omitnan');
        RFmap(elIdx, azIdx) = meanResp;
    end

    if applySmoothing
        RFmap = imgaussfilt(RFmap, 1);
    end

    % --- Plot both panels in a single figure
    figure('Position', [100 100 1100 400]);

    subplot(1, 2, 1);
    imagesc(normActivity);
    caxis([0 1]); colormap(gca, flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlabel('Position (cm)');
    ylabel('Lap #');
    title(sprintf('VR Activity - ROI %d', roiIdx));
    colorbar;

    subplot(1, 2, 2);
    imagesc(azimuth, elevation, RFmap);
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 12, 'Box', 'off');
    xline(0, 'k--', 'LineWidth', 1.2);
    yline(0, 'k--', 'LineWidth', 1.2);
    xlabel('Azimuth (°)');
    ylabel('Elevation (°)');
    title(sprintf('RF Map - ROI %d', roiIdx));
    colormap(gca, turbo);
    colorbar;

    % --- Save
    figName = sprintf('%s_%s_ROI%d_VR_RF.pdf', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, roiIdx);
    figPath = fullfile(figSaveDir, figName);
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperPosition', [0 0 11 4], ...
             'PaperOrientation', 'landscape');
    print(gcf, figPath, '-dpdf', '-r300');
    close(gcf);
end
end
