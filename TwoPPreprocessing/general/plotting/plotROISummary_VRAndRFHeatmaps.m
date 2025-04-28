function plotROISummary_VRAndRFHeatmaps(sessionFileInfo, VRCorrStimName, RFMappingStimName, applySmoothing)
% Plot VR heatmaps + RF peak maps for each ROI individually. Metric used
% for VR and RF is median. 

if nargin < 4
    applySmoothing = true;
end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'PerROI_VRPSTHandHeatmaps_PeakRFGrid');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

%% Load VR data
isVR = find(strcmp(VRCorrStimName, {sessionFileInfo.stimFiles.name}), 1);
VRResponse = load(sessionFileInfo.stimFiles(isVR).Response, 'response');

%% Load RF mapping data
isRF = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
RFResponse = load(sessionFileInfo.stimFiles(isRF).Response, 'response');
RFpsthData = RFResponse.psthData;

%% Build azimuth / elevation grids
stimValues = vertcat(RFpsthData.stimValue);
azimuth = unique(stimValues(:,1));
elevation = unique(stimValues(:,2));

nAz = numel(azimuth);
nEl = numel(elevation);

nROIs = size(VRResponse.lapPositionActivity, 1);

% Parameters for RF calculation
prestimWindow = [-0.4 0]; % seconds
poststimWindow = [0 1]; % seconds to calculate the peak; currenly only including the peri-stimulus period. 
tVec = RFpsthData(1).timeVector(:)'; % assume all PSTH share same time vector
preMask = tVec >= prestimWindow(1) & tVec <= prestimWindow(2);
postMask = tVec >= poststimWindow(1) & tVec <= poststimWindow(2);

%% Compute RF peak value maps per ROI
RFmaps = nan(nEl, nAz, nROIs); % elevation x azimuth x ROI

for roiIdx = 1:nROIs
    tempRF = nan(nEl, nAz);
    for iStim = 1:numel(RFpsthData)
        stimPos = RFpsthData(iStim).stimValue;
        azIdx = find(azimuth == stimPos(1));
        elIdx = find(elevation == stimPos(2));
        if isempty(azIdx) || isempty(elIdx)
            continue;
        end
        % Extract this ROI's PSTH trace
        roiTrace = RFpsthData(iStim).alignedResponses(roiIdx,:,:);
        roiTrace = squeeze(roiTrace);
        if size(roiTrace,1)==1
            roiTrace = roiTrace';
        end
        % Changed to median 
        medianTrace = median(roiTrace,2,'omitnan');
        F0 = median(medianTrace(preMask),'omitnan');
        zn = medianTrace - F0;
        peakResp = max(zn(postMask),[],'omitnan');
        tempRF(elIdx, azIdx) = peakResp;
    end
    RFmaps(:,:,roiIdx) = tempRF;
end

%% Smooth VR lap activity if requested
lapActivityFull = VRResponse.lapPositionActivity;
if applySmoothing
    w = gausswin(9);
    w = w / sum(w);
    for iCell = 1:nROIs
        for iLap = 1:size(lapActivityFull, 2)
            trace = squeeze(lapActivityFull(iCell,iLap,:));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapActivityFull(iCell,iLap,:) = smoothed;
        end
    end
end

%% Loop through ROIs to plot
for roiIdx = 1:nROIs
    roiActivity = squeeze(lapActivityFull(roiIdx,:,:)); % laps x bins
    if all(isnan(roiActivity),'all'), continue; end
    normActivity = normalize(roiActivity, 2, 'range');

    RFmap = RFmaps(:,:,roiIdx);

    figure('Position', [100 100 1100 400]);

    % --- VR Panel
    subplot(1,2,1);
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

    % --- RF Panel
    subplot(1,2,2);
    imagesc(azimuth, elevation, RFmap);
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 12, 'Box', 'off');
    xline(0, 'k--', 'LineWidth', 1.2);
    yline(0, 'k--', 'LineWidth', 1.2);
    xlabel('Azimuth (°)');
    ylabel('Elevation (°)');
    title(sprintf('RF Peak Map - ROI %d', roiIdx));
    colormap(gca, turbo);
    colorbar;

    % --- Save
    figName = sprintf('%s_%s_ROI%d_VR_RFpeak.pdf', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, roiIdx);
    figPath = fullfile(figSaveDir, figName);
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperPosition', [0 0 11 4], ...
             'PaperOrientation', 'landscape');
    print(gcf, figPath, '-dpdf', '-r300');
    close(gcf);
end

end
