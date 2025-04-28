function plotAllROIHeatmaps_VRandRFV2(sessionFileInfo, VRCorrStimName, RFMappingStimName, applySmoothing)

if nargin < 4
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

% Get psthData from inside RF response
psthData = RFResponse.psthData;

%% Prepare figure output
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end

pdfName = sprintf('%s_%s_AllROIs_LapPositionActivity_and_RFMapping.pdf', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name);
pdfPath = fullfile(figSaveDir, pdfName);

lapActivityFull = VRResponse.lapPositionActivity;
nROIs = size(lapActivityFull, 1);

stimValues = vertcat(psthData.stimValue);
azimuth = unique(stimValues(:, 1));
elevation = unique(stimValues(:, 2));

for roiIdx = 1:nROIs
    roiActivity = squeeze(lapActivityFull(roiIdx, :, :));
    if all(isnan(roiActivity), 'all'), continue; end

    % Smoothing
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

    % RF Map
    RFmap = nan(length(elevation), length(azimuth));
    trialTraces = cell(length(elevation), length(azimuth));
    meanRFTraces = cell(length(elevation), length(azimuth));
    semRFTraces = cell(length(elevation), length(azimuth));
    for i = 1:length(psthData)
        stimPos = psthData(i).stimValue;
        azIdx = find(azimuth == stimPos(1));
        elIdx = find(elevation == stimPos(2));
        if isempty(azIdx) || isempty(elIdx), continue; end

        roiTrace = psthData(i).alignedResponses(roiIdx, :, :);
        roiTrace = squeeze(roiTrace);  % trials x time
        if isempty(roiTrace), continue; end

        tVec = psthData(i).timeVector;
        if isempty(tVec), continue; end

        trialTraces{elIdx, azIdx} = roiTrace;

        % Truncate to min size to avoid indexing errors
        minLength = min(length(tVec), size(roiTrace, 2));
        roiTrace = roiTrace(:, 1:minLength);
        tVec = tVec(1:minLength);

        timeMask = tVec >= 0 & tVec <= 0.5;
        meanResp = mean(mean(roiTrace(:, timeMask), 2, 'omitnan'), 'omitnan');
        RFmap(elIdx, azIdx) = meanResp;

        meanRFTraces{elIdx, azIdx} = mean(roiTrace, 1, 'omitnan');
        semRFTraces{elIdx, azIdx} = std(roiTrace, [], 1, 'omitnan') ./ sqrt(size(roiTrace, 1));
    end

    % --- Plot
    fig = figure('Visible', 'off', 'Position', [100 100 1600 900]);

    % Mean ± SEM (VR)
    subplot(2, 2, 1); hold on;
    x = 1:size(meanActivity, 2);
    fill([x fliplr(x)], [meanActivity + semActivity, fliplr(meanActivity - semActivity)], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, meanActivity, 'k', 'LineWidth', 2);
    xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
    xticks([0 50 70 90 110 140]); xticklabels({'0','50','70','90','110','140'});
    xlabel('Position (cm)'); ylabel('Mean Delta F/F');
    title(sprintf('%s - ROI %d (VR)', sessionFileInfo.animal_name, roiIdx));
    box off; set(gca, 'FontSize', 12);

    % Heatmap (VR)
    subplot(2, 2, 2);
    imagesc(normActivity);
    caxis([0 1]); colormap(gca, flipud(gray));
    xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
    xticks([0 50 70 90 110 140]); xticklabels({'0','50','70','90','110','140'});
    xlabel('Position (cm)'); ylabel('Lap #');
    title('Lap-by-position activity');
    colorbar; ylabel(colorbar, 'deltaf/f (normalised)');
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');

    % RF Map
    subplot(2, 2, 3);
    imagesc(azimuth, elevation, RFmap);
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 12, 'Box', 'off');
    xline(0, 'k--'); yline(0, 'k--');
    xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
    title(sprintf('RF Map - ROI %d', roiIdx));
    colormap(gca, redWhiteBlue);
    colorbar;

    % RF Traces Grid
    subplot(2, 2, 4); hold on;
    tVec = psthData(1).timeVector;
    tLen = length(tVec);
    for el = 1:length(elevation)
        for az = 1:length(azimuth)
            meanTrace = meanRFTraces{el, az};
            semTrace = semRFTraces{el, az};
            if isempty(meanTrace) || length(meanTrace) ~= tLen, continue; end

            offsetX = (az - 1) * (max(tVec) + 0.1);
            offsetY = (length(elevation) - el) * 1.5;

            fill([tVec fliplr(tVec)] + offsetX, ...
                [meanTrace + semTrace, fliplr(meanTrace - semTrace)] + offsetY, ...
                [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            plot(tVec + offsetX, meanTrace + offsetY, 'k');
            xline(offsetX + 0, 'k:', 'Alpha', 0.4);
        end
    end
    title('RF Timecourses by Grid');
    xlabel('Time (s) per location'); ylabel('Offset traces');
    set(gca, 'FontSize', 12);

    % Save
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', [16 9]);
    set(fig, 'PaperPosition', [0 0 16 9]);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'Append', true);
    close(fig);
end

end
