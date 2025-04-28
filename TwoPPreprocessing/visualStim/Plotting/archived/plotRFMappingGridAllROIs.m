function plotRFMappingGridAllROIs(sessionFileInfo, RFMappingStimName, FOV_label)
% Plots RF line profile maps (not heatmaps) for all ROIs in a tiled layout.
% Each page contains ~100 ROIs with traces in a grid matching azimuth × elevation.

%% Load response data
iStim = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
RFResponseStruct = load(sessionFileInfo.stimFiles(iStim).Response, 'response');
RFResponse = RFResponseStruct.response;
psthData = RFResponse.psthData;

%% Azimuth & Elevation grid
stimValues = vertcat(psthData.stimValue);
azimuth = unique(stimValues(:, 1));
elevation = unique(stimValues(:, 2));

%% Output setup
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pdfName = sprintf('%s_%s_RFMapping_LineProfiles.pdf', sessionFileInfo.animal_name, sessionFileInfo.session_name);
pdfPath = fullfile(figSaveDir, pdfName);

nROIs = size(psthData(1).meanResponse, 1);
plotsPerPage = 100;

for startIdx = 1:plotsPerPage:nROIs
    endIdx = min(startIdx + plotsPerPage - 1, nROIs);
    currentROIs = startIdx:endIdx;

    % Create figure and tile layout
    fig = figure('Visible', 'off', 'Position', [100 100 1800 1200]);
    t = tiledlayout(10, 10, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Use sgtitle for correct rendering
    sgtitle(sprintf('%s | FOV %s | ROIs %d to %d', ...
        sessionFileInfo.animal_name, FOV_label, startIdx, endIdx), ...
        'FontSize', 16, 'FontWeight', 'bold');

    for roiIdx = currentROIs
        fprintf('Plotting ROI %d...\n', roiIdx);

        % Collect line traces per RF position
        RFTraces = cell(length(elevation), length(azimuth));
        timeVec = psthData(1).timeVector;

        for i = 1:length(psthData)
            stimPos = psthData(i).stimValue;
            azIdx = find(azimuth == stimPos(1));
            elIdx = find(elevation == stimPos(2));
            if isempty(azIdx) || isempty(elIdx), continue; end

            aligned = squeeze(psthData(i).alignedResponses(roiIdx, :, :));
            if size(aligned, 1) > size(aligned, 2), aligned = aligned'; end
            meanTrace = mean(aligned, 1, 'omitnan');
            semTrace = std(aligned, 0, 1, 'omitnan') ./ sqrt(size(aligned, 1));
            RFTraces{elIdx, azIdx} = struct('mean', meanTrace, 'sem', semTrace);
        end

        ax = nexttile;
        hold on;
        for el = 1:length(elevation)
            for az = 1:length(azimuth)
                traceData = RFTraces{el, az};
                if isempty(traceData), continue; end

                offsetX = (az - 1) * (max(timeVec) + 0.1);
                offsetY = (length(elevation) - el) * 1.5;

                m = traceData.mean;
                s = traceData.sem;
                fill([timeVec fliplr(timeVec)] + offsetX, ...
                     [m + s, fliplr(m - s)] + offsetY, ...
                     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                plot(timeVec + offsetX, m + offsetY, 'k', 'LineWidth', 1);
                xline(offsetX, 'k:', 'Alpha', 0.3);
            end
        end
        axis tight off;
        title(sprintf('ROI %d', roiIdx), 'FontSize', 7);
    end

    % Add overall label
    annotation(fig, 'textbox', [0.4 0.01 0.2 0.03], ...
        'String', 'RF Temporal Response Profiles', ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'EdgeColor', 'none');

    % Save page
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');
    close(fig);
end

end
