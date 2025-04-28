function plotRFTemporalTracesByPosition(sessionFileInfo, RFMappingStimName, azRange, elRange)
% Plots all ROI traces for each azimuth/elevation in the target range.
% Shows tiny normalised traces in a tiled layout.

%% Load data
iStim = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
RFResponseStruct = load(sessionFileInfo.stimFiles(iStim).Response, 'response');
RFResponse = RFResponseStruct.response;
psthData = RFResponse.psthData;

%% Get azimuth and elevation grid
stimValues = vertcat(psthData.stimValue);
azimuth = unique(stimValues(:, 1));
elevation = unique(stimValues(:, 2));

azMask = azimuth >= azRange(1) & azimuth <= azRange(2);
elMask = elevation >= elRange(1) & elevation <= elRange(2);

targetAz = azimuth(azMask);
targetEl = elevation(elMask);

%% Output setup
fig = figure('Position', [100 100 1400 800]);
t = tiledlayout(length(targetEl), length(targetAz), 'TileSpacing', 'compact', 'Padding', 'compact');
t.Title.String = sprintf('%s | RF Normalised Tiny Traces (%d ROIs)', sessionFileInfo.animal_name, size(psthData(1).meanResponse, 1));
t.Title.FontSize = 16;

nROIs = size(psthData(1).meanResponse, 1);

%% Loop over stimulus positions
for elIdx = 1:length(targetEl)
    for azIdx = 1:length(targetAz)
        az = targetAz(azIdx);
        el = targetEl(elIdx);

        % Find corresponding entry in psthData
        matchIdx = find(arrayfun(@(x) isequal(x.stimValue, [az el]), psthData));
        if isempty(matchIdx), continue; end

        d = psthData(matchIdx);
        tVec = d.timeVector;
        aligned = d.alignedResponses;  % [nROIs x time x trials]

        ax = nexttile;
        hold on;
        title(sprintf('Az %d | El %d', az, el), 'FontSize', 8);
        xline(0, 'k--', 'LineWidth', 0.8);

        for roi = 1:nROIs
            this = squeeze(aligned(roi, :, :));
            if size(this, 1) > size(this, 2), this = this'; end
            if isempty(this) || size(this, 2) ~= length(tVec), continue; end

            meanTrace = mean(this, 1, 'omitnan');
            if all(isnan(meanTrace)), continue; end

            % Normalise trace for compact viewing
            meanTrace = meanTrace - min(meanTrace);
            maxVal = max(meanTrace);
            if maxVal > 0
                meanTrace = meanTrace / maxVal;
            end

            % Offset in Y slightly for plotting multiple if needed
            plot(tVec, meanTrace, 'k-', 'LineWidth', 0.5);
        end

        xlim([min(tVec), max(tVec)]);
        ylim([0 1]);
        set(gca, 'FontSize', 6);
    end
end

%% Save
saveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
fname = sprintf('%s_RFNormalisedTinyTraces_Az%dto%d_El%dto%d.pdf', ...
    sessionFileInfo.animal_name, azRange(1), azRange(2), elRange(1), elRange(2));
exportgraphics(fig, fullfile(saveDir, fname), 'ContentType', 'vector');
close(fig);

end
