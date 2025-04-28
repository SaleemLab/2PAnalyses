function plotRFResponses_ByPosition_AllROIs(sessionFileInfo, RFMappingStimName)
% Plots normalized per-ROI RF temporal responses, one page per stimulus position.
% Each page: small subplots for each ROI’s mean response at that position.

%% Load RF mapping response
iStim = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(iStim)
    error('Stimulus "%s" not found in sessionFileInfo.', RFMappingStimName);
end
respStruct = load(sessionFileInfo.stimFiles(iStim).Response, 'response');
RFResponse = respStruct.response;
psthData = RFResponse.psthData;

%% Prepare output PDF path
figDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
pdfName = sprintf('%s_%s_RF_by_position.pdf', sessionFileInfo.animal_name, RFMappingStimName);
pdfPath = fullfile(figDir, pdfName);

%% Get basic dimensions
nPositions = numel(psthData);
tVec = psthData(1).timeVector;
nROIs = size(psthData(1).alignedResponses, 1);

% Compute grid for ROI subplots
gridRows = ceil(sqrt(nROIs));
gridCols = ceil(nROIs/gridRows);

%% Loop over each stimulus position
for posIdx = 1:nPositions
    stimPos = psthData(posIdx).stimValue;
    % Compute per-ROI mean trace
    aligned = psthData(posIdx).alignedResponses;  % [ROIs x time x trials]
    meanByROI = squeeze(mean(aligned, 3, 'omitnan')); % [ROIs x time]

    % Normalize each ROI trace to [0 1]
    normTraces = nan(size(meanByROI));
    for r = 1:nROIs
        tr = meanByROI(r, :);
        if all(isnan(tr)), continue; end
        mn = min(tr);
        mx = max(tr);
        if mx > mn
            normTraces(r, :) = (tr - mn) / (mx - mn);
        else
            normTraces(r, :) = zeros(size(tr));
        end
    end

    % Create figure for this position
    fig = figure('Visible','off', 'Units','normalized','Position',[0 0 1 1]);
    t = tiledlayout(gridRows, gridCols, 'TileSpacing','compact','Padding','compact');
    title(t, sprintf('Stim pos [Az %d, El %d]', stimPos(1), stimPos(2)), 'FontSize', 16, 'FontWeight', 'bold');

    % Plot each ROI
    for r = 1:nROIs
        ax = nexttile;
        tr = normTraces(r, :);
        if ~all(isnan(tr))
            plot(ax, tVec, tr, 'k', 'LineWidth', 0.5);
        end
        xline(ax, 0, 'k--', 'LineWidth', 0.8);
        axis(ax, 'off');
        title(ax, sprintf('ROI %d', r), 'FontSize', 6);
    end

    % Export page to PDF
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType','vector');
    close(fig);
end
end
