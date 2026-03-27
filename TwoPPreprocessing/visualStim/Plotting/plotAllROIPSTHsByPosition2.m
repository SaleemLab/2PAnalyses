function plotAllROIPSTHsByPosition2(sessionFileInfo, response, timeWindow, baselineWindow, peakTimeWindow)
% plotAllROIPSTHsByPosition  Loop over every ROI and save one page per ROI
%   into a single multi-page PDF. Tiles ROI's normalized, per-trial 
%   baseline-subtracted PSTHs across spatial positions.

%% 1. Defaults and Data Extraction
if nargin<3 || isempty(timeWindow),      timeWindow = [0 2]; end
if nargin<4 || isempty(baselineWindow),  baselineWindow = [-0.75 0]; end
if nargin<5 || isempty(peakTimeWindow),  peakTimeWindow = [0.25 2]; end

psth   = response.psthData;
nPos   = numel(psth);
stimVs = vertcat(psth.stimValue);

% Extract unique coordinates for the grid layout
azimuth   = sort(unique(stimVs(:,1)), 'descend'); % Visual field: Left is negative
elevation = sort(unique(stimVs(:,2)), 'ascend'); 
nAz = numel(azimuth);
nEl = numel(elevation);

% Prepare output folder & PDF
outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir, 'All_ROI_SpatialGrid_Flat.pdf');
if exist(pdfPath,'file'), delete(pdfPath); end

nROI = size(psth(1).alignedResponses, 1);

%% 2. Loop over each ROI
for roiIdx = 1:nROI
    fprintf('Plotting ROI %d/%d...\n', roiIdx, nROI);
    
    % --- Step A: Find Global Max for Normalization (across all positions) ---
    peaks = nan(nPos, 1);
    for k = 1:nPos
        tVec = psth(k).timeVector(:);
        data2d = squeeze(psth(k).alignedResponses(roiIdx, :, :)); % [Time x Trials]
        
        % Calculate F0 for normalization check
        preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        trialF0 = median(data2d(preMask, :), 1, 'omitnan');
        zeroed = data2d - trialF0; 
        
        meanTr = mean(zeroed, 2, 'omitnan');
        postMask = tVec >= peakTimeWindow(1) & tVec <= peakTimeWindow(2);
        if any(postMask)
            peaks(k) = max(meanTr(postMask), [], 'omitnan');
        else
            peaks(k) = 0;
        end
    end
    
    absMax = max(peaks);
    if isempty(absMax) || absMax <= 0, absMax = 1; end
    [~, bestIdx] = max(peaks);

    % --- Step B: Create Figure ---
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1400 950]);
    til = tiledlayout(nEl, nAz, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(sprintf('ROI %d Spatial Tuning (Per-Trial Median Subtracted)', roiIdx), 'FontWeight', 'bold');

    % --- Step C: Plot each spatial position ---
    for k = 1:nPos
        pos = psth(k).stimValue;
        col = find(azimuth == pos(1), 1);
        row = find(elevation == pos(2), 1);
        
        % Tile index: (row-1)*nAz + col
        idx = (row-1)*nAz + col;
        ax = nexttile(til, idx); hold(ax, 'on');
        
        tVec = psth(k).timeVector(:);
        data2d = squeeze(psth(k).alignedResponses(roiIdx, :, :));
        
        % PER-TRIAL MEDIAN SUBTRACTION (Ensures flat baseline)
        preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        trialF0 = median(data2d(preMask, :), 1, 'omitnan');
        zeroedTrials = data2d - trialF0; 
        
        % Calculate stats
        meanTr = mean(zeroedTrials, 2, 'omitnan');
        semTr  = std(zeroedTrials, 0, 2, 'omitnan') ./ sqrt(size(zeroedTrials, 2));
        
        % Normalize
        normMean = meanTr / absMax;
        normSEM  = semTr / absMax;
        
        % Plot window
        plotMask = tVec >= baselineWindow(1) & tVec <= timeWindow(2);
        t_plot = tVec(plotMask);
        m_plot = normMean(plotMask);
        s_plot = normSEM(plotMask);
        
        % Shaded SEM cloud 
        fill(ax, [t_plot; flipud(t_plot)], [m_plot + s_plot; flipud(m_plot - s_plot)], ...
            [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        
        % Mean PSTH line
        plot(ax, t_plot, m_plot, 'k', 'LineWidth', 1.5);
        
        % Reference markers
        xline(ax, 0, 'r--', 'LineWidth', 1); % Stim Onset
        yline(ax, 0, 'k:', 'Alpha', 0.4);    % True Baseline Anchor
        
        % Highlight the best position 
        if k == bestIdx
            set(ax, 'LineWidth', 2, 'Box', 'on', 'XColor', 'r', 'YColor', 'r');
        end
        
        % Labels and Formatting
        title(ax, sprintf('Az %d, El %.1f', pos(1), pos(2)), 'FontSize', 8);
        xlim(ax, [baselineWindow(1), timeWindow(2)]);
        ylim(ax, [-0.6, 1.2]); 
        set(ax, 'FontSize', 7);
        
        if col ~= 1, set(ax, 'YTickLabel', []); else ylabel(ax, 'norm \DeltaF/F'); end
        if row ~= nEl, set(ax, 'XTickLabel', []); else xlabel(ax, 'Time (s)'); end
    end

    % --- Step D: Save page ---
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');
    close(fig);
end

fprintf('Finished! Grid plots saved to: %s\n', pdfPath);
end