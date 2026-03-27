function plotAllROITuningCurves(sessionFileInfo, response, timeWindow, baselineWindow)
% plotAllROITuningCurves  Loop over every ROI and save one page per ROI
%   into a single multi-page PDF with per-trial baseline subtraction.
%
% Usage:
%   plotAllROITuningCurves(sessionFileInfo, response)
%   plotAllROITuningCurves(sessionFileInfo, response, [0 3], [-0.5 0])

%% 1. Defaults and Setup
if nargin < 3 || isempty(timeWindow)
    timeWindow = [0 3];
end
if nargin < 4 || isempty(baselineWindow)
    baselineWindow = [-0.5 0];
end

psth   = response.psthData;
nStim  = numel(psth);
stimVs = vertcat(psth.stimValue);

% Prepare output folder & PDF
outDir  = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir, 'All_ROI_TuningCurves.pdf');
if exist(pdfPath, 'file'), delete(pdfPath); end

% Number of ROIs
nROI = size(psth(1).alignedResponses, 1);

%% 2. Loop over each ROI
for roiIdx = 1:nROI
    fprintf('Plotting ROI %d/%d...\n', roiIdx, nROI);
    
    % --- Step A: Calculate Best Response (Peak) for Normalization ---
    peaks = nan(nStim, 1);
    for k = 1:nStim
        tVec = psth(k).timeVector(:);
        % AR3d is [ROI x Time x Trials]
        data2d = squeeze(psth(k).alignedResponses(roiIdx, :, :)); % [Time x Trials]
        
        % Per-trial baseline subtraction
        preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        trialF0 = mean(data2d(preMask, :), 1, 'omitnan'); % Mean per trial
        zeroedTrials = data2d - trialF0; % Subtract F0 from each trial
        
        % Calculate Mean Trace for peak detection
        meanTr = mean(zeroedTrials, 2, 'omitnan');
        postMask = tVec >= 0.5 & tVec <= 2.0; % Using your specified window
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
    nCols = ceil(sqrt(nStim));
    nRows = ceil(nStim / nCols);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 800]);
    til = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    sgtitle(sprintf('ROI %d PSTHs (Baseline Subtracted)', roiIdx), 'FontWeight', 'bold');

    % --- Step C: Plot each condition ---
    for k = 1:nStim
        ax = nexttile(til); hold(ax, 'on');
        
        tVec = psth(k).timeVector(:);
        data2d = squeeze(psth(k).alignedResponses(roiIdx, :, :));
        
        % PER-TRIAL SUBTRACTION (The fix for the flat baseline)
        preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        trialF0 = mean(data2d(preMask, :), 1, 'omitnan');
        zeroedTrials = data2d - trialF0; 
        
        % Calculate mean and SEM from zeroed trials
        meanTr = mean(zeroedTrials, 2, 'omitnan');
        semTr  = std(zeroedTrials, 0, 2, 'omitnan') ./ sqrt(size(zeroedTrials, 2));
        
        % Normalize by the best condition's peak
        normMean = meanTr / absMax;
        normSEM  = semTr / absMax;
        
        % Plotting
        maskPlot = tVec >= baselineWindow(1) & tVec <= timeWindow(2);
        t_plot = tVec(maskPlot);
        m_plot = normMean(maskPlot);
        s_plot = normSEM(maskPlot);
        
        % Shaded SEM
        fill(ax, [t_plot; flipud(t_plot)], [m_plot + s_plot; flipud(m_plot - s_plot)], ...
            [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        
        % Mean Trace
        plot(ax, t_plot, m_plot, 'k', 'LineWidth', 1.5);
        
        % Reference Lines
        xline(ax, 0, 'r--', 'LineWidth', 1); % Stimulus Onset
        yline(ax, 0, 'k:', 'Alpha', 0.3);    % Zero Baseline
        
        % Highlight the best condition
        if k == bestIdx
            set(ax, 'LineWidth', 2, 'Box', 'on', 'XColor', 'r', 'YColor', 'r');
        end
        
        % Labels
        title(ax, mat2str(psth(k).stimValue), 'FontSize', 8);
        xlim(ax, [baselineWindow(1), timeWindow(2)]);
        ylim(ax, [-0.5, 1.2]); % Adjusted for normalized view
        set(ax, 'FontSize', 8);
        if mod(k-1, nCols) ~= 0, ylabel(ax, ''); else ylabel(ax, 'norm \DeltaF/F'); end
        if k <= (nRows-1)*nCols, xlabel(ax, ''); else xlabel(ax, 'Time (s)'); end
    end

    % --- Step D: Append to PDF ---
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');
    close(fig);
end

fprintf('Finished! Saved to: %s\n', pdfPath);

end