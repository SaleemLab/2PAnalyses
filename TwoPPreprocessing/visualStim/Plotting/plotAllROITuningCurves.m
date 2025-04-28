function plotAllROITuningCurves(sessionFileInfo, response, timeWindow, baselineWindow)
% plotAllROITuningCurves  Loop over every ROI and save one page per ROI
%   into a single multi-page PDF. Each page tiles that ROI’s normalized,
%   baseline-subtracted PSTHs across stimulus conditions.
%
% Usage:
%   plotAllROITuningCurves(sessionFileInfo, response)
%   plotAllROITuningCurves(sessionFileInfo, response, [0 2.8], [-0.25 0])
%
% Inputs:
%   sessionFileInfo - struct with .Directories.save_folder
%   response        - struct with .psthData (from getTrialResponsePSTHsV4)
%   timeWindow      - [t0 t1] post-stimulus window in s (default [0 3])
%   baselineWindow  - [t0 t1] pre-stimulus window in s (default [−0.4 0])

% Defaults
if nargin<3 || isempty(timeWindow)
    timeWindow = [0 3];
end
if nargin<4 || isempty(baselineWindow)
    baselineWindow = [-0.4 0];
end

psth   = response.psthData;
nStim  = numel(psth);
stimVs = vertcat(psth.stimValue);

% Output folder & PDF
outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir,'All_ROI_TuningCurves.pdf');
if exist(pdfPath,'file'), delete(pdfPath); end

% Number of ROIs
nROI = size(psth(1).alignedResponses, 1);

% Loop over ROIs
for roiIdx = 1:nROI
    disp(['Plotting ROIIdx: ', num2str(roiIdx)])

    % Determine best response
    peaks = nan(nStim,1);
    for k = 1:nStim
        tVec   = psth(k).timeVector(:)';
        AR3d   = squeeze(psth(k).alignedResponses(roiIdx,:,:));
        data2d = reshape(AR3d, numel(tVec), []);
        meanTr = mean(data2d,2,'omitnan');
        preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        F0 = any(preMask) * mean(meanTr(preMask),'omitnan');
        zeroed = meanTr - F0;
        postMask = tVec >= timeWindow(1) & tVec <= timeWindow(2);
        peaks(k) = max(zeroed(postMask),[],'omitnan');
    end
    [~, bestIdx] = max(peaks);
    absMax = max(peaks);
    if isempty(absMax) || absMax == 0, absMax = 1; end

    % Create figure
    nCols = ceil(sqrt(nStim));
    nRows = ceil(nStim / nCols);
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 800]);
    til = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
    sgtitle(sprintf('ROI %d PSTHs', roiIdx), 'FontWeight','bold');

    for k = 1:nStim
        ax = nexttile(til); hold(ax, 'on');

        % Extract and baseline-subtract
        tVec   = psth(k).timeVector(:)';
        AR3d   = squeeze(psth(k).alignedResponses(roiIdx,:,:));
        data2d = reshape(AR3d, numel(tVec), []);
        meanTr = mean(data2d,2,'omitnan');
        semTr  = std(data2d,0,2,'omitnan') ./ sqrt(size(data2d,2));
        preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
        F0 = mean(meanTr(preMask),'omitnan');
        zeroed = meanTr - F0;

        % Normalize
        normMean = zeroed / absMax;
        normSEM  = semTr / absMax;

        % Plot window
        maskPlot = (tVec>=baselineWindow(1)&tVec<=0) | ...
                   (tVec>=timeWindow(1)&tVec<=timeWindow(2));
        t_plot = tVec(maskPlot);
        m_plot = normMean(maskPlot);
        s_plot = normSEM(maskPlot);

        % Plot shaded error
        fill(ax, [t_plot fliplr(t_plot)], ...
                  [m_plot'+s_plot' fliplr(m_plot'-s_plot')], ...
                  [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.5);
        plot(ax, t_plot, m_plot, 'k', 'LineWidth', 1.5);
        xline(ax, 0, 'r--', 'LineWidth', 1);

        for v = 1
            if v >= baselineWindow(1) && v <= timeWindow(2)
                xline(ax, v, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            end
        end

        % Highlight best
        if k == bestIdx
            set(ax, 'LineWidth', 2, 'Box', 'on', 'XColor', 'r', 'YColor', 'r');
        end

        stimLabel = mat2str(psth(k).stimValue); % Convert stimValue to string
        title(ax, stimLabel, 'Interpreter', 'none', 'FontSize', 9);
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'norm ΔF/F');
        xlim(ax, [baselineWindow(1), timeWindow(2)]);
        ylim(ax, [-1, 1]);
        box(ax, 'off'); set(ax, 'FontSize', 8);
    end

    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');
    close(fig);
end
end
