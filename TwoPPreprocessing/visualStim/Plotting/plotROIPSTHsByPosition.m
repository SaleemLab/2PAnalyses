function plotROIPSTHsByPosition(sessionFileInfo, response, roiIdx, timeWindow, baselineWindow)
% plotROIPSTHsByPosition  Plot one ROI’s PSTH across spatial positions,
% normalized and baseline‑subtracted, with the strongest post‑stimulus peak
% highlighted in red. 
%
% Usage:
%   plotROIPSTHsByPosition(sessionFileInfo,response, roiIdx, timeWindow, baselineWindow)
%   plotROIPSTHsByPosition(sessionFileInfo,response, roiIdx)
%
% Inputs:
%   response       - struct with field .psthData (from getTrialResponsePSTHsV4)
%   roiIdx         - index of the ROI to plot
%   timeWindow     - [t0 t1] post‑stimulus window in seconds (default [0 1.5])
%   baselineWindow - [t0 t1] pre‑stimulus window in seconds (default [-0.4 0])

% --- Defaults
if nargin<3 || isempty(timeWindow)
    timeWindow = [0 1.5]; % changed from 3 to 1.5
end
if nargin<4 || isempty(baselineWindow)
    baselineWindow = [-0.4 0];
end

psth   = response.psthData;
nPos   = numel(psth);
stimVs = vertcat(psth.stimValue);
azimuth   = unique(stimVs(:,1));  % columns
elevation = unique(stimVs(:,2));  % rows
nAz = numel(azimuth);
nEl = numel(elevation);

% --- Find which position has the largest positive peak post‑stimulus
peaks = nan(nPos,1);
for k = 1:nPos
    tVec   = psth(k).timeVector(:)';
    AR3d   = squeeze(psth(k).alignedResponses(roiIdx,:,:));     % [1 x T x trials]
    data2d = reshape(AR3d, numel(tVec), []);                    % [T x trials]
    meanTr = mean(data2d, 2, 'omitnan');                        % [T x 1]
    % baseline subtract
    preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
    if any(preMask)
        F0 = mean(meanTr(preMask), 'omitnan');
    else
        F0 = 0;
    end
    zeroed = meanTr - F0;
    % post‑stim mask
    postMask = tVec >= timeWindow(1) & tVec <= timeWindow(2);
    if any(postMask)
        peaks(k) = max(zeroed(postMask), [], 'omitnan');
    else
        peaks(k) = -Inf;
    end
end
[~, bestIdx] = max(peaks);

% --- Normalize to the global maximum across post‑stim peaks
absMax = max(peaks);
if isempty(absMax) || absMax == 0
    absMax = 1;
end

% --- Create figure
fig = figure('Color','w','Position',[100 100 1200 800]);
til = tiledlayout(nEl, nAz, 'TileSpacing','compact', 'Padding','compact');
sgtitle(sprintf('ROI %d PSTHs (normalized, best pos highlighted)', roiIdx), ...
        'FontWeight','bold');

% --- Plot each position in its spatial tile
for k = 1:nPos
    pos = psth(k).stimValue;
    col = find(azimuth   == pos(1), 1);
    row = find(elevation == pos(2), 1);
    idx = (row-1)*nAz + col;
    ax  = nexttile(til, idx);
    hold(ax,'on');

    % Extract and process trace
    tVec   = psth(k).timeVector(:)';
    AR3d   = squeeze(psth(k).alignedResponses(roiIdx,:,:));
    data2d = reshape(AR3d, numel(tVec), []);
    meanTr = mean(data2d,2,'omitnan');
    semTr  = std(data2d,0,2,'omitnan') ./ sqrt(size(data2d,2));

    % Baseline subtract
    preMask = tVec >= baselineWindow(1) & tVec <= baselineWindow(2);
    if any(preMask)
        F0 = mean(meanTr(preMask), 'omitnan');
    else
        F0 = 0;
    end
    zeroed = meanTr - F0;

    % Normalise
    normMean = zeroed / absMax;
    normSEM  = semTr   / absMax;

    % Mask for plotting: baselineWindow then full post‑stim
    maskPlot = (tVec >= baselineWindow(1) & tVec <= 0) | ...
               (tVec >= timeWindow(1)     & tVec <= timeWindow(2));
    t_plot    = tVec(maskPlot);
    m_plot    = normMean(maskPlot);
    s_plot    = normSEM(maskPlot);

    % Shaded SEM
    xPatch = [t_plot, fliplr(t_plot)];
    yPatch = [m_plot'+s_plot', fliplr(m_plot'-s_plot')];
    fill(ax, xPatch, yPatch, [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.5);

    % Mean line
    plot(ax, t_plot, m_plot, 'k', 'LineWidth',1.5);

    % Zero and mark lines
    xline(ax, 0,    'r--', 'LineWidth',1);
    for v = [0.25, 0.5, 0.75]
        if v >= baselineWindow(1) && v <= timeWindow(2)
            xline(ax, v, ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1);
        end
    end

    % Highlight best position
    if k == bestIdx
        set(ax, 'LineWidth',2, 'Box','on', 'XColor','r', 'YColor','r');
    end

    % Labels and formatting
    title(ax, sprintf('Az %d°, El %d°', pos(1), pos(2)), 'FontSize',9);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'norm ΔF/F');
    xlim(ax, [baselineWindow(1), timeWindow(2)]);
    ylim(ax, [-1 1]);
    box(ax,'off'); set(ax,'FontSize',8);
end

% --- Save figure to session folder ---
outDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'ROI_PSTH_byPosition');
if ~exist(outDir,'dir'), mkdir(outDir); end
outPath = fullfile(outDir, sprintf('ROI%d_normPSTH_highlighted.pdf', roiIdx));
exportgraphics(fig, outPath, 'ContentType','vector');

end
