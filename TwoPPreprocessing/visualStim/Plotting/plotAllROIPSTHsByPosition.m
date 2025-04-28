function plotAllROIPSTHsByPosition(sessionFileInfo, response, timeWindow, baselineWindow, peakTimeWindow)
% plotAllROIPSTHsByPosition  Loop over every ROI and save one page per ROI
%   into a single multi‐page PDF. Each page tiles that ROI’s normalized,
%   baseline‑subtracted PSTHs across spatial positions, with the best
%   post‑stimulus peak highlighted in red.
%
% Usage:
%   plotAllROIPSTHsByPosition(sessionFileInfo, response)
%   plotAllROIPSTHsByPosition(sessionFileInfo, response, [0 2.8], [-0.25 0])
%
% Inputs:
%   sessionFileInfo - struct with .Directories.save_folder
%   response        - struct with .psthData (from getTrialResponsePSTHsV4)
%   timeWindow      - [t0 t1] post‑stimulus window in s (default [0 3])
%   baselineWindow  - [t0 t1] pre‑stimulus window in s (default [−0.4 0])

% Defaults
if nargin<3 || isempty(timeWindow)
    timeWindow = [0 3];
end
if nargin<4 || isempty(baselineWindow)
    baselineWindow = [-0.4 0];
end
if nargin<4 || isempty(baselineWindow)
    peakTimeWindow = [0 1];
end

psth   = response.psthData;
nPos   = numel(psth);
stimVs = vertcat(psth.stimValue);
azimuth   = unique(stimVs(:,1));  % columns
elevation = unique(stimVs(:,2));  % rows
nAz = numel(azimuth);
nEl = numel(elevation);

% Prepare output folder & PDF
outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir,'All_ROI_PSTH_byPosition.pdf');
if exist(pdfPath,'file'), delete(pdfPath); end

% Number of ROIs
nROI = size( psth(1).alignedResponses, 1 );

% Loop over each ROI
for roiIdx = 1:nROI
    % --- Find best position for this ROI ---
    peaks = nan(nPos,1);
    for k = 1:nPos
        tVec   = psth(k).timeVector(:)';
        AR3d   = squeeze(psth(k).alignedResponses(roiIdx,:,:));
        data2d = reshape(AR3d, numel(tVec), []);
        medianTr = median(data2d,2,'omitnan');
        preMask = tVec>=baselineWindow(1)&tVec<=baselineWindow(2);
        F0 = any(preMask)*median(medianTr(preMask),'omitnan');
        zeroed = medianTr - F0;
%         postMask = tVec>=timeWindow(1)&tVec<=timeWindow(2);
        postMask = tVec>=peakTimeWindow(1)&tVec<=peakTimeWindow(2);
        if any(postMask)
            peaks(k) = max(zeroed(postMask),[], 'omitnan');
        else
            peaks(k) = -Inf;
        end
    end
    [~, bestIdx] = max(peaks);
    absMax = max(peaks);
    if isempty(absMax)||absMax==0, absMax=1; end

    % --- Create invisible figure for this ROI ---
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 800]);
    til = tiledlayout(nEl, nAz, 'TileSpacing','compact','Padding','compact');
    sgtitle(sprintf('ROI %d PSTHs', roiIdx), 'FontWeight','bold');

    % --- Plot each spatial position ---
    for k = 1:nPos
        pos = psth(k).stimValue;
        col = find(azimuth==pos(1),1);
        row = find(elevation==pos(2),1);
        idx = (row-1)*nAz + col;
        ax = nexttile(til, idx); hold(ax,'on');

        % extract & baseline‐subtract
        tVec   = psth(k).timeVector(:)';
        AR3d   = squeeze(psth(k).alignedResponses(roiIdx,:,:));
        data2d = reshape(AR3d, numel(tVec), []);
        medianTr = median(data2d,2,'omitnan');
        semTr  = std(data2d,0,2,'omitnan') ./ sqrt(size(data2d,2));
        preMask = tVec>=baselineWindow(1)&tVec<=baselineWindow(2);
        if any(preMask)
            F0 = median(medianTr(preMask),'omitnan');
        else
            F0 = 0;
        end
        zeroed = medianTr - F0;

        % normalize
        normMedian = zeroed / absMax;
        normSEM  = semTr   / absMax;

        % mask for plotting
        maskPlot = (tVec>=baselineWindow(1)&tVec<=0) | ...
                   (tVec>=timeWindow(1)&tVec<=timeWindow(2));
        t_plot    = tVec(maskPlot);
        m_plot    = normMedian(maskPlot);
        s_plot    = normSEM(maskPlot);

        % shaded SEM
        xPatch = [t_plot,         fliplr(t_plot)];
        yPatch = [m_plot'+s_plot', fliplr(m_plot'-s_plot')];
        fill(ax, xPatch, yPatch, [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.5);

        % median line
        plot(ax, t_plot, m_plot, 'k', 'LineWidth',1.5);

        % mark lines
        xline(ax, 0,    'r--','LineWidth',1);
        for v=[0.25,0.5,0.75]
            if v>=baselineWindow(1)&&v<=timeWindow(2)
                xline(ax, v,':','Color',[0.5 0.5 0.5],'LineWidth',1);
            end
        end

        % highlight best
        if k==bestIdx
            set(ax,'LineWidth',2,'Box','on','XColor','r','YColor','r');
        end

        % formatting
        title(ax, sprintf('Az %d°, El %d°', pos(1), pos(2)), 'FontSize',9);
        xlabel(ax,'Time (s)'); ylabel(ax,'norm ΔF/F');
        xlim(ax,[baselineWindow(1), timeWindow(2)]);
        ylim(ax,[-1,1]);
        box(ax,'off'); set(ax,'FontSize',8);
    end

    % --- Append this ROI’s page ---
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType','vector');
    close(fig);
end
end
