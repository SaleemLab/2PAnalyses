function plotPopulationPSTHs(response, baselineWindow)
% plotPopulationPSTHs  Plot mean±SD PSTHs across all ROIs×trials per position
%
% Usage:
%   plotPopulationPSTHs(response)
%   plotPopulationPSTHs(response, [-0.5 0])
%
% Inputs:
%   response       - struct with .psthData (as produced by V4)
%   baselineWindow - optional [t0 t1] in seconds for prestim shading (default: [-0.5 0])

if nargin<2, baselineWindow = [-0.5 0]; end

psth = response.psthData;
nPos = numel(psth);

% layout as close to square as possible
nCols = ceil(sqrt(nPos));
nRows = ceil(nPos/nCols);

figure('Color','w','Position',[100 100 900 700]);
t = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
sgtitle('Population PSTH (mean \pm SD) by Position','FontWeight','bold');

for k = 1:nPos
    % extract alignedResponses: [nROI x nTime x nTrials]
    AR = psth(k).alignedResponses;
    tVec = psth(k).timeVector;              % 1 x nTime
    [nROI, nTime, nTrials] = size(AR);
    
    % flatten ROI×trials → rows
    data2d = reshape(AR, nROI*nTrials, nTime);
    
    % compute mean & SD
    popMean = nanmean(data2d, 1);           % 1 x nTime
    popSD   = nanstd(data2d, 0, 1);         % 1 x nTime
    
    ax = nexttile(t, k);
    hold(ax,'on');
    
    % shade mean ± SD
    fill([tVec fliplr(tVec)], ...
         [popMean+popSD, fliplr(popMean-popSD)], ...
         [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.5, 'Parent',ax);
    
    % plot mean trace
    plot(ax, tVec, popMean, 'k-', 'LineWidth',1.5);
    
    % mark stimulus onset
    xline(ax, 0, 'r--', 'LineWidth',1);
    
    % optional baseline shading
    bw = baselineWindow;
    yl = get(ax,'YLim');
    patch([bw(1) bw(2) bw(2) bw(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.3, 'Parent',ax);
    
    % labels & title
    pos = psth(k).stimValue;
    title(ax, sprintf('Az %d°, El %d°', pos(1), pos(2)), 'FontSize',10);
    xlabel(ax,'Time (s)');
    ylabel(ax,'\DeltaF/F');
    box(ax,'off');
    set(ax,'FontSize',9);
    ylim(ax, yl);  % keep Y‐limits consistent
end

end
