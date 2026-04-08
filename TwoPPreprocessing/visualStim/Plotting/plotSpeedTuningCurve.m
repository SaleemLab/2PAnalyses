function ax = plotSpeedTuningCurve(statMean, statStd, statCount, moveMean, moveStd, moveCount, speedBins, roiIdx, modulationValue, label, targetAx)
    
% plotSpeedTuningCurve(response.tuningCurve.statMean, response.tuningCurve.statStd, response.tuningCurve.statCount, response.tuningCurve.moveMean, response.tuningCurve.moveStd, response.tuningCurve.moveCount, response.tuningCurve.speedBins, 355, response.tuningCurve.modulationValue, 'dff', targetAx)
if nargin < 11 || isempty(targetAx)
        figure('Color', 'w');
        ax = gca;
else
        ax = targetAx;
end

    yStat = statMean(roiIdx);
    nStat = statCount(roiIdx);
    errStat = statStd(roiIdx) ./ sqrt(nStat);
    
    yMove = moveMean(roiIdx, :);
    nMove = moveCount(roiIdx, :);
    errMove = moveStd(roiIdx, :) ./ sqrt(nMove);
    
    xMove = speedBins(1:end-1) + diff(speedBins)/2; 
    xStat = 0;
    
    hold(ax, 'on');
    
    errorbar(ax, xMove, yMove, errMove, 'k-', 'LineWidth', 1.5, 'CapSize', 0);
    plot(ax, xMove, yMove, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    
    errorbar(ax, xStat, yStat, errStat, 'k', 'LineWidth', 1.5, 'CapSize', 0);
    plot(ax, xStat, yStat, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    title(ax, sprintf('%s\nROI %d: %.1f%% Mod', label, roiIdx, modulationValue));
    xlabel(ax, 'Speed (cm/s)');
    ylabel(ax, 'dFF (Neuropil Corr)');
    
    %defaultAxesProperties(gca, true);
    
    allY = [yStat, yMove];
    allErr = [errStat, errMove];
    ymin = min(allY - allErr);
    ymax = max(allY + allErr);
    ylim(ax, [min(0, ymin*1.1), ymax*1.2]);
end

%% 
