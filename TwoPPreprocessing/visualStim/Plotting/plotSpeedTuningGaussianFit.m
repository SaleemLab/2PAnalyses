function plotSpeedTuningGaussianFit(speeds, meanVals, semVals, gaussParams, char, R2, pval, fitR2, preStimBaselineMean, ax)
% plotSpeedTuningGaussianFit(speeds, meanVals, semVals, gaussParams, char, R2, pval, fitR2, preStimBaselineMean, ax)
%
% R2/pval  = CROSS-VALIDATED values (e.g. unit.runR2, unit.runR2_pval) --
%            the reliability criterion, shown with its p-value.
% fitR2    = DESCRIPTIVE Gaussian fit R2 (unit.gaussR2_run) -- shown
%            alongside for reference/write-up later, no p-value exists
%            for this one.

if nargin < 10 || isempty(ax)
    ax = gca;
end

categoryNames  = {'Lowpass', 'Highpass', 'Bandpass', 'Trough'};
categoryColors = { ...
    hex2rgb('#3288bd'), ...  % lowpass  - dark blue
    hex2rgb('#a6cee3'), ...  %  highpass - light blue
    hex2rgb('#e6ab02'), ...  % bandpass - orange/gold
    hex2rgb('#66a61e'), ...  % trough   - green
    };

thisColor = categoryColors{char};
thisName  = categoryNames{char};

gaussFun = @(params, xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));

nSpeeds = numel(speeds);

% Ed's exact method: log2-transform speed values; force 0 speed to sit
% one step below the slowest real speed (log2(0) is undefined)
speeds = speeds(:)';
xVals = nan(1, nSpeeds);
nonZeroMask = speeds > 0;
xVals(nonZeroMask) = log2(speeds(nonZeroMask));
if any(~nonZeroMask)
    xVals(~nonZeroMask) = min(xVals(nonZeroMask)) - 1;
end

xFine = linspace(min(xVals), max(xVals), 200);
fitIdxFine = interp1(xVals, 1:nSpeeds, xFine, 'linear', 'extrap');
yFine = feval(gaussFun, gaussParams, fitIdxFine);

axes(ax); 
hold(ax, 'on');

if ~isempty(preStimBaselineMean) && ~isnan(preStimBaselineMean)
    yline(ax, preStimBaselineMean, 'k--', 'LineWidth', 1);
end

plot(ax, xFine, yFine, '-', 'Color', thisColor, 'LineWidth', 2.5);

errorbar(ax, xVals, meanVals, semVals, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 5, 'LineWidth', 1, 'CapSize', 0);

xlabel(ax, 'Visual speed (\circ/s)');
ylabel(ax, '\DeltaF/F');
title(ax, thisName, 'FontWeight', 'normal');

set(ax, 'XTick', xVals, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false), ...
    'XTickLabelRotation', 0);

xl = xlim(ax); yl = ylim(ax);
if pval < 0.001
    pvalStr = 'p < 0.001';
else
    pvalStr = sprintf('p = %.3f', pval);
end
text(ax, xl(2), yl(2), sprintf('CV R^2 = %.3f, %s\nFit R^2 = %.3f', R2, pvalStr, fitR2), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontAngle', 'italic', 'FontSize', 9);

box(ax, 'off');
set(ax, 'TickDir', 'out');
defaultAxesProperties(ax, true);

end

%% local helper
function rgb = hex2rgb(hexStr)
    hexStr = strrep(hexStr, '#', '');
    rgb = [hex2dec(hexStr(1:2)), hex2dec(hexStr(3:4)), hex2dec(hexStr(5:6))] / 255;
end
