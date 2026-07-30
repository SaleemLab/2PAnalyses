function plotSpeedTuningGaussianFit_dual(speeds, meanVals_stat, meanVals_run, ...
    semVals_stat, semVals_run, gaussParams_stat, gaussParams_run, ...
    char_stat, char_run, R2_stat, R2_run, pval_stat, pval_run, ...
    fitR2_stat, fitR2_run, preStimBaselineMean_stat, preStimBaselineMean_run, ax)
% plotSpeedTuningGaussianFit_dual(...)
%
% Same style as plotSpeedTuningGaussianFit, but overlays Stationary
% (black) and Locomotion (red) fits on one axis.
%
% R2/pval  = CROSS-VALIDATED values (e.g. unit.statR2/runR2, ..._pval) --
%            the reliability criterion, shown with its p-value.
% fitR2    = DESCRIPTIVE Gaussian fit R2 (unit.gaussR2_stat/run) -- shown
%            alongside for reference/write-up later, no p-value exists
%            for this one.
% preStimBaselineMean_stat/_run = per-condition pre-stim baseline
%            (unit.preStimBaselineMean_stat / _run), each drawn as a
%            dashed line in its condition's color.

if nargin < 18 || isempty(ax)
    ax = gca;
end

categoryNames = {'Lowpass', 'Highpass', 'Bandpass', 'Trough'};
name_stat = categoryNames{char_stat};
name_run  = categoryNames{char_run};

color_stat = [0 0 0];   % black -- Stationary
color_run  = [1 0 0];   % red   -- Locomotion

gaussFun = @(params, xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));

nSpeeds = numel(speeds);
speeds = speeds(:)';
xVals = nan(1, nSpeeds);
nonZeroMask = speeds > 0;
xVals(nonZeroMask) = log2(speeds(nonZeroMask));
if any(~nonZeroMask)
    xVals(~nonZeroMask) = min(xVals(nonZeroMask)) - 1;
end

xFine = linspace(min(xVals), max(xVals), 200);
fitIdxFine = interp1(xVals, 1:nSpeeds, xFine, 'linear', 'extrap');
yFine_stat = feval(gaussFun, gaussParams_stat, fitIdxFine);
yFine_run  = feval(gaussFun, gaussParams_run, fitIdxFine);

axes(ax);
hold(ax, 'on');

% --- pre-stim baselines, one dashed line per condition, condition-colored ---
if ~isempty(preStimBaselineMean_stat) && ~isnan(preStimBaselineMean_stat)
    yline(ax, preStimBaselineMean_stat, '--', 'Color', color_stat, 'LineWidth', 1);
end
if ~isempty(preStimBaselineMean_run) && ~isnan(preStimBaselineMean_run)
    yline(ax, preStimBaselineMean_run, '--', 'Color', color_run, 'LineWidth', 1);
end

% Stationary
plot(ax, xFine, yFine_stat, '-', 'Color', color_stat, 'LineWidth', 2.5);
errorbar(ax, xVals, meanVals_stat, semVals_stat, 'o', ...
    'Color', color_stat, 'MarkerFaceColor', color_stat, ...
    'MarkerSize', 5, 'LineWidth', 1, 'CapSize', 0);

% Locomotion
plot(ax, xFine, yFine_run, '-', 'Color', color_run, 'LineWidth', 2.5);
errorbar(ax, xVals, meanVals_run, semVals_run, 'o', ...
    'Color', color_run, 'MarkerFaceColor', color_run, ...
    'MarkerSize', 5, 'LineWidth', 1, 'CapSize', 0);

xlabel(ax, 'Visual speed (\circ/s)');
ylabel(ax, '\DeltaF/F');
title(ax, sprintf('Stat: %s | Run: %s', name_stat, name_run), 'FontWeight', 'normal');
set(ax, 'XTick', xVals, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false), ...
    'XTickLabelRotation', 0);

xl = xlim(ax); yl = ylim(ax);

pvalStr_stat = tern(pval_stat < 0.001, 'p < 0.001', sprintf('p = %.3f', pval_stat));
pvalStr_run  = tern(pval_run  < 0.001, 'p < 0.001', sprintf('p = %.3f', pval_run));

text(ax, xl(2), yl(2), sprintf('Stat: CV R^2=%.3f, %s, Fit R^2=%.3f', R2_stat, pvalStr_stat, fitR2_stat), ...
    'Color', color_stat, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontAngle', 'italic', 'FontSize', 8);
text(ax, xl(2), yl(2) - 0.08*range(yl), sprintf('Run: CV R^2=%.3f, %s, Fit R^2=%.3f', R2_run, pvalStr_run, fitR2_run), ...
    'Color', color_run, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontAngle', 'italic', 'FontSize', 8);

box(ax, 'off');
set(ax, 'TickDir', 'out');
defaultAxesProperties(ax, true);

end

%% local helper
function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end