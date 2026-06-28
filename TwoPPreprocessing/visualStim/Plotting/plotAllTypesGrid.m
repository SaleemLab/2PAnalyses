function plotAllTypesGrid(response, targetStruct, useField, nPerType)
% PLOTALLTYPESGRID
% Plots a random sample of each tuning type in a grid for visual inspection.
%
% USAGE:
%   plotAllTypesGrid(response, 'tuningCurve', 'spks', 9)

    if nargin < 4, nPerType = 9; end

    cls      = response.(targetStruct).(useField).classification;
    pval     = response.(targetStruct).(useField).pValFull;
    y        = response.(targetStruct).(useField).moveMean;
    edges    = response.(targetStruct).speedBins;
    movingCenters = (edges(1:end-1) + diff(edges)/2)';

    r2_thresh = 0.1;
    p_thresh  = 0.01;

    types      = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};
    typeLabels = {'Low-Pass', 'High-Pass', 'Band-Pass', 'Trough-Inverted'};
    colors     = {[0.85 0.33 0.1], [0 0.45 0.74], [0 0.45 0.74], [0.85 0.33 0.1]};

    for t = 1:numel(types)
        % Get passing cells of this type
        idx = find(strcmp(cls.tuningType, types{t}) & ...
                   cls.R2 >= r2_thresh & pval <= p_thresh);

        if isempty(idx)
            fprintf('No %s cells passing thresholds\n', types{t});
            continue;
        end

        % Random sample
        nPlot   = min(nPerType, numel(idx));
        plotIdx = idx(randperm(numel(idx), nPlot));
        nCols   = 3;
        nRows   = ceil(nPlot / nCols);

        figure('Name', sprintf('%s — %s', typeLabels{t}, useField), ...
               'Color', 'w', 'Position', [50 50 900 300*nRows]);

        for k = 1:nPlot
            r = plotIdx(k);
            subplot(nRows, nCols, k);
            plot(movingCenters, y(r,:), 'ko-', ...
                'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerSize', 4);
            hold on;

            % Plot fit if params are valid
            p = cls.fitParams(r,:);
            if all(isfinite(p))
                gaussFun = @(params, xdata) params(1) + params(2) .* ...
                    exp(-(((xdata - params(3)).^2) ./ (2 * (params(4).^2))));
                xDense = linspace(min(movingCenters), max(movingCenters), 200);
                yFit   = gaussFun(p, xDense);
                plot(xDense, yFit, '-', 'Color', colors{t}, 'LineWidth', 2);
            end

            title(sprintf('ROI %d | R2=%.2f | p=%.3f', r, cls.R2(r), pval(r)), ...
                'FontSize', 7);
            xlabel('Speed (cm/s)', 'FontSize', 7);
            box off; set(gca, 'TickDir', 'out', 'FontSize', 7);
        end

        sgtitle(sprintf('%s — %s (n=%d shown of %d passing)', ...
            typeLabels{t}, useField, nPlot, numel(idx)), ...
            'FontWeight', 'bold', 'FontSize', 10);
    end
end
