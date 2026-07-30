function fig = plotPooledSMI_CDF(pooledData, targetArea, compareBy, daysToPlot)
% plotPooledSMI_CDF: Plots CDF of SMI values for a given area.
% compareBy: 'Total' (one curve) or 'Days' (separate curve per day).
% daysToPlot (optional, only used when compareBy = 'Days'): vector of day
%            numbers to include, e.g. [1 3 5]. Defaults to all days found
%            in pooledData.(targetArea).Days.

if ~isfield(pooledData, targetArea)
    error('Area %s not found in pooled data.', targetArea);
end
if nargin < 4, daysToPlot = []; end % empty = plot all days

% Fixed palette so a given day always gets the same color, regardless of
% how many/which days end up being plotted (unlike lines(), which just
% assigns colors sequentially based on however many curves are present).
dayColorPalette = {'r', 'm', [0.4 0.7 0.2], 'k', 'b'};

fig = figure('Color', 'w', 'Position', [200 200 500 450]);
hold on;

if strcmpi(compareBy, 'Total')
    % Plot single curve for the entire area
    [f, x] = ecdf(pooledData.(targetArea).AllSMI);
    plot(x, f, 'k', 'LineWidth', 2.5);
    legendLabels = {sprintf('%s (n=%d)', targetArea, length(pooledData.(targetArea).AllSMI))};

elseif strcmpi(compareBy, 'Days')
    % Plot separate curve for each Day of Experience
    dayFields = fieldnames(pooledData.(targetArea).Days);
    legendLabels = {};
    for d = 1:length(dayFields)
        % Extract the actual day number from the field name (e.g. 'Day3' -> 3)
        % so the color is tied to the day identity, not plot order.
        dayNum = str2double(regexp(dayFields{d}, '\d+', 'match', 'once'));
        if isnan(dayNum), dayNum = d; end % fallback if field name has no digits

        % Skip this day entirely if the user asked for a specific subset
        % and this day isn't in it.
        if ~isempty(daysToPlot) && ~ismember(dayNum, daysToPlot)
            continue;
        end

        daySMI = pooledData.(targetArea).Days.(dayFields{d});
        if isempty(daySMI), continue; end

        thisColor = dayColorPalette{mod(dayNum - 1, length(dayColorPalette)) + 1};

        [f, x] = ecdf(daySMI);
        plot(x, f, 'Color', thisColor, 'LineWidth', 2);
        legendLabels{end+1} = sprintf('%s: %s (n=%d)', targetArea, dayFields{d}, length(daySMI));
    end
end

% Aesthetics matching Diamanti et al.
line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Zero line
line([-1 1], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); % Median guide
grid on;
xlim([-1 1]); ylim([0 1]);
xlabel('Spatial Modulation Index (SMI)');
ylabel('Cumulative Probability');
title(sprintf('Spatial Modulation in %s', targetArea));
legend(legendLabels, 'Location', 'southeast', 'Box', 'off');
end