function figHandle = plotModulation_CDF(allData, varargin)
% plotModulation_CDF Plots superimposed CDFs of Modulation Index for specified days.

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addParameter(p, 'DaysToPlot', [1 3 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @ischar);
    addParameter(p, 'ExcludeMice', {'M24043', 'M24046', 'M24048', 'M24049'}, @iscellstr);
    addParameter(p, 'SavePath', '', @ischar);
    parse(p, allData, varargin{:});
    
    daysToPlot = p.Results.DaysToPlot;
    typeToPlot = p.Results.TypeToPlot;
    excludeMice = p.Results.ExcludeMice;

    uMice = setdiff(unique({allData.MouseID}), excludeMice, 'stable');
    nMice = length(uMice);
    colors = lines(max(daysToPlot)); % Consistent colors for days

    figHandle = figure('Position', [100 100 350*nMice 300]);
    t = tiledlayout(figHandle, 1, nMice, 'TileSpacing', 'compact', 'Padding', 'compact');
    axList = gobjects(0);

    for m = 1:nMice
        ax = nexttile; axList(end+1) = ax;
        hold(ax, 'on');
        hasData = false;

        for d = daysToPlot
            % Find data
            session = allData(strcmp({allData.MouseID}, uMice{m}) & ...
                              [allData.Day] == d & ...
                              strcmpi({allData.Type}, typeToPlot));
            if isempty(session), continue; end
            session = session(1);
            
            % --- PLOT CDF ---
            mi = session.Modulation(~isnan(session.Modulation));
            if ~isempty(mi)
                [f, x] = ecdf(mi);
                plot(ax, x, f, 'Color', colors(d,:), 'LineWidth', 2, ...
                     'DisplayName', sprintf('Day %d', d));
                hasData = true;
            end
        end

        % Formatting
        title(ax, uMice{m}, 'Interpreter', 'none');
        set(ax, 'TickDir', 'out', 'box', 'off');
        if m==1, ylabel(ax, 'Cumul. Probability'); end
        if m==round(nMice/2), xlabel(ax, 'Modulation Index'); end
        if hasData
             xline(ax, 1, 'k:', 'HandleVisibility', 'off');
             yline(ax, 0.5, 'k:', 'HandleVisibility', 'off');
        end
    end
    
    % Add legend to last plot
    if isgraphics(ax), legend(ax, 'Location', 'best'); end
    linkaxes(axList(isgraphics(axList)), 'x');
    
    title(t, sprintf('Modulation CDF (%s) - Days %s', typeToPlot, num2str(daysToPlot)), ...
        'FontSize', 14);

    if ~isempty(p.Results.SavePath), saveas(figHandle, p.Results.SavePath); end
end