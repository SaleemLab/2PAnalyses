function figHandle = plotModulation_Histogram(allData, varargin)
% plotModulation_Histogram Plots histograms of (Max-Min)/Mean for Mice x Days.

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addParameter(p, 'DaysToPlot', [1 2 3 4 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @ischar);
    addParameter(p, 'ExcludeMice', {}, @iscellstr);
    addParameter(p, 'SavePath', '', @ischar);
    parse(p, allData, varargin{:});
    
    daysToPlot = p.Results.DaysToPlot;
    typeToPlot = p.Results.TypeToPlot;
    uMice = setdiff(unique({allData.MouseID}), p.Results.ExcludeMice, 'stable');
    [nMice, nDays] = deal(length(uMice), length(daysToPlot));

    figHandle = figure('Position', [100 100 300*nDays 250*nMice]);
    t = tiledlayout(figHandle, nMice, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
    axList = gobjects(0);
    
    allValues = []; % Store all data for global scaling

    for m = 1:nMice
        for d = 1:nDays
            ax = nexttile; axList(end+1) = ax;
            
            session = allData(strcmp({allData.MouseID}, uMice{m}) & ...
                              [allData.Day] == daysToPlot(d) & ...
                              strcmpi({allData.Type}, typeToPlot));
            if isempty(session)
                axis(ax, 'off'); 
                if m==1, title(ax, ['Day ' num2str(daysToPlot(d))]); end
                continue;
            end
            session = session(1);
       
            mi = session.Modulation(~isnan(session.Modulation));
            
            % --- FIX PART 1: Collect all data ---
            allValues = [allValues; mi(:)]; 
            
            if ~isempty(mi)
                 % Use specific bin edges for consistency across plots
                 % e.g., 0 to 3 in steps of 0.1
                 histogram(ax, mi, 'BinWidth', 0.1, 'Normalization', 'count', ...
                     'FaceColor', 'k', 'EdgeColor', 'none');
                     
                 text(ax, 0.95, 0.95, sprintf('%d laps', session.NumLaps), ...
                     'Units', 'normalized', 'Horiz', 'right', 'Vert', 'top', 'Color', [.3 .3 .3]);
            else
                 text(ax, 0.5, 0.5, 'No Valid Cells', 'Horiz', 'center');
            end

            set(ax, 'TickDir', 'out', 'box', 'off');
            if m==1, title(ax, ['Day ' num2str(daysToPlot(d))]); end
            if d==1, ylabel(ax, uMice{m}, 'Interpreter', 'none', 'FontSize', 10); end
            if m==nMice && d==round(nDays/2), xlabel(ax, 'Modulation Index'); end
        end
    end
    
    % --- FIX PART 2: Robust Global Limits ---
    if ~isempty(allValues)
        % Find 99th percentile to ignore extreme outliers
        xMax = prctile(allValues, 99);
        % Ensure it's at least 1.0 so empty/low plots don't look weird
        xMax = max(xMax, 1.0); 
        
        % Apply to all valid axes
        set(axList(isgraphics(axList)), 'XLim', [0, xMax]);
    end
    
    title(t, sprintf('Modulation Distribution (%s)', typeToPlot), 'FontSize', 14);
    if ~isempty(p.Results.SavePath), saveas(figHandle, p.Results.SavePath); end
end