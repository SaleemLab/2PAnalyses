function figHandle = plotTuningHeatmaps(allData, varargin)
% plotGridTuning Plots a grid (Mice x Days of learning) of position tuning.
% Allows filtering by specific days and data type (Boutons/Somas).

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addParameter(p, 'ExcludeMice', {'M24043', 'M24046', 'M24048', 'M24049'}, @iscellstr);
    addParameter(p, 'SavePath', '', @ischar);
    % New filtering parameters
    addParameter(p, 'DaysToPlot', [1 2 3 4 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @ischar);
    parse(p, allData, varargin{:});
    
    excludeMice = p.Results.ExcludeMice;
    daysToPlot = p.Results.DaysToPlot;
    typeToPlot = p.Results.TypeToPlot;

    % Identify Mice (exclude specifically requested ones)
    uniqueMice = setdiff(unique({allData.MouseID}), excludeMice, 'stable');
    
    [nMice, nDays] = deal(length(uniqueMice), length(daysToPlot));

    figHandle = figure('Position', [100 100 300*nDays 250*nMice]);
    t = tiledlayout(figHandle, nMice, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');

    for m = 1:nMice
        for d = 1:nDays
            ax = nexttile; 
            
            % --- EXACT FILTERING ---
            % Find data matching: Mouse AND Day AND Type
            thisSession = allData(strcmp({allData.MouseID}, uniqueMice{m}) & ...
                                  [allData.Day] == daysToPlot(d) & ...
                                  strcmpi({allData.Type}, typeToPlot));

            if isempty(thisSession)
                axis(ax, 'off'); 
                % Still label the top row even if empty, so grid makes sense
                if m==1, title(ax, ['Day ' num2str(daysToPlot(d))]); end
                continue;
            end
            
            % Handle rare duplicates (warn if it happens, take first)
            if length(thisSession) > 1
                warning('Duplicate data found for %s Day %d (%s). Using first.', ...
                    uniqueMice{m}, daysToPlot(d), typeToPlot);
                thisSession = thisSession(1); 
            end

            % Pass explicit 'ax' handle to helper
            plotTuningInAxes(ax, thisSession.OddMean, thisSession.EvenMean);
            
            % Formatting Labels
            if m==1, title(ax, ['Day ' num2str(daysToPlot(d))]); end
            if d==1, ylabel(ax, uniqueMice{m}, 'Interpreter', 'none', 'FontSize', 10); end
            
            % Only show X-labels on the bottom-most plots
            if m ~= nMice
                set(ax, 'XTickLabel', {}); 
                xlabel(ax, ''); % Remove 'Position (cm)' label from non-bottom plots
            end
        end
    end
    
    % Global Title to indicate what type we are looking at
    title(t, sprintf('%s Tuning Across Days', typeToPlot), 'FontSize', 14);

    if ~isempty(p.Results.SavePath)
        saveas(figHandle, p.Results.SavePath);
    end
end