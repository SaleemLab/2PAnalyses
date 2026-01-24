function figHandle = plotPooledPopulation(allData, targetArea, varargin)
% plotPooledPopulation Plots activity pooled across ALL mice for a specific area.
%
% Usage: 
%   plotPooledPopulation(allData, 'RSP', 'DaysToPlot', [1 3 5])
%   plotPooledPopulation(allData, 'RSP', 'TypeToPlot', 'Somas')

    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    addRequired(p, 'targetArea', @ischar);
    addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'TypeToPlot', 'Boutons', @ischar);
    addParameter(p, 'SavePath', '', @ischar);
    parse(p, allData, targetArea, varargin{:});

    daysToPlot = p.Results.DaysToPlot;
    typeRequested = p.Results.TypeToPlot;
    nDays = length(daysToPlot);

    figHandle = figure('Position', [100 100 300*nDays 400]);
    t = tiledlayout(figHandle, 1, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');

    for d = 1:nDays
        day = daysToPlot(d);
        ax = nexttile;

        % FIND ALL sessions that match: Day AND Area AND Type
        % Note: Changed 'TargetArea' to 'Area' to match getTuningData.m
        daySessions = allData([allData.Day] == day & ...
                              strcmpi({allData.TargetArea}, targetArea) & ...
                              strcmpi({allData.Type}, typeRequested));

        if isempty(daySessions)
            title(ax, ['Day ' num2str(day) ' - No Data']); 
            axis(ax, 'off'); 
            continue;
        end

        % POOL THEM: Concatenate every mouse's data for this day
        pooledOdd = vertcat(daySessions.OddMean);
        pooledEven = vertcat(daySessions.EvenMean);

        plotTuningInAxes(ax, pooledOdd, pooledEven);

        % Labels
        title(ax, ['Day ' num2str(day)]);
        if d==1
             ylabel(ax, {['Pooled ' targetArea]; ['(' typeRequested ')']}, ...
                 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if d==nDays
            colorbar; 
        end
    end

    % Global title for context
    title(t, sprintf('Pooled Population: %s (%s)', targetArea, typeRequested), ...
        'Interpreter', 'none', 'FontSize', 14);

    if ~isempty(p.Results.SavePath)
        saveas(figHandle, p.Results.SavePath);
    end
end
