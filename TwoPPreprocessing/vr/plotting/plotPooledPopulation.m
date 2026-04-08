% function figHandle = plotPooledPopulation(allData, targetArea, varargin)
% % plotPooledPopulation Plots activity pooled across ALL mice for a specific area.
% %
% % Usage: 
% %   plotPooledPopulation(allData, 'RSP', 'DaysToPlot', [1 3 5])
% %   plotPooledPopulation(allData, 'RSP', 'TypeToPlot', 'Somas')
% 
%     p = inputParser;
%     addRequired(p, 'allData', @isstruct);
%     addRequired(p, 'targetArea', @ischar);
%     addParameter(p, 'DaysToPlot', [1, 2, 3, 4, 5], @isnumeric);
%     addParameter(p, 'TypeToPlot', 'Boutons', @ischar);
%     addParameter(p, 'SavePath', '', @ischar);
%     parse(p, allData, targetArea, varargin{:});
% 
%     daysToPlot = p.Results.DaysToPlot;
%     typeRequested = p.Results.TypeToPlot;
%     nDays = length(daysToPlot);
% 
%     figHandle = figure('Position', [100 100 300*nDays 400]);
%     t = tiledlayout(figHandle, 1, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
%     for d = 1:nDays
%         day = daysToPlot(d);
%         ax = nexttile;
% 
%         % FIND ALL sessions that match: Day AND Area AND Type
%         % Note: Changed 'TargetArea' to 'Area' to match getTuningData.m
%         daySessions = allData([allData.Day] == day & ...
%                               strcmpi({allData.TargetArea}, targetArea) & ...
%                               strcmpi({allData.Type}, typeRequested));
% 
%         if isempty(daySessions)
%             title(ax, ['Day ' num2str(day) ' - No Data']); 
%             axis(ax, 'off'); 
%             continue;
%         end
% 
%         % POOL THEM: Concatenate every mouse's data for this day
%         pooledOdd = vertcat(daySessions.OddMean);
%         pooledEven = vertcat(daySessions.EvenMean);
% 
%         plotTuningInAxes(ax, pooledOdd, pooledEven);
% 
%         % Labels
%         title(ax, ['Day ' num2str(day)]);
%         if d==1
%              ylabel(ax, {['Pooled ' targetArea]; ['(' typeRequested ')']}, ...
%                  'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
%         end
%         if d==nDays
%             colorbar; 
%         end
%     end
% 
%     % Global title for context
%     title(t, sprintf('Pooled Population: %s (%s)', targetArea, typeRequested), ...
%         'Interpreter', 'none', 'FontSize', 14);
% 
%     if ~isempty(p.Results.SavePath)
%         saveas(figHandle, p.Results.SavePath);
%     end
% end

function figHandle = plotPooledPopulation(allData, targetArea, varargin)
% plotPooledPopulation Plots activity pooled across ALL mice for a specific area.
% Uses Odd-lap peaks to sort Even-lap activity for stable ROIs.
    
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
    
    figHandle = figure('Position', [100 100 300*nDays 450], 'Color', 'w');
    t = tiledlayout(figHandle, 1, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for d = 1:nDays
        day = daysToPlot(d);
        ax = nexttile;
        
        % 1. Find sessions matching criteria
        daySessions = allData([allData.Day] == day & ...
                              strcmpi({allData.TargetArea}, targetArea) & ...
                              strcmpi({allData.Type}, typeRequested));
        
        if isempty(daySessions)
            title(ax, ['Day ' num2str(day) ' - No Data']); 
            axis(ax, 'off'); 
            continue;
        end
        
        % 2. Collect Stable ROIs from both halves
        allStableOdd = [];
        allStableEven = [];
        
        for s = 1:length(daySessions)
            thisSession = daySessions(s);
            
            % Use your stability index to filter
            if isfield(thisSession, 'lapCorr_HalvesStableIdx') && ...
               ~isempty(thisSession.lapCorr_HalvesStableIdx)
                
                idx = thisSession.lapCorr_HalvesStableIdx;
                
                % Pool the specific ROIs
                allStableOdd = vertcat(allStableOdd, thisSession.OddMean(idx, :));
                allStableEven = vertcat(allStableEven, thisSession.EvenMean(idx, :));
            end
        end
        
        % 3. Apply Cross-Validated Sorting
        if ~isempty(allStableEven)
            % Find peaks in ODD to define the order
            [~, peakBins] = max(allStableOdd, [], 2);
            [~, sortIdx] = sort(peakBins);
            
            % Apply that order to BOTH for the plotting function
            sortedOdd = allStableOdd(sortIdx, :);
            sortedEven = allStableEven(sortIdx, :);
            
            % 4. Use your usual plotting function
            % Passing sortedEven as both arguments if you only want to see Even,
            % or sortedOdd/sortedEven if your function plots them side-by-side.
            plotTuningInAxes(ax, sortedOdd, sortedEven);
            
            % Labels
            title(ax, sprintf('Day %d (n=%d)', day, size(sortedEven, 1)));
        else
            title(ax, ['Day ' num2str(day) ' - 0 Stable']);
        end
        
        if d == 1
             ylabel(ax, {['Pooled ' targetArea]; ['(' typeRequested ')']}, ...
                 'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
        end
        
        if d == nDays
            colorbar; 
        end
    end
    
    % Global title
    title(t, sprintf('Pooled %s (%s): Even Laps Sorted by Odd Peaks', targetArea, typeRequested), ...
        'Interpreter', 'none', 'FontSize', 14);
    
    if ~isempty(p.Results.SavePath)
        saveas(figHandle, p.Results.SavePath);
    end
end