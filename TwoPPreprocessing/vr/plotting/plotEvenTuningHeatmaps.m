function figHandle = plotEvenTuningHeatmaps(allData, varargin)
    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    
    % --- Optional Inclusion/Exclusion Parameters ---
    addParameter(p, 'ExcludeMice', {}, @(x) iscell(x) || ischar(x));
    addParameter(p, 'DaysToPlot', [], @isnumeric); 
    addParameter(p, 'TypeToPlot', {}, @(x) iscellstr(x) || ischar(x) || isstring(x));
    
    % --- Optional Filtering Parameters ---
    addParameter(p, 'UseStabilityFilter', false, @islogical); % Now default false
    addParameter(p, 'StableThresh', 0.4, @isnumeric);         % Custom threshold
    addParameter(p, 'RatioThreshold', -Inf, @isnumeric);      % Default -Inf (off)
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x)); 
    
    parse(p, allData, varargin{:});
    
    % Handle Days
    daysToPlot = p.Results.DaysToPlot;
    if isempty(daysToPlot), daysToPlot = unique([allData.Day]); end
    
    % Get Types
    typesToPlot = setify(p.Results.TypeToPlot);
    if isempty(typesToPlot)
        rawTypes = {allData.Type};
        cleanTypes = cellfun(@char, rawTypes, 'UniformOutput', false);
        typesToPlot = unique(cleanTypes(~cellfun(@isempty, cleanTypes)));
    end
    
    excludeMice = setify(p.Results.ExcludeMice);
    allMice = unique({allData.MouseID});
    uniqueMice = setdiff(allMice, excludeMice, 'stable');
    
    if isempty(uniqueMice), figHandle = []; return; end
    
    typeToPlot = typesToPlot{1}; 
    [nMice, nDays] = deal(length(uniqueMice), length(daysToPlot));
    
    figHandle = figure('Position', [100 100 400*nDays 300*nMice], 'Color', 'w');
    t = tiledlayout(figHandle, nMice, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for m = 1:nMice
        for d = 1:nDays
            ax = nexttile; 
            
            thisSession = allData(strcmp({allData.MouseID}, uniqueMice{m}) & ...
                                  [allData.Day] == daysToPlot(d) & ...
                                  strcmpi({allData.Type}, typeToPlot));
            
            if isempty(thisSession) || isempty(thisSession(1).OddMean)
                axis(ax, 'off'); continue;
            end
            thisSession = thisSession(1); 
            
            % --- ROI Selection & Filtering ---
            roisToKeepIdx = 1:size(thisSession.OddMean, 1);
            
            % Optional Stability Filter with custom threshold
            if p.Results.UseStabilityFilter && isfield(thisSession, 'lapCorr_HalvesRho')
                % Use the rho value and the threshold passed in varargin
                stableIdx = find(thisSession.lapCorr_HalvesRho >= p.Results.StableThresh);
                roisToKeepIdx = intersect(roisToKeepIdx, stableIdx);
            end
            
            % Optional Ratio Filter
            if p.Results.RatioThreshold ~= -Inf && isfield(thisSession, 'ratioVarToTuningRange')
                ratioIdx = find(thisSession.ratioVarToTuningRange < p.Results.RatioThreshold); 
                roisToKeepIdx = intersect(roisToKeepIdx, ratioIdx);
            end
            
            if isempty(roisToKeepIdx)
                title(ax, sprintf('Day %d (No Cells)', daysToPlot(d)));
                axis(ax, 'off'); continue;
            end 
            
            % --- Sorting & Plotting ---
            OddMeanPlot = thisSession.OddMean(roisToKeepIdx, :);
            EvenMeanPlot = thisSession.EvenMean(roisToKeepIdx, :);
            
            normOdd = normalize(OddMeanPlot, 2, 'range');
            [~, peakIdx] = max(normOdd, [], 2);
            [~, sortIdx] = sort(peakIdx);
            
            normEvenSorted = normalize(EvenMeanPlot(sortIdx, :), 2, 'range');
            normEvenSorted(isnan(normEvenSorted)) = 0;
            
            imagesc(ax, normEvenSorted);
            colormap(ax, flipud(gray));
            set(ax, 'CLim', [0.25 0.75], 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 8);
            
            % Horizontal grid/spatial lines
            for xVal = 40:40:160, xline(ax, xVal, 'k--', 'LineWidth', 1); end
            
            % --- TITLE: Day Number and Number of Laps ---
            nTotal = thisSession.NumLaps;
            nEven = floor(nTotal/2);
            nOdd = ceil(nTotal/2);
            
            titleStr = {sprintf('Day %d', daysToPlot(d)), ...
                        sprintf('(%d Odd / %d Even)', nOdd, nEven)};
            title(ax, titleStr, 'FontSize', 9);
            
            if d == 1, ylabel(ax, uniqueMice{m}, 'FontWeight', 'bold', 'Interpreter', 'none'); end
        end
    end
    
    % Main figure title
    title(t, sprintf('Even Laps (Sorted by Odd Peaks) - %s', typeToPlot), 'FontSize', 12, 'FontWeight', 'bold');
    
    if ~isempty(p.Results.SavePath)
        exportgraphics(figHandle, p.Results.SavePath, 'Resolution', 300);
    end
end

function cellArray = setify(input)
    if isempty(input), cellArray = {};
    elseif ischar(input), cellArray = {input};
    elseif isstring(input), cellArray = cellstr(input);
    elseif iscell(input), cellArray = input;
    else, cellArray = {input};
    end
end