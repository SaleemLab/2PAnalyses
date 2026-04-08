function figHandle = plotEvenTuningHeatmaps(allData, varargin)
% PLOTEVENTUNINGHEATMAPS Plots even laps sorted by the peak of odd laps.
    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    
    % --- Optional Inclusion/Exclusion Parameters ---
    addParameter(p, 'ExcludeMice', {}, @(x) iscell(x) || ischar(x));
    addParameter(p, 'DaysToPlot', unique([allData.Day]), @isnumeric);
    addParameter(p, 'TypeToPlot', unique({allData.Type}), @iscellstr);
    addParameter(p, 'TargetArea', unique({allData.TargetArea}), @iscellstr);
    
    % --- Optional Filtering Parameters ---
    addParameter(p, 'RatioThreshold', -Inf, @isnumeric); 
    addParameter(p, 'UseStabilityFilter', false, @islogical);
    addParameter(p, 'SavePath', "Z:\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\Tuning_EvenLaps_SortedByOdd_RatioStable.png", @ischar); 
    
    parse(p, allData, varargin{:});
    
    % --- Helper Function Calls (Defined below) ---
    excludeMice = setify(p.Results.ExcludeMice);
    daysToPlot  = p.Results.DaysToPlot;
    typesToPlot = setify(p.Results.TypeToPlot);
    areasToPlot = setify(p.Results.TargetArea);
    ratioThreshold = p.Results.RatioThreshold;
    useRatioFilter = isfinite(ratioThreshold);
    useStabilityFilter = p.Results.UseStabilityFilter;
    
    BASE_SAVE_DIR = 'Z:\ibn-vision\USERS\Sonali\Figures';
    
    allMice = unique({allData.MouseID});
    uniqueMice = setdiff(allMice, excludeMice, 'stable');
    
    if isempty(uniqueMice), figHandle = []; return; end
    
    % Use first available type/area if multiple exist
    typeToPlot = typesToPlot{1};
    areaToPlot = areasToPlot{1};
    
    [nMice, nDays] = deal(length(uniqueMice), length(daysToPlot));
    figHandle = figure('Position', [100 100 400*nDays 300*nMice], 'Color', 'w');
    t = tiledlayout(figHandle, nMice, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for m = 1:nMice
        for d = 1:nDays
            ax = nexttile; 
            
            thisSession = allData(strcmp({allData.MouseID}, uniqueMice{m}) & ...
                                  [allData.Day] == daysToPlot(d) & ...
                                  strcmpi({allData.Type}, typeToPlot) & ...
                                  strcmpi({allData.TargetArea}, areaToPlot));
            
            if isempty(thisSession)
                axis(ax, 'off'); continue;
            end
            thisSession = thisSession(1); 
            
            % --- ROI Selection & Filtering ---
            roisToKeepIdx = 1:size(thisSession.OddMean, 1);
            if useStabilityFilter && isfield(thisSession, 'lapCorr_HalvesStableIdx')
                roisToKeepIdx = intersect(roisToKeepIdx, thisSession.lapCorr_HalvesStableIdx);
            end
            if useRatioFilter && isfield(thisSession, 'ratioVarToTuningRange')
                ratioIdx = find(thisSession.ratioVarToTuningRange < ratioThreshold); 
                roisToKeepIdx = intersect(roisToKeepIdx, ratioIdx);
            end
            
            if isempty(roisToKeepIdx)
                axis(ax, 'off'); continue;
            end 
            
            % --- SORTING LOGIC: Use ODD to sort EVEN ---
            OddMeanPlot = thisSession.OddMean(roisToKeepIdx, :);
            EvenMeanPlot = thisSession.EvenMean(roisToKeepIdx, :);
            
            % Normalize Odd to find peaks
            normOdd = normalize(OddMeanPlot, 2, 'range');
            [~, peakIdx] = max(normOdd, [], 2);
            [~, sortIdx] = sort(peakIdx);
            
            % Normalize Even and apply Odd's sort order
            normEvenSorted = normalize(EvenMeanPlot(sortIdx, :), 2, 'range');
            
            % --- PLOTTING ---
            imagesc(ax, normEvenSorted);
            colormap(ax, flipud(gray));
            set(ax, 'CLim', [0 1]); 
            set(ax, 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 9);
            xline(40, 'k--', 'LineWidth', 2.5);
            xline(80, 'k--', 'LineWidth', 2.5);
            xline(120, 'k--', 'LineWidth', 2.5);
            xline(160, 'k--', 'LineWidth', 2.5);
            xticks([1 40 80 120 160]);
            xticklabels({'1', '40', '80', '120', '160', '200'});
            % Even laps are usually floor(Total/2)
            totalLaps = thisSession.NumLaps;
            evenLapCount = floor(totalLaps / 2);
            lapText = sprintf('(%d Even Laps)', evenLapCount);
            
            if m==1
                title(ax, {['Day ' num2str(daysToPlot(d))], lapText}, 'FontSize', 10); 
            else
                title(ax, lapText, 'FontSize', 9); 
            end
            
            if d==1, ylabel(ax, uniqueMice{m}, 'Interpreter', 'none', 'FontWeight', 'bold'); end
            if m ~= nMice, set(ax, 'XTickLabel', {}); end
            
            % Scale Bar (20cm)
            numNeurons = size(normEvenSorted, 1);
            line(ax, [5 25], [numNeurons+5 numNeurons+5], 'Color', 'r', 'LineWidth', 2);
        end
    end
    
    filterText = getFilterSummary(useStabilityFilter, useRatioFilter, ratioThreshold);
    title(t, sprintf('Even Laps Sorted by Odd Peaks: %s %s [%s]', typeToPlot, areaToPlot, filterText), 'FontSize', 14);
    
    % Save
    savePath = p.Results.SavePath;
    if ~isempty(savePath)
        if isempty(fileparts(savePath)), savePath = fullfile(BASE_SAVE_DIR, savePath); end
        print(figHandle, savePath, '-dpng', '-r300');
    end
end

% --- HELPER FUNCTIONS ---

function cellArray = setify(input)
    if ischar(input) || (isstring(input) && isscalar(input))
        cellArray = {char(input)};
    elseif isnumeric(input) && isempty(input)
        cellArray = {};
    elseif iscell(input)
        cellArray = input;
    else
        cellArray = {input};
    end
end

function summary = getFilterSummary(useStability, useRatio, ratioThresh)
    parts = {};
    if useStability, parts{end+1} = 'Stable'; end
    if useRatio, parts{end+1} = sprintf('Ratio < %.2g', ratioThresh); end
    if isempty(parts), summary = 'Unfiltered'; else, summary = strjoin(parts, ' & '); end
end