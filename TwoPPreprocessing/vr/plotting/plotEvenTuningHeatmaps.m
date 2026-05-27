function figHandle = plotEvenTuningHeatmaps(allData, varargin)
    p = inputParser;
    addRequired(p, 'allData', @isstruct);
    
    % --- Optional Inclusion/Exclusion Parameters ---
    addParameter(p, 'ExcludeMice', {}, @(x) iscell(x) || ischar(x));
    addParameter(p, 'DaysToPlot', [], @isnumeric); 
    addParameter(p, 'TypeToPlot', {}, @(x) iscellstr(x) || ischar(x) || isstring(x));
    
    % --- Optional Filtering Parameters ---
    addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x)); 
    
    parse(p, allData, varargin{:});
    
    % Handle Days
    daysToPlot = p.Results.DaysToPlot;
    if isempty(daysToPlot), daysToPlot = unique([allData.Day]); end
    
    % Get Types
    typesToPlot = setify(p.Results.TypeToPlot);
    if isempty(typesToPlot)
        rawTypes = {allData.TypeImaged}; 
        cleanTypes = cellfun(@char, rawTypes, 'UniformOutput', false);
        typesToPlot = unique(cleanTypes(~cellfun(@isempty, cleanTypes))); % <-- FIXED HERE
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
                                  strcmpi({allData.TypeImaged}, typeToPlot));
            
            if isempty(thisSession) || ~isfield(thisSession(1), 'ConditionData') || isempty(thisSession(1).ConditionData)
                axis(ax, 'off'); continue;
            end
            thisSession = thisSession(1); 
            
            % Get the name of the first available condition block to extract traces
            condNames = fieldnames(thisSession.ConditionData);
            if isempty(condNames), axis(ax, 'off'); continue; end
            targetCond = condNames{1}; 
            
            % Extract full raw activity matrix: [ROIs x Laps x PositionBins]
            fullActivity = thisSession.ConditionData.(targetCond).LapActivity;
            
            % -----------------------------------------------------------------
            % Grab pre-computed FilteredROIs field instead of calculating here
            % -----------------------------------------------------------------
            if isfield(thisSession, 'FilteredROIs') && ~isempty(thisSession.FilteredROIs)
                roisToKeepIdx = thisSession.FilteredROIs;
            else
                roisToKeepIdx = 1:size(fullActivity, 1);
                warning('No FilteredROIs field found for Mouse %s Day %d. Using all ROIs.', ...
                        thisSession.MouseID, thisSession.Day);
            end
            
            if isempty(roisToKeepIdx)
                title(ax, sprintf('Day %d (No Cells)', daysToPlot(d)));
                axis(ax, 'off'); continue;
            end 
            
            % -----------------------------------------------------------------
            % Compute Odd and Even spatial means from LapActivity on the fly
            % -----------------------------------------------------------------
            nTotalLaps = size(fullActivity, 2);
            oddLapIndices = 1:2:nTotalLaps;
            evenLapIndices = 2:2:nTotalLaps;
            
            % Slice by your filtered ROIs, then take the mean across laps (dim 2)
            OddMeanPlot  = squeeze(mean(fullActivity(roisToKeepIdx, oddLapIndices, :), 2, 'omitnan'));
            EvenMeanPlot = squeeze(mean(fullActivity(roisToKeepIdx, evenLapIndices, :), 2, 'omitnan'));
            
            % Ensure 2D matrices if only 1 ROI survives squeeze execution
            if size(OddMeanPlot, 2) == 1 && length(roisToKeepIdx) == 1
                OddMeanPlot = OddMeanPlot';
                EvenMeanPlot = EvenMeanPlot';
            end
            
            % --- Sorting & Plotting (Kept Exactly As Is) ---
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
            
            % --- TITLE: Day Number and Number of Laps (Kept Exactly As Is) ---
            nEven = length(evenLapIndices);
            nOdd = length(oddLapIndices);
            
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

% function figHandle = plotEvenTuningHeatmaps(allData, varargin)
%     p = inputParser;
%     addRequired(p, 'allData', @isstruct);
% 
%     % --- Optional Inclusion/Exclusion Parameters ---
%     addParameter(p, 'ExcludeMice', {}, @(x) iscell(x) || ischar(x));
%     addParameter(p, 'DaysToPlot', [], @isnumeric); 
%     addParameter(p, 'TypeToPlot', {}, @(x) iscellstr(x) || ischar(x) || isstring(x));
% 
%     % --- Optional Filtering Parameters ---
%     addParameter(p, 'UseStabilityFilter', false, @islogical); % Now default false
%     addParameter(p, 'HalvesRhoThreshold_boutons', 0.6, @isnumeric); 
%     addParameter(p, 'HalvesRhoThreshold_somas', 0.8, @isnumeric); 
%     addParameter(p, 'SavePath', '', @(x) ischar(x) || isstring(x)); 
% 
%     parse(p, allData, varargin{:});
% 
%     % Handle Days
%     daysToPlot = p.Results.DaysToPlot;
%     if isempty(daysToPlot), daysToPlot = unique([allData.Day]); end
% 
%     % Get Types
%     typesToPlot = setify(p.Results.TypeToPlot);
%     if isempty(typesToPlot)
%         rawTypes = {allData.Type};
%         cleanTypes = cellfun(@char, rawTypes, 'UniformOutput', false);
%         typesToPlot = unique(cleanTypes(~cellfun(@isempty, cleanTypes)));
%     end
% 
%     excludeMice = setify(p.Results.ExcludeMice);
%     allMice = unique({allData.MouseID});
%     uniqueMice = setdiff(allMice, excludeMice, 'stable');
% 
%     if isempty(uniqueMice), figHandle = []; return; end
% 
%     typeToPlot = typesToPlot{1}; 
%     [nMice, nDays] = deal(length(uniqueMice), length(daysToPlot));
% 
%     figHandle = figure('Position', [100 100 400*nDays 300*nMice], 'Color', 'w');
%     t = tiledlayout(figHandle, nMice, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
%     for m = 1:nMice
%         for d = 1:nDays
%             ax = nexttile; 
% 
%             thisSession = allData(strcmp({allData.MouseID}, uniqueMice{m}) & ...
%                                   [allData.Day] == daysToPlot(d) & ...
%                                   strcmpi({allData.Type}, typeToPlot));
% 
%             if isempty(thisSession) || isempty(thisSession(1).OddMean)
%                 axis(ax, 'off'); continue;
%             end
%             thisSession = thisSession(1); 
% 
%             % 
%             roisToKeepIdx = 1:size(thisSession.OddMean, 1);
% 
%             if strcmpi(thisSession.TypeImaged, 'Boutons')
%             % Optional Stability Filter with 
%                 if p.Results.UseStabilityFilter && isfield(thisSession, 'lapCorr_HalvesRho') ...
%                     && isfield(thisSession, 'cvExpVar') && isfield(thisSession, 'uniqueBoutonIdx')
% 
%                     % Use the rho value and the threshold passed in varargin
%                     halvesStableIdx = find(thisSession.lapCorr_HalvesRho >= HalvesRhoThreshold_boutons);
%                     roisToKeepIdx = intersect(roisToKeepIdx, halvesStableIdx);
%                 end
%             elseif strcmpi(thisSession.TypeImaged, 'Somas')
%                 if p.Results.UseStabilityFilter && isfield(thisSession, 'lapCorr_HalvesRho') ...
%                     && isfield(thisSession, 'cvExpVar')
% 
%                     cvExpVar_nFoldMedian = median(thisSession.cvExpVar, 'omitnan');
% 
%                     halvesStableIdx = find(thisSession.lapCorr_HalvesRho >= HalvesRhoThreshold_somas);
% 
%                     roisToKeepIdx = intersect(roisToKeepIdx, halvesStableIdx);
%                 end
%             end 
% 
%             if isempty(roisToKeepIdx)
%                 title(ax, sprintf('Day %d (No Cells)', daysToPlot(d)));
%                 axis(ax, 'off'); continue;
%             end 
% 
%             % --- Sorting & Plotting ---
%             OddMeanPlot = thisSession.OddMean(roisToKeepIdx, :);
%             EvenMeanPlot = thisSession.EvenMean(roisToKeepIdx, :);
% 
%             normOdd = normalize(OddMeanPlot, 2, 'range');
%             [~, peakIdx] = max(normOdd, [], 2);
%             [~, sortIdx] = sort(peakIdx);
% 
%             normEvenSorted = normalize(EvenMeanPlot(sortIdx, :), 2, 'range');
%             normEvenSorted(isnan(normEvenSorted)) = 0;
% 
%             imagesc(ax, normEvenSorted);
%             colormap(ax, flipud(gray));
%             set(ax, 'CLim', [0.25 0.75], 'TickDir', 'out', 'YDir', 'normal', 'FontSize', 8);
% 
%             % Horizontal grid/spatial lines
%             for xVal = 40:40:160, xline(ax, xVal, 'k--', 'LineWidth', 1); end
% 
%             % --- TITLE: Day Number and Number of Laps ---
%             nTotal = thisSession.NumLaps;
%             nEven = floor(nTotal/2);
%             nOdd = ceil(nTotal/2);
% 
%             titleStr = {sprintf('Day %d', daysToPlot(d)), ...
%                         sprintf('(%d Odd / %d Even)', nOdd, nEven)};
%             title(ax, titleStr, 'FontSize', 9);
% 
%             if d == 1, ylabel(ax, uniqueMice{m}, 'FontWeight', 'bold', 'Interpreter', 'none'); end
%         end
%     end
% 
%     % Main figure title
%     title(t, sprintf('Even Laps (Sorted by Odd Peaks) - %s', typeToPlot), 'FontSize', 12, 'FontWeight', 'bold');
% 
%     if ~isempty(p.Results.SavePath)
%         exportgraphics(figHandle, p.Results.SavePath, 'Resolution', 300);
%     end
% end
% 
% function cellArray = setify(input)
%     if isempty(input), cellArray = {};
%     elseif ischar(input), cellArray = {input};
%     elseif isstring(input), cellArray = cellstr(input);
%     elseif iscell(input), cellArray = input;
%     else, cellArray = {input};
%     end
% end