function figHandle = plotTuningHeatmaps(allData, varargin)
% plotTuningHeatmaps Plots a grid (Mice x Days) of position tuning using flexible filters.
%
% Usage: 
%   1. Plot everything: plotTuningHeatmaps(allData)
%   2. Exclude specific mice: plotTuningHeatmaps(allData, 'ExcludeMice', {'M1', 'M2'})
%   3. Specify Days/Type: plotTuningHeatmaps(allData, 'DaysToPlot', [1 3 5], 'TypeToPlot', 'Somas')
%   4. Apply Ratio Filter: plotTuningHeatmaps(allData, 'RatioThreshold', 0.5)
%
% NOTE: This function assumes the existence of 'plotTuningInAxes(ax, OddMean, EvenMean)'
%       to perform the core plotting (sort by Odd, plot Even).

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
    addParameter(p, 'SavePath', "Z:\ibn-vision\USERS\Sonali\Figures\Tuning_dataclub2026_22.png", @ischar); 
    
    parse(p, allData, varargin{:});
    
    % Get the results, ensuring cell arrays are handled
    excludeMice = setify(p.Results.ExcludeMice);
    daysToPlot  = p.Results.DaysToPlot;
    typesToPlot = setify(p.Results.TypeToPlot);
    areasToPlot = setify(p.Results.TargetArea);
    
    ratioThreshold = p.Results.RatioThreshold;
    useRatioFilter = isfinite(ratioThreshold);
    useStabilityFilter = p.Results.UseStabilityFilter;
    
    % --- Define Base Save Directory ---
    BASE_SAVE_DIR = 'Z:\ibn-vision\USERS\Sonali\Figures';
    
    % --- 1. Identify Mice (Exclude specified ones) ---
    allMice = unique({allData.MouseID});
    uniqueMice = setdiff(allMice, excludeMice, 'stable');
    
    if isempty(uniqueMice)
        warning('No mice remaining after exclusion filter.');
        figHandle = []; return;
    end
    
    % --- 2. Determine Plotting Grid ---
    if length(typesToPlot) ~= 1 || length(areasToPlot) ~= 1
        error('This function must be run with exactly ONE TypeToPlot and ONE TargetArea for the grid layout.');
    end
    typeToPlot = typesToPlot{1};
    areaToPlot = areasToPlot{1};
    
    [nMice, nDays] = deal(length(uniqueMice), length(daysToPlot));
    
    % Use a generous figure size to improve on-screen look before saving
    figHandle = figure('Position', [100 100 400*nDays 300*nMice]);
    t = tiledlayout(figHandle, nMice, nDays, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- 3. Iterate and Plot ---
    for m = 1:nMice
        for d = 1:nDays
            ax = nexttile; 
            
            % --- SESSION FILTERING (Input Metadata) ---
            thisSession = allData(strcmp({allData.MouseID}, uniqueMice{m}) & ...
                                  [allData.Day] == daysToPlot(d) & ...
                                  strcmpi({allData.Type}, typeToPlot) & ...
                                  strcmpi({allData.TargetArea}, areaToPlot));
            
            if isempty(thisSession)
                axis(ax, 'off'); 
                if m==1, title(ax, ['Day ' num2str(daysToPlot(d))]); end
                continue;
            end
            
            % Get the correct session (handle duplicates)
            thisSession = thisSession(1); 
            
            % --- ROI DATA INITIALIZATION ---
            OddMeanPlot = thisSession.OddMean;
            EvenMeanPlot = thisSession.EvenMean;
            roisToKeepIdx = 1:size(OddMeanPlot, 1);
            
            % --- FILTERING LOGIC 1: STABILITY INDEX ---
            if useStabilityFilter && isfield(thisSession, 'lapCorr_HalvesStableIdx') && ~isempty(thisSession.lapCorr_HalvesStableIdx)
                stableRoiIndices = thisSession.lapCorr_HalvesStableIdx;
                roisToKeepIdx = intersect(roisToKeepIdx, stableRoiIndices);
            end
            
            % --- FILTERING LOGIC 2: RATIO THRESHOLD ---
            if useRatioFilter
                if isfield(thisSession, 'ratioVarToTuningRange') && ~isempty(thisSession.ratioVarToTuningRange)
                    ratioValues = thisSession.ratioVarToTuningRange;
                    
                    % Keep indices where the ratio IS LESS THAN the specified threshold
                    ratioIdx = find(ratioValues < ratioThreshold); 
                    
                    roisToKeepIdx = intersect(roisToKeepIdx, ratioIdx);
                else
                    warning('Session %s Day %d: Ratio field is missing. Skipping ratio filter.', ...
                        uniqueMice{m}, daysToPlot(d));
                end 
            end
            
            % --- FINAL CHECK AND PLOTTING ---
            if isempty(roisToKeepIdx)
                warning('Session %s Day %d: No ROIs satisfy ALL active filters. Skipping plot.', ...
                    uniqueMice{m}, daysToPlot(d));
                axis(ax, 'off');
                continue;
            end 
            
            % Apply the final combined indices
            OddMeanPlot = OddMeanPlot(roisToKeepIdx, :);
            EvenMeanPlot = EvenMeanPlot(roisToKeepIdx, :);
            
            % Call helper function 
            plotTuningInAxes(ax, OddMeanPlot, EvenMeanPlot); 
            
            % --- POSITIONAL SCALE BAR ADDITION ---
            numNeurons = size(EvenMeanPlot, 1);
            scaleLength = 20; % 20 cm, assuming 1 bin = 1 cm
            % Position the scale bar slightly below the last ROI (adjust +10 margin as needed)
            yPosition = numNeurons + 5; 
            xStart = 5; % 5 cm from the left edge
            xEnd = xStart + scaleLength;

            % Draw the white line segment
            line(ax, [xStart xEnd], [yPosition yPosition], 'Color', 'w', 'LineWidth', 5);
            
            % Add text label
            text(ax, xStart + scaleLength/2, yPosition + 10, '20 cm', 'Color', 'w', 'FontSize', 10, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            
            % Formatting Labels
            lapText = sprintf(' (Laps: %d)', thisSession.NumLaps);
            if m==1
                title(ax, ['Day ' num2str(daysToPlot(d)) lapText]); 
            else
                title(ax, lapText); 
            end
            
            if d==1, ylabel(ax, uniqueMice{m}, 'Interpreter', 'none', 'FontSize', 10); end
            
            % Only show X-labels on the bottom-most plots
            if m ~= nMice
                set(ax, 'XTickLabel', {}); 
                xlabel(ax, ''); 
            end
        end
    end
    
    % Global Title
    filterText = getFilterSummary(useStabilityFilter, useRatioFilter, ratioThreshold);
    title(t, sprintf('%s %s Tuning Heatmaps (%s)', typeToPlot, areaToPlot, filterText), 'FontSize', 14);
    
    % --- SAVE LOGIC MODIFICATION (Using 'print' for high DPI) ---
    savePath = p.Results.SavePath;
    if ~isempty(savePath)
        % Check if the provided path is just a filename (no directory components)
        if isempty(fileparts(savePath))
            % If it's just a filename, prepend the base directory
            finalSavePath = fullfile(BASE_SAVE_DIR, savePath);
        else
            % If a full or partial path was provided, use it as is
            finalSavePath = savePath;
        end
        
        % Ensure the directory exists before saving
        saveDir = fileparts(finalSavePath);
        if ~isempty(saveDir) && ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        
        try
            % Use 'print' command for high DPI (e.g., 300 DPI for PNG)
            print(figHandle, finalSavePath, '-dpng', '-r300'); 
            fprintf('Figure saved successfully to: %s (300 DPI)\n', finalSavePath);
        catch ME
            warning('Failed to save figure to %s. Error: %s', finalSavePath, ME.message);
        end
    end
end 

% --- HELPER FUNCTIONS ---

function cellArray = setify(input)
    % Converts single char/string input to cell array if not already one.
    if ischar(input) || (isstring(input) && isscalar(input))
        cellArray = {char(input)};
    elseif isnumeric(input) && isempty(input)
        cellArray = {};
    elseif iscell(input)
        cellArray = input;
    else
        error('Invalid input type for inclusion/exclusion parameter.');
    end
end

function summary = getFilterSummary(useStability, useRatio, ratioThresh)
    % Generates a summary string of active filters for the title.
    parts = {};
    if useStability
        parts{end+1} = 'Stable ROIs';
    end
    if useRatio
        parts{end+1} = sprintf('Ratio < %.2g', ratioThresh);
    end
    if isempty(parts)
        summary = 'Unfiltered';
    else
        summary = strjoin(parts, ' & ');
    end
end