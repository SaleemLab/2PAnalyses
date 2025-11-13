function figHandle = computeModulationDistribution_PlotGridAcrossMice(filteredTable, varargin)
%
%   figHandle = computeModulationDistribution_PlotGridAcrossMice(filteredTable, ...)
%   Computes a modulation metric ((max-min)/mean) for each cell's median
%   spatial tuning curve and plots the distribution for each session.
%   The grid is structured as Days (rows) x Mice (columns).
%
%   Returns a handle to the generated figure.
%
%   Additional parameters:
%       'signalToUse'         - String, (default: 'dFFNeuropilCorrected')
%       'applySmoothing'      - Logical, (default: true)
%       'ExcludeMice'         - Cell array of strings, (default: {})
%       'SavePath'            - String, full file path (e.g., 'C:\Fig.png')
%                               If provided, saves the figure. (default: '')
%       'GCaMPInjectionSite'  - Column name for injection site (default:
%       'GCaMPInjectionSite'): RSPp, RSPa, ENTm
% 
    %% Parse inputs
    p = inputParser;
    addRequired(p, 'filteredTable', @istable);
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    addParameter(p, 'ExcludeMice', {}, @iscellstr);
    addParameter(p, 'SavePath', '', @ischar);
    addParameter(p, 'GCaMPInjectionSiteColumn', 'GCaMPInjectionSite', @ischar);
    parse(p, filteredTable, varargin{:});
    signalToUse = p.Results.signalToUse;
    applySmoothing = p.Results.applySmoothing;
    miceToExclude = p.Results.ExcludeMice;
    savePath = p.Results.SavePath;
    injectionSiteColumn = p.Results.GCaMPInjectionSiteColumn;
    %% Identify grid layout from table 
    requiredCols = {'MouseID', 'Session', 'DayOfExperience'};
    if ~all(ismember(requiredCols, filteredTable.Properties.VariableNames))
        error('filteredTable must contain: %s', strjoin(requiredCols, ', '));
    end
    hasInjectionSite = ismember(injectionSiteColumn, filteredTable.Properties.VariableNames);
    if ~hasInjectionSite
        warning('Column ''%s'' not found. Will not plot injection site.', injectionSiteColumn);
    end
    uniqueMice = unique(filteredTable.MouseID, 'stable');
    uniqueMice = setdiff(uniqueMice, miceToExclude, 'stable');
    numMice = length(uniqueMice);
    uniqueDays = unique(filteredTable.DayOfExperience, 'stable');
    numDays = length(uniqueDays);
    if numMice == 0 || numDays == 0
        warning('filteredTable, after exclusions, is empty or has no valid mice/days.');
        figHandle = [];
        return;
    end
    % Grid is Days (rows) x Mice (cols)
    figHandle = figure('Position', [100 100 300*numMice 250*numDays]);
    tl = tiledlayout(figHandle, numDays, numMice, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    if applySmoothing
        w = gausswin(6);
        w = w / sum(w);
    end
    % Store axes handles for resizing
    ax = gobjects(numDays, numMice);
    allMetrics = cell(numDays, numMice); % To store metrics for global scaling
    %% Iterate through grid and plot
    for thisDay = 1:numDays
        dayNum = uniqueDays(thisDay);
        
        for iMouse = 1:numMice
            mouseID = uniqueMice{iMouse};
            
            ax(thisDay, iMouse) = nexttile; % Get and store axis handle
            
            % Find session row for this specific day and mouse
            sessionRow = filteredTable( ...
                strcmp(filteredTable.MouseID, mouseID) & ...
                filteredTable.DayOfExperience == dayNum, :);
            if isempty(sessionRow)
                title(''); axis off; continue;
            end
            
            % Get injection site info
            injectionSite = '';
            if hasInjectionSite
                mouseRows = filteredTable(strcmp(filteredTable.MouseID, mouseID), :);
                if ~isempty(mouseRows)
                    injectionSite = char(mouseRows.(injectionSiteColumn)(1));
                end
            end
            if height(sessionRow) > 1 
                warning('Multiple entries for %s - Day %d. Using first.', mouseID, dayNum);
                sessionRow = sessionRow(1, :);
            end
            
            
            sessionString_char = char(sessionRow.Session);
            clear sessionFileInfo response
            try
                sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID, sessionString_char);
                if ~isfile(sessionFileInfoFilePath)
                    title(''); axis off;
                    warning('Info File Missing for %s, %s', mouseID, sessionString_char);
                    continue;
                end
                load(sessionFileInfoFilePath, 'sessionFileInfo');
            catch ME
                title(''); axis off;
                warning('Error loading sessionFileInfo for %s, %s: %s', mouseID, sessionString_char, ME.message);
                continue;
            end
            
            stimNames = string({sessionFileInfo.stimFiles.name});
            responseFileIdx = [];
            idx_Combined = find(contains(stimNames, "VRCorr") & contains(stimNames, "CombinedRuns"));
            idx_VrOnly = find(contains(stimNames, "VRCorr") & ~contains(stimNames, "CombinedRuns"));
            if ~isempty(idx_Combined)
                responseFileIdx = idx_Combined(1);
            elseif ~isempty(idx_VrOnly)
                responseFileIdx = idx_VrOnly(end);
            end
            if isempty(responseFileIdx)
                title(''); axis off;
                warning('No valid VRCorr file found for %s, %s', mouseID, sessionString_char);
                continue;
            end
            
             if ~isfield(sessionFileInfo.stimFiles, 'Response') || isempty(sessionFileInfo.stimFiles(responseFileIdx).Response)
                 warning('Response field missing or empty for selected stimFile in %s, %s', mouseID, sessionString_char);
                 title(''); axis off; continue;
             end
            responseFilePath = sessionFileInfo.stimFiles(responseFileIdx).Response;
            
            try
                if ~isfile(responseFilePath)
                     title(''); axis off;
                     warning('Response file not found at path: %s', responseFilePath);
                     continue;
                end
                load(responseFilePath, 'response');
            catch ME
                title(''); axis off;
                warning('Error loading response file for %s, %s: %s', mouseID, sessionString_char, ME.message);
                continue;
            end
            if ~exist('response', 'var')
                title(''); axis off;
                warning('File loaded, but ''response'' variable not found in: %s', responseFilePath);
                continue;
            end
            
            if ~isfield(response, 'lapPositionActivity') || ...
               ~isfield(response.lapPositionActivity, signalToUse) || ...
               isempty(response.lapPositionActivity.(signalToUse))
                 title(''); axis off;
                 warning('Data Missing in response file for %s, %s', mouseID, sessionString_char);
                 continue;
            end
            lapActivity = response.lapPositionActivity.(signalToUse);
            
            %% --- ADDED ---
            % Get the number of laps
            numLaps = size(lapActivity, 2);

            %% --- MODIFIED ---
            % Use the numLaps variable in the check
            if numLaps < 2
                title(''); axis off;
                warning('(<2 laps) for %s, %s', mouseID, sessionString_char);
                continue;
            end

            % Apply Spatial Smoothing
            if applySmoothing
                for iCell = 1:size(lapActivity, 1)
                    for iLap = 1:size(lapActivity, 2)
                        trace = squeeze(lapActivity(iCell, iLap, :));
                        if all(isnan(trace)), continue; end
                        nanMask = isnan(trace);
                        trace(nanMask) = 0;
                        smoothed = filtfilt(w, 1, trace);
                        smoothed(nanMask) = NaN;
                        lapActivity(iCell, iLap, :) = smoothed;
                    end
                end
            end
            % Calculate Modulation Metric
            medianAllLaps = squeeze(mean(lapActivity, 2, 'omitnan')); 
            if size(lapActivity, 1) == 1
                medianAllLaps = reshape(medianAllLaps, 1, []);
            end
            cellMax = max(medianAllLaps, [], 2, 'omitnan');
            cellMin = min(medianAllLaps, [], 2, 'omitnan');
            cellMean = mean(medianAllLaps, 2, 'omitnan');
            modulationMetric = (cellMax - cellMin) ./ cellMean;
            
            modulationMetric(cellMean == 0) = NaN;
            modulationMetric(isnan(cellMean)) = NaN; 
            
            modulationMetric = modulationMetric(~isnan(modulationMetric)); 
            allMetrics{thisDay, iMouse} = modulationMetric; % Store for global scaling
            % Plotting (Histogram)
            if ~isempty(modulationMetric)
                histogram(modulationMetric, 20, 'Normalization', 'count', ...
                    'EdgeColor', 'none', 'FaceColor', 'k');
                
                %% --- ADDED ---
                % Add text for the number of laps to the top-right corner
                lapString = sprintf('%d laps', numLaps);
                text(0.95, 0.95, lapString, ...
                     'Units', 'normalized', ...
                     'HorizontalAlignment', 'right', ...
                     'VerticalAlignment', 'top', ...
                     'FontSize', 8, ...
                     'Color', [0.3 0.3 0.3]); % Dark gray
            else
                text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
                axis off;
            end
           
            % Format Axes (Standard L-shape)
            set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 10);
            
            % Add Labels only on the Grid Edges
            if thisDay == 1 % Top row: Mouse ID & Injection Site
                if hasInjectionSite && ~isempty(injectionSite)
                    title({mouseID; injectionSite}, 'Interpreter', 'none', 'FontSize', 10);
                else
                    title(mouseID, 'Interpreter', 'none', 'FontSize', 10);
                end
            else
                title('');
            end
            
            % X-axis label and tick labels
            if thisDay == numDays % Bottom row
                if iMouse == 4 % 4th mouse
                    xlabel('Modulation Index');
                end
            else % Not bottom row
                xticklabels({});
            end
            
            % Y-axis label and tick labels
            if iMouse == 1 % First column
                ylabel(['Day ' num2str(dayNum)], 'FontSize', 10);
            end
            % All plots show y-tick labels since Y-axis is NOT scaled
            
        end % end mouse loop
    end % end day loop
    %% Set consistent X-limits and add axis padding for visual gap
    allValues = vertcat(allMetrics{:});
    
    if ~isempty(allValues)
        % Find a robust maximum (e.g., 99th percentile) for X
        globalXMax = prctile(allValues, 99); 
        
        % Calculate 5% padding for the X-axis
        xPadding = (globalXMax - 0) * 0.05;
        if xPadding == 0, xPadding = 0.1; end
        
        % Apply limits AND clean ticks to all axes
        for i = 1:numel(ax)
            if isgraphics(ax(i)) % Check if the plot exists
                % Get this axis's local Y-limit
                yL = get(ax(i), 'YLim');
                yMax = yL(2);
                
                % Calculate local 5% Y padding
                yPadding = (yMax - 0) * 0.05;
                if yPadding == 0 % Handle case of flat or empty plot
                    yPadding = 0.1; 
                end
                % Set LIMITS to be negative (this creates the gap)
                set(ax(i), 'XLim', [-xPadding, globalXMax]);
                set(ax(i), 'YLim', [-yPadding, yMax]);
                
                % Set TICKS to start at 0 (this hides the negative)
                xTicks = get(ax(i), 'XTick');
                xTicks(xTicks < 0) = 0;
                set(ax(i), 'XTick', unique(xTicks));
                
                yTicks = get(ax(i), 'YTick');
                yTicks(yTicks < 0) = 0;
                set(ax(i), 'YTick', unique(yTicks));
            end
        end
    end
    %% Align subplot widths
    drawnow;
    minWidth = Inf;
    validAx = ax(isgraphics(ax)); 
    if ~isempty(validAx)
        lastColIdx = numMice; 
        lastColAx = ax(:, lastColIdx);
        lastColValidAx = lastColAx(isgraphics(lastColAx));
        
        if ~isempty(lastColValidAx)
            pos = get(lastColValidAx(1), 'Position'); 
            minWidth = pos(3);
        else 
             pos = get(validAx(1), 'Position');
             minWidth = pos(3);
        end
        
        if isfinite(minWidth)
            for i = 1:numel(ax)
                if isgraphics(ax(i))
                    pos = get(ax(i), 'Position');
                    pos(3) = minWidth; 
                    set(ax(i), 'Position', pos);
                end
            end
        end
    end
    %% Add Title and Save
    title(tl, sprintf('Modulation Metric Distribution - Signal: %s', signalToUse), ...
           'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
    if ~isempty(savePath)
        try
            disp(['Saving figure to: ' savePath]);
            parentFolder = fileparts(savePath);
            if ~isempty(parentFolder) && ~isfolder(parentFolder)
                mkdir(parentFolder);
            end
            saveas(figHandle, savePath);
        catch ME
            warning('Could not save figure to path: %s', savePath);
            warning('Error: %s', ME.message);
        end
    end
end % end function