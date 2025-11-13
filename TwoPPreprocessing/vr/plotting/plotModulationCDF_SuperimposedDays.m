function figHandle = plotModulationCDF_SuperimposedDays(filteredTable, varargin)
    %% Parse inputs
    p = inputParser;
    addRequired(p, 'filteredTable', @istable);
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    addParameter(p, 'ExcludeMice', {}, @iscellstr);
    addParameter(p, 'DaysToSuperimpose', [1, 3, 5], @isnumeric);
    addParameter(p, 'SavePath', '', @ischar);
    addParameter(p, 'GCaMPInjectionSiteColumn', 'GCaMPInjectionSite', @ischar);
    parse(p, filteredTable, varargin{:});

    signalToUse = p.Results.signalToUse;
    applySmoothing = p.Results.applySmoothing;
    miceToExclude = p.Results.ExcludeMice;
    daysToPlot = p.Results.DaysToSuperimpose;
    numDaysToPlot = length(daysToPlot);
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
    
    if numMice == 0
        warning('filteredTable, after exclusions, has no valid mice.');
        figHandle = [];
        return;
    end
    
    colors = lines(numDaysToPlot);
    figHandle = figure('Position', [100 100 350*numMice 300]);
    tl = tiledlayout(figHandle, 1, numMice, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    
    if applySmoothing
        w = gausswin(10);
        w = w / sum(w);
    end
    
    ax = gobjects(1, numMice);
    allMetrics = cell(numMice, numDaysToPlot); 

    %% Iterate through mice, then days
    for thisMouse = 1:numMice
        mouseID = uniqueMice{thisMouse};
        ax(thisMouse) = nexttile; 
        hold(ax(thisMouse), 'on'); 
        
        hasData = false;
        
        injectionSite = '';
        if hasInjectionSite
            mouseRows = filteredTable(strcmp(filteredTable.MouseID, mouseID), :);
            if ~isempty(mouseRows)
                injectionSite = char(mouseRows.(injectionSiteColumn)(1));
            end
        end
        
        for thisDay = 1:numDaysToPlot
            dayNum = daysToPlot(thisDay);
            thisColor = colors(thisDay, :);
            
            sessionRow = filteredTable( ...
                strcmp(filteredTable.MouseID, mouseID) & ...
                filteredTable.DayOfExperience == dayNum, :);
            
            if isempty(sessionRow)
                continue; 
            end
            
            if height(sessionRow) > 1 
                warning('Multiple entries for %s - Day %d. Using first.', mouseID, dayNum);
                sessionRow = sessionRow(1, :);
            end
            
            sessionString_char = char(sessionRow.Session);
            clear sessionFileInfo response 
            try
                sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID, sessionString_char);
                if ~isfile(sessionFileInfoFilePath), continue; end
                load(sessionFileInfoFilePath, 'sessionFileInfo');
            catch ME
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
            
            if isempty(responseFileIdx), continue; end
            if ~isfield(sessionFileInfo.stimFiles, 'Response') || isempty(sessionFileInfo.stimFiles(responseFileIdx).Response), continue; end
            
            responseFilePath = sessionFileInfo.stimFiles(responseFileIdx).Response;
            try
                if ~isfile(responseFilePath), continue; end
                load(responseFilePath, 'response');
            catch ME
                warning('Error loading response file for %s, %s: %s', mouseID, sessionString_char, ME.message);
                continue;
            end
            
            if ~exist('response', 'var'), continue; end
            if ~isfield(response, 'lapPositionActivity') || ...
               ~isfield(response.lapPositionActivity, signalToUse) || ...
               isempty(response.lapPositionActivity.(signalToUse))
                 continue;
            end
            
            lapActivity = response.lapPositionActivity.(signalToUse);
            if size(lapActivity, 2) < 2, continue; end
            
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
            
            allMetrics{thisMouse, thisDay} = modulationMetric; 
            
            if ~isempty(modulationMetric)
                [f, x] = ecdf(modulationMetric);
                plot(ax(thisMouse), x, f, 'Color', thisColor, 'LineWidth', 2, ...
                    'DisplayName', ['Day ' num2str(dayNum)]);
                hasData = true;
            end
            
        end 
        
        if hasData
            yline(0.5, 'k--', 'LineWidth', 1, 'HandleVisibility','off'); 
            xline(1, 'k--', 'LineWidth', 1, 'HandleVisibility','off'); 
 
            if hasInjectionSite && ~isempty(injectionSite)
                title({mouseID; injectionSite}, 'Interpreter', 'none', 'FontSize', 10);
            else
                title(mouseID, 'Interpreter', 'none', 'FontSize', 10);
            end
            
            if thisMouse == 1
                ylabel('Cumul. Probability');
            else
                yticklabels({});
            end
            
            if thisMouse == 4
                xlabel('Modulation Index');
            end 
            
        else
            title({mouseID; 'No Data'}, 'Interpreter', 'none');
            axis off;
        end
        
    end 
    
    %% Set consistent X-limits and add axis padding for visual gap
    allValues = vertcat(allMetrics{:});
    
    if ~isempty(allValues)
        globalXMax = prctile(allValues, 99); 
        if globalXMax == 0, globalXMax = 1; end 
        
        xPadding = (globalXMax - 0) * 0.05;
        yPadding = (1.0 - 0) * 0.05; 
        
        for i = 1:numel(ax)
            if isgraphics(ax(i)) 
                
                % This is the fix: Apply these settings to ALL axes
                set(ax(i), 'TickDir', 'out', 'box', 'off', 'FontSize', 10);
                
                set(ax(i), 'XLim', [-xPadding, globalXMax]);
                set(ax(i), 'YLim', [-yPadding, 1.0]); 
                
                xTicks = get(ax(i), 'XTick');
                xTicks(xTicks < 0) = 0;
                set(ax(i), 'XTick', unique(xTicks));
                
                set(ax(i), 'YTick', [0, 0.25, 0.5, 0.75, 1.0]);
            end
        end
    end
    
    % Add legend to the last mouse plot
    if isgraphics(ax(numMice))
         legend(ax(numMice), 'Location', 'best', 'Interpreter', 'none');
    end

    %% Add Title and Save
    titleStr = sprintf('Modulation Metric CDF (Days %s) - Signal: %s', ...
        strjoin(string(daysToPlot), ', '), signalToUse);
    title(tl, titleStr, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
    
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
end