function [figHandle, dataStats] = computeModulationIndex_AcrossLearningDays(varargin)
% computeModulationIndex Computes modulation across space ((max-min)/mean) 
% across days and mice, with optional plotting.
%
% Usage:
%   [fig, data] = computeModulation_Unified('Plotting', true, 'ImagedType', 'Boutons');
%   [fig, data] = computeModulation_Unified(myFilteredTable, 'SavePath', 'C:\fig.png');
%
% Outputs:
%   figHandle - Handle to the generated figure (empty if Plotting is false).
%   dataStats - Structure containing:
%       .metricsGrid  - Cell array (Days x Mice) of modulation values for every cell.
%       .lapsGrid     - Matrix (Days x Mice) of lap counts per session.
%       .mouseIDList  - Cell array of mouse IDs (columns of the grid).
%       .dayList      - Vector of days (rows of the grid).

    %% 1. Parse Inputs
    p = inputParser;
    p.KeepUnmatched = true; % Allow flexibility if you have other params to pass blindly

    % -- Core Data Inputs --
    % Optional: User can pass a table directly as the first argument
    addOptional(p, 'inputTable', [], @istable);
    
    % -- Filtering Defaults (from your script) --
    addParameter(p, 'DaysOfInterest', [1, 2, 3, 4, 5], @isnumeric);
    addParameter(p, 'ImagedType', 'Boutons', @ischar); 
    addParameter(p, 'ExcludeMice', {'M24046', 'M24043', 'M24048', 'M24049'}, @iscellstr);

    % -- Analysis/Plotting Options --
    addParameter(p, 'Plotting', true, @islogical);
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    addParameter(p, 'SavePath', '', @ischar);
    addParameter(p, 'GCaMPInjectionSiteColumn', 'GCaMPInjectionSite', @ischar);

    parse(p, varargin{:});
    
    % Extract parsed results for easier access
    filteredTable = p.Results.inputTable;
    doPlotting = p.Results.Plotting;
    signalToUse = p.Results.signalToUse;
    applySmoothing = p.Results.applySmoothing;
    miceToExclude = p.Results.ExcludeMice;
    savePath = p.Results.SavePath;
    injectionSiteColumn = p.Results.GCaMPInjectionSiteColumn;

    %% 2. Data Preparation (Filter Table if needed)
    if isempty(filteredTable)
        % If no table was passed as arg 1, we need to build it.
        if exist('filterMasterTable', 'file') ~= 2
             error('filterMasterTable not found in path. Please pass a table directly or ensure this function is available.');
        end
        
        fprintf('Constructing master table internally (Type: %s, Days: %s)...\n', ...
            p.Results.ImagedType, mat2str(p.Results.DaysOfInterest));
            
        filteredTable = filterMasterTable(...
            'Exclude', 0, ...
            'Suite2PPreprocessing', 1, ...
            'DayOfExperience', p.Results.DaysOfInterest, ...
            'ImagedType', p.Results.ImagedType);
            
        if isempty(filteredTable)
             warning('Internal filtering resulted in an empty table.');
             figHandle = []; dataStats = []; return;
        end
    end

    %% 3. Identify Grid Layout
    requiredCols = {'MouseID', 'Session', 'DayOfExperience'};
    if ~all(ismember(requiredCols, filteredTable.Properties.VariableNames))
        error('Table missing required columns: %s', strjoin(requiredCols, ', '));
    end

    uniqueMice = unique(filteredTable.MouseID, 'stable');
    uniqueMice = setdiff(uniqueMice, miceToExclude, 'stable');
    numMice = length(uniqueMice);
    
    uniqueDays = unique(filteredTable.DayOfExperience, 'stable');
    % If user specified days and they are in the table, stick to that order/subset
    if ~isempty(p.Results.DaysOfInterest)
         % Intersect keeps them if they exist, but maybe we want to force the user list?
         % Let's just use what's actually IN the table after filtering.
         uniqueDays = intersect(p.Results.DaysOfInterest, uniqueDays, 'stable');
    end
    numDays = length(uniqueDays);

    if numMice == 0 || numDays == 0
        warning('No valid mice/days found after exclusions.');
        figHandle = []; dataStats = []; return;
    end

    % Initialize Data Output Structures
    dataStats.metricsGrid = cell(numDays, numMice);
    dataStats.lapsGrid = nan(numDays, numMice);
    dataStats.mouseIDList = uniqueMice;
    dataStats.dayList = uniqueDays;

    %% 4. Initialize Figure (if plotting)
    figHandle = [];
    ax = gobjects(numDays, numMice);
    if doPlotting
        % Scale figure size based on grid
        figHandle = figure('Position', [100 100 300*numMice 250*numDays], ...
            'Name', ['Modulation Dist: ' signalToUse]);
        tl = tiledlayout(figHandle, numDays, numMice, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Smoothing window (pre-calculated)
        if applySmoothing
            w = gausswin(6); w = w / sum(w);
        end
    end

    hasInjectionSite = ismember(injectionSiteColumn, filteredTable.Properties.VariableNames);

    %% 5. Main Loop (Days x Mice)
    fprintf('Processing %d Days across %d Mice...\n', numDays, numMice);
    
    for thisDay = 1:numDays
        dayNum = uniqueDays(thisDay);

        for iMouse = 1:numMice
            mouseID = uniqueMice{iMouse};

            % --- A. Retrieve Session Data ---
            sessionRow = filteredTable(strcmp(filteredTable.MouseID, mouseID) & ...
                                       filteredTable.DayOfExperience == dayNum, :);

            if doPlotting
                ax(thisDay, iMouse) = nexttile(tl);
            end

            if isempty(sessionRow)
                if doPlotting, axis off; end
                continue;
            end
            
            if height(sessionRow) > 1
                 sessionRow = sessionRow(1,:); % Take first if duplicate
            end

            % --- B. Load Files ---
            sessionString = char(sessionRow.Session);
            % (Helper function to safely load data and return empty if failed)
            [response, errorMsg] = loadSessionData(mouseID, sessionString, signalToUse);
            
            if isempty(response)
                if doPlotting
                    text(0.5, 0.5, errorMsg, 'HorizontalAlignment', 'center', 'Interpreter', 'none', 'FontSize', 8);
                    axis off; 
                end
                continue;
            end

            lapActivity = response.lapPositionActivity.(signalToUse);
            numLaps = size(lapActivity, 2);
            dataStats.lapsGrid(thisDay, iMouse) = numLaps; % Store lap count

            if numLaps < 2
                if doPlotting
                    text(0.5, 0.5, '< 2 Laps', 'HorizontalAlignment', 'center'); axis off;
                end
                continue;
            end

            % --- C. Compute Metric ---
            % 1. Smooth if requested
            if applySmoothing
                 lapActivity = smoothLapActivity(lapActivity, w);
            end
            
            % 2. Calculate (Max - Min) / Mean
            % Average across laps first to get median tuning curve
            medianTuning = squeeze(mean(lapActivity, 2, 'omitnan'));
            % Handle 1-cell case to ensure correct orientation
            if isvector(medianTuning), medianTuning = medianTuning(:)'; end

            cellMax = max(medianTuning, [], 2, 'omitnan');
            cellMin = min(medianTuning, [], 2, 'omitnan');
            cellMean = mean(medianTuning, 2, 'omitnan');
            
            modMetric = (cellMax - cellMin) ./ cellMean;
            
            % Clean NaNs/Infs
            modMetric(cellMean == 0 | isnan(cellMean) | isinf(modMetric)) = NaN;
            modMetric = modMetric(~isnan(modMetric));
            
            % Store data for output
            dataStats.metricsGrid{thisDay, iMouse} = modMetric;

            % --- D. Plotting (if enabled) ---
            if doPlotting && ~isempty(modMetric)
                histogram(ax(thisDay, iMouse), modMetric, 20, ...
                    'Normalization', 'count', 'EdgeColor', 'none', 'FaceColor', 'k');
                
                % Add Lap Count Label
                text(ax(thisDay, iMouse), 0.95, 0.95, sprintf('%d laps', numLaps), ...
                     'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
                 
                formatAxis(ax(thisDay, iMouse), thisDay, iMouse, numDays, dayNum, ...
                    mouseID, hasInjectionSite, injectionSiteColumn, filteredTable);
            elseif doPlotting
                 text(0.5, 0.5, 'No Valid Cells', 'HorizontalAlignment', 'center'); axis off;
            end

        end % mouse loop
    end % day loop

    %% 6. Finalize Plots
    if doPlotting
        finalizeFigureLayout(ax, dataStats.metricsGrid, numDays, numMice);
        title(tl, sprintf('Modulation ((Max-Min)/Mean) - %s', signalToUse), 'Interpreter', 'none');

        if ~isempty(savePath)
            saveFigure(figHandle, savePath);
        end
    end
    fprintf('Done.\n');
end

%% --- Local Helper Functions ---

function [response, msg] = loadSessionData(mouseID, sessionStr, signalReq)
    response = []; msg = '';
    try
        infoPath = findSessionFileInfoFilePath(mouseID, sessionStr);
        if ~isfile(infoPath), msg = 'Info Missing'; return; end
        loadedInfo = load(infoPath, 'sessionFileInfo');
        sfi = loadedInfo.sessionFileInfo;

        stimNames = string({sfi.stimFiles.name});
        idx = find(contains(stimNames, "VRCorr") & contains(stimNames, "CombinedRuns"), 1);
        if isempty(idx)
             idx = find(contains(stimNames, "VRCorr") & ~contains(stimNames, "CombinedRuns"), 1, 'last');
        end

        if isempty(idx) || ~isfield(sfi.stimFiles(idx), 'Response') || isempty(sfi.stimFiles(idx).Response)
             msg = 'No VRCorr/Response'; return;
        end

        respPath = sfi.stimFiles(idx).Response;
        if ~isfile(respPath), msg = 'Resp File Missing'; return; end
        
        loadedResp = load(respPath, 'response');
        if ~isfield(loadedResp, 'response'), msg = 'Var Missing'; return; end
        
        if ~isfield(loadedResp.response, 'lapPositionActivity') || ...
           ~isfield(loadedResp.response.lapPositionActivity, signalReq)
             msg = 'Signal Missing'; return;
        end
        response = loadedResp.response;
    catch
        msg = 'Load Error';
    end
end

function smoothed = smoothLapActivity(activity, w)
    smoothed = activity;
    for c = 1:size(activity, 1)
        for l = 1:size(activity, 2)
            trace = squeeze(activity(c, l, :));
            mask = isnan(trace);
            if all(mask), continue; end
            trace(mask) = 0;
            filt = filtfilt(w, 1, trace);
            filt(mask) = NaN;
            smoothed(c, l, :) = filt;
        end
    end
end

function formatAxis(axH, dayIdx, mouseIdx, nDays, dayVal, mID, hasInj, injCol, table)
    set(axH, 'TickDir', 'out', 'box', 'off', 'FontSize', 10);
    % Top Row Labels
    if dayIdx == 1
        titleStr = mID;
        if hasInj
             mRow = table(strcmp(table.MouseID, mID), :);
             if ~isempty(mRow), titleStr = {mID; char(mRow.(injCol)(1))}; end
        end
        title(axH, titleStr, 'Interpreter', 'none', 'FontSize', 10);
    else
        title(axH, '');
    end
    % X Labels (bottom row only, 4th col preferred for main label)
    if dayIdx ~= nDays
        xticklabels(axH, {});
    elseif mouseIdx == min(4, size(table,2)) % try to put it near middle
        xlabel(axH, 'Modulation Index');
    end
    % Y Labels (first col only)
    if mouseIdx == 1
        ylabel(axH, ['Day ' num2str(dayVal)], 'FontSize', 10);
    end
end

function finalizeFigureLayout(ax, metricsGrid, nDays, nMice)
    allVals = vertcat(metricsGrid{:});
    if isempty(allVals), return; end
    
    % Robust global X-limit
    xMax = prctile(allVals, 99);
    xPad = max(0.1, xMax * 0.05);
    
    validAx = ax(isgraphics(ax));
    if isempty(validAx), return; end
    
    % Apply uniform limits and clean 0-ticks
    for i = 1:numel(validAx)
        yL = get(validAx(i), 'YLim');
        yPad = max(0.1, yL(2) * 0.05);
        set(validAx(i), 'XLim', [-xPad, xMax], 'YLim', [-yPad, yL(2)]);
        
        % Hide negative ticks created by padding
        xt = get(validAx(i), 'XTick'); set(validAx(i), 'XTick', xt(xt>=0));
        yt = get(validAx(i), 'YTick'); set(validAx(i), 'YTick', yt(yt>=0));
    end
    
    % Align widths (basic attempt)
    drawnow;
    posAll = get(validAx, 'Position');
    if iscell(posAll), posMat = vertcat(posAll{:}); else, posMat = posAll; end
    minW = min(posMat(:,3));
    for i = 1:numel(validAx)
        p = get(validAx(i), 'Position'); p(3) = minW; set(validAx(i), 'Position', p);
    end
end

function saveFigure(h, pathStr)
    try
        folder = fileparts(pathStr);
        if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
        saveas(h, pathStr);
        fprintf('Figure saved: %s\n', pathStr);
    catch ME
        warning('Failed to save figure: %s', ME.message);
    end
end