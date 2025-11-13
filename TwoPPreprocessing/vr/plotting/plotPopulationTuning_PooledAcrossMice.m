function figHandle = plotPopulationTuning_PooledAcrossMice(filteredTable, TargetArea, varargin)
% plotPopulationTuning Generates a single-row heatmap of concatenated
% population activity from a pre-filtered table.
%
%   figHandle = plotPopulationTuning(filteredTable, TargetArea, ...)
%   Returns a handle to the generated figure.
%
%   Required Argument:
%       filteredTable     - Table *already filtered* for the desired area.
%       TargetArea        - String, the name of the area for plot labels (e.g., 'RSP')
%
%   Additional parameters:
%       'signalToUse'     - String, (default: 'dFFNeuropilCorrected')
%       'applySmoothing'  - Logical, (default: true)
%       'SavePath'        - String, full file path (e.g., 'C:\Fig.png')
%                           If provided, saves the figure. (default: '')

    %% Parse inputs ---
    p = inputParser;
    addRequired(p, 'filteredTable', @istable);
    addRequired(p, 'TargetArea', @ischar);

    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    addParameter(p, 'SavePath', '', @ischar);

    parse(p, filteredTable, TargetArea, varargin{:});
    signalToUse = p.Results.signalToUse;
    applySmoothing = p.Results.applySmoothing;
    savePath = p.Results.SavePath;

    %% Identify grid layout from table

    requiredCols = {'MouseID', 'Session', 'DayOfExperience'};
    if ~all(ismember(requiredCols, filteredTable.Properties.VariableNames))
        error('filteredTable must contain: %s', strjoin(requiredCols, ', '));
    end
    uniqueDays = unique(filteredTable.DayOfExperience, 'stable');
    numDays = length(uniqueDays);
    if numDays == 0
        warning('filteredTable contains no data for the specified days.');
        figHandle = [];
        return;
    end
    figHandle = figure('Position', [100 100 300*numDays 400]);
    tl = tiledlayout(figHandle, 1, numDays, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    if applySmoothing
        w = gausswin(6);
        w = w / sum(w);
    end

    ax = gobjects(1, numDays);

    %% Iterate through days and plot

    for iDay = 1:numDays
        dayNum = uniqueDays(iDay);

        ax(iDay) = nexttile;

        dayRows = filteredTable( ...
            filteredTable.DayOfExperience == dayNum, :);

        if isempty(dayRows)
            title(['Day ' num2str(dayNum)]);
            axis off;
            continue;
        end

        all_medianOdd_cell = {}; 
        all_medianEven_cell = {}; 

        for iMouse = 1:height(dayRows)
            mouseRow = dayRows(iMouse, :);
            mouseID = mouseRow.MouseID{1}; 
            sessionString_char = char(mouseRow.Session);

            %  Load sessionFileInfo 
            clear sessionFileInfo response

            try
                sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID, sessionString_char);
                if ~isfile(sessionFileInfoFilePath)
                    warning('Info File Missing for %s, %s', mouseID, sessionString_char);
                    continue;
                end
                load(sessionFileInfoFilePath, 'sessionFileInfo');
            catch ME
                warning('Error loading sessionFileInfo for %s, %s: %s', mouseID, sessionString_char, ME.message);
                continue;
            end

            %  Find Response
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
                warning('No valid VRCorr file found for %s, %s', mouseID, sessionString_char);
                continue;
            end

            % Check if the selected stimFile entry has the 'Response' field 
            if ~isfield(sessionFileInfo.stimFiles, 'Response') || isempty(sessionFileInfo.stimFiles(responseFileIdx).Response)
                 warning('Response field missing or empty for selected stimFile in %s, %s', mouseID, sessionString_char);
                 continue; % Skip this mouse session if field is missing
            end

            responseFilePath = sessionFileInfo.stimFiles(responseFileIdx).Response;

            %  Load Response Data
            try
                if ~isfile(responseFilePath)
                     warning('Response file path found, but file not found at: %s', responseFilePath);
                     continue;
                end
                load(responseFilePath, 'response');
            catch ME
                warning('Error loading response file for %s, %s: %s', mouseID, sessionString_char, ME.message);
                continue;
            end

            if ~exist('response', 'var')
                warning('File loaded, but ''response'' variable not found in: %s', responseFilePath);
                continue;
            end

            %  Check for valid signal data 
            if ~isfield(response, 'lapPositionActivity') || ...
               ~isfield(response.lapPositionActivity, signalToUse) || ...
               isempty(response.lapPositionActivity.(signalToUse))
                 warning('Data Missing in response file for %s, %s', mouseID, sessionString_char);
                 continue;
            end

            lapActivity_mouse = response.lapPositionActivity.(signalToUse);
            if size(lapActivity_mouse, 2) < 2
                warning('(<2 laps) for %s, %s', mouseID, sessionString_char);
                continue;
            end


            % Apply spatial smoothing
            if applySmoothing
                for iCell = 1:size(lapActivity_mouse, 1)
                    for iLap = 1:size(lapActivity_mouse, 2)
                        trace = squeeze(lapActivity_mouse(iCell, iLap, :));
                        if all(isnan(trace)), continue; end
                        nanMask = isnan(trace);
                        trace(nanMask) = 0;
                        smoothed = filtfilt(w, 1, trace);
                        smoothed(nanMask) = NaN;
                        lapActivity_mouse(iCell, iLap, :) = smoothed;
                    end
                end
            end

            % Process Laps (Using MEDIAN) 
            oddLaps = lapActivity_mouse(:, 1:2:end, :);
            evenLaps = lapActivity_mouse(:, 2:2:end, :);
            medianOdd_mouse = squeeze(median(oddLaps, 2, 'omitnan')); % Changed name
            medianEven_mouse = squeeze(median(evenLaps, 2, 'omitnan'));% Changed name

            if size(lapActivity_mouse, 1) == 1
                medianOdd_mouse = reshape(medianOdd_mouse, 1, []);
                medianEven_mouse = reshape(medianEven_mouse, 1, []);
            end

            % Store the 2D matrices
            all_medianOdd_cell{end+1} = medianOdd_mouse;
            all_medianEven_cell{end+1} = medianEven_mouse;

        end 

        % Check if we collected any data at all for this day 
        if isempty(all_medianOdd_cell)
            title(['Day ' num2str(dayNum)]);
            axis off;
            continue;
        end

        medianOdd = vertcat(all_medianOdd_cell{:}); % Changed name
        medianEven = vertcat(all_medianEven_cell{:});% Changed name

        normOdd = normalize(medianOdd, 2, 'range');
        normEven = normalize(medianEven, 2, 'range');
        [~, peakIdx] = max(normOdd, [], 2);
        [~, sortIdx] = sort(peakIdx);

        imagesc(normEven(sortIdx, :));
        caxis([0 1]); colormap(flipud(gray));

        set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 10, 'YDir', 'normal');
        xline(50, 'k--', 'LineWidth', 1.5);
        xline(70, 'k--', 'LineWidth', 1.5);
        xline(90, 'k--', 'LineWidth', 1.5);
        xline(110, 'k--', 'LineWidth', 1.5);
        xticks([0 50 70 90 110 140]);
        xticklabels({'0', '50', '70', '90', '110', '140'});

        title(['Day ' num2str(dayNum)]);

        if iDay == 1
            ylabel({[TargetArea ' Population']; '(All ROIs Concatenated)'}, ...
                   'Interpreter', 'none', 'FontSize', 9);
        % else % Removed this 'else' block
        %    yticklabels({}); % No longer hide y-tick labels
        end
        % ---

        if iDay == numDays
            cb = colorbar;
            ylabel(cb, 'Norm. Activity');
            %cb.Ticks = [0, 0.25, 0.5, 0.75, 1];
            %cb.TickLabels = {'', '0.25', '', '0.75', ''};
        end

    end 

    %% Align Subplot Widths @Chat
    drawnow;

    minWidth = Inf;
    % Check if last axis exists before getting position
    if isgraphics(ax(numDays))
        pos = get(ax(numDays), 'Position');
        minWidth = pos(3);
    else
        % Fallback if last plot is empty: find last valid plot
        validAx = ax(isgraphics(ax));
        if ~isempty(validAx)
            pos = get(validAx(end), 'Position');
            minWidth = pos(3);
        end
    end


    if isfinite(minWidth) && numDays > 1
        for i = 1:numDays
            % Only align axes that were actually created
            if isgraphics(ax(i))
                pos = get(ax(i), 'Position');
                pos(3) = minWidth;
                set(ax(i), 'Position', pos);
            end
        end
    end

    %% Add title and save

    title(tl, sprintf('%s Population Tuning (Sorted by Odd Laps) - Signal: %s', ...
           TargetArea, signalToUse), ...
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

end 