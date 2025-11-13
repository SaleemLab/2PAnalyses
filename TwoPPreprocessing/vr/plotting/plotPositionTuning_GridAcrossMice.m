function figHandle = plotPositionTuning_GridAcrossMice(filteredTable, varargin) % Renamed back
% plotPositionTuningGridAcrossMice Generates a grid of position tuning heatmaps.
%
%   figHandle = plotPositionTuningGridAcrossMice(filteredTable, ...)
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

    %% Parse inputs ---
    p = inputParser;
    addRequired(p, 'filteredTable', @istable);
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    addParameter(p, 'ExcludeMice', {'M24043', 'M24046', 'M24048', 'M24049'}, @iscellstr);
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
    figHandle = figure('Position', [100 100 300*numDays 250*numMice]);

    tl = tiledlayout(figHandle, numMice, numDays, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    if applySmoothing
        w = gausswin(6);
        w = w / sum(w);
    end

    % Store axes handles for resizing
    ax = gobjects(numMice, numDays);

    %% Iterate through grid and plot

    for iMouse = 1:numMice
        mouseID = uniqueMice{iMouse};

        injectionSite = '';
        if hasInjectionSite
            mouseRows = filteredTable(strcmp(filteredTable.MouseID, mouseID), :);
            if ~isempty(mouseRows)
                injectionSite = char(mouseRows.(injectionSiteColumn)(1));
            end
        end

        for thisDay = 1:numDays
            dayNum = uniqueDays(thisDay);

            ax(iMouse, thisDay) = nexttile; % Get and store axis handle

            sessionRow = filteredTable( ...
                strcmp(filteredTable.MouseID, mouseID) & ...
                filteredTable.DayOfExperience == dayNum, :);

            if isempty(sessionRow)
                title(''); axis off; continue;
            end
            
            % Multiple sessions issue only arises because of two gains used
            % in the same sessions; TODO: split sessions into two and
            % treat seperately. 
            if height(sessionRow) > 1 
                warning('Multiple entries for %s - Day %d. Using first.', mouseID, dayNum);
                sessionRow = sessionRow(1, :);
            end
            % Removes "" to '' 
            sessionString_char = char(sessionRow.Session);
            % Load sessionFileInfo 
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
            % Find Response File Path 
            stimNames = string({sessionFileInfo.stimFiles.name});
            responseFileIdx = [];
            
            % If two or more runs are combined; select the combined run
            % response file 
            idx_Combined = find(contains(stimNames, "VRCorr") & contains(stimNames, "CombinedRuns"));
            % If two runs exist and have not been combined - select the
            % second one; This again is due to two VR gains used in the
            % same sessions. 
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

            % Check if Response field exists and is not empty 
             if ~isfield(sessionFileInfo.stimFiles, 'Response') || isempty(sessionFileInfo.stimFiles(responseFileIdx).Response)
                 warning('Response field missing or empty for selected stimFile in %s, %s', mouseID, sessionString_char);
                 title(''); axis off; continue;
             end

            responseFilePath = sessionFileInfo.stimFiles(responseFileIdx).Response;
            % Load Response Data 
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

            % Check for valid signal data 
            if ~isfield(response, 'lapPositionActivity') || ...
               ~isfield(response.lapPositionActivity, signalToUse) || ...
               isempty(response.lapPositionActivity.(signalToUse))
                 title(''); axis off;
                 warning('Data Missing in response file for %s, %s', mouseID, sessionString_char);
                 continue;
            end

            lapActivity = response.lapPositionActivity.(signalToUse);
            if size(lapActivity, 2) < 2
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

            % Process Laps (Using MEDIAN) @Aman 
            oddLaps = lapActivity(:, 1:2:end, :);
            evenLaps = lapActivity(:, 2:2:end, :);
            medianOdd = squeeze(median(oddLaps, 2, 'omitnan')); % Changed name
            medianEven = squeeze(median(evenLaps, 2, 'omitnan'));% Changed name

            if size(lapActivity, 1) == 1
                medianOdd = reshape(medianOdd, 1, []);
                medianEven = reshape(medianEven, 1, []);
            end
            normOdd = normalize(medianOdd, 2, 'range');
            normEven = normalize(medianEven, 2, 'range');
            [~, peakIdx] = max(normOdd, [], 2);
            [~, sortIdx] = sort(peakIdx);

            % Plotting 
            imagesc(normEven(sortIdx, :));
            caxis([0 1]); colormap(flipud(gray));

            %  Format Axes 
            set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 10, 'YDir', 'normal');
            xline(50, 'k--', 'LineWidth', 1.5);
            xline(70, 'k--', 'LineWidth', 1.5);
            xline(90, 'k--', 'LineWidth', 1.5);
            xline(110, 'k--', 'LineWidth', 1.5);
            xticks([0 50 70 90 110 140]);
            % Add Labels only on the Grid Edges 
            if iMouse == 1
                title(['Day ' num2str(dayNum)]);
            else
                title('');
            end

            if iMouse == numMice
                xticklabels({'0', '50', '70', '90', '110', '140'});
                xlabel(tl, 'Position (cm)', 'FontSize', 10);
            else
                xticklabels({});
            end

            % Show y-ticks on all, label text only on first column 
            if thisDay == 1
                if ~isempty(injectionSite)
                    ylabel({mouseID; injectionSite}, 'Interpreter', 'none', 'FontSize', 10);
                else
                    ylabel(mouseID, 'Interpreter', 'none', 'FontSize', 10);
                end
            % else % Removed this 'else' block
                % yticklabels({}); % 
            end
            % ---

            if thisDay == numDays
                cb = colorbar;
                ylabel(cb, 'Norm. Activity');
                % Standard ticks/labels for colorbar
                % cb.Ticks = [0, 0.25, 0.5, 0.75, 1];
                % cb.TickLabels = {'', '0.25', '', '0.75', ''};
            end

        end % end day loop
    end % end mouse loop
    %% Align subplot widths @ChatGPT
    drawnow;

    minWidth = Inf;
    % Find minimum width *among plots that have content*
    validAx = ax(isgraphics(ax)); % Get handles of plots that were actually drawn
    if ~isempty(validAx)
        % Check width of the last column plots that were drawn
        lastColIdx = numDays;
        lastColAx = ax(:, lastColIdx);
        lastColValidAx = lastColAx(isgraphics(lastColAx));
        if ~isempty(lastColValidAx)
            pos = get(lastColValidAx(1), 'Position'); % Use first valid axis in last col
            minWidth = pos(3);
        else % Fallback: Use any valid axis if last column is empty
             pos = get(validAx(1), 'Position');
             minWidth = pos(3);
        end


        if isfinite(minWidth)
            for i = 1:numel(ax)
                % Only adjust axes that were actually created
                if isgraphics(ax(i))
                    pos = get(ax(i), 'Position');
                    pos(3) = minWidth; % Set all widths to the minimum width
                    set(ax(i), 'Position', pos);
                end
            end
        end
    end

    %% Add Title and Save @ChatGPT

    % Your code's title said "Even sorted by Even", plotting normOdd
    % Updated title to match plotting normEven sorted by normOdd
    title(tl, sprintf('Even Lap Tuning (Sorted by Odd Laps) - Signal: %s', signalToUse), ...
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