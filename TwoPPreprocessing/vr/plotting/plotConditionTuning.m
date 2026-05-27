function fig = plotConditionTuning(allData, targetArea, days, savePath)
    % plotConditionTuning: Individual plot normalization, fixed Baseline-Odd sorting.
    
    % 1. Identify manipulation conditions
    allConds = {};
    for i = 1:length(allData)
        if isfield(allData(i), 'ConditionData') && ~isempty(allData(i).ConditionData)
            names = fieldnames(allData(i).ConditionData);
            manips = names(~contains(lower(names), 'baseline') & ~contains(lower(names), 'norm'));
            allConds = [allConds, manips'];
        end
    end
    manipConds = unique(allConds, 'stable');
    plotLabels = [{'Baseline Odd', 'Baseline Even'}, manipConds];
    numCols = length(plotLabels);
    numDays = length(days);
    
    % Colors
    omitColor = [0.9 0.2 0.2]; % Red
    swapColor = [0.2 0.4 0.9]; % Blue
    baseColor = [0 0 0];       % Black
    
    fig = figure('Color', 'w', 'Position', [50 100 300*numCols + 120 400*numDays]);
    t = tiledlayout(numDays, numCols, 'TileSpacing', 'compact', 'Padding', 'normal');
    
    for d = 1:numDays
        currentDay = days(d);
        targetDayIdx = [allData.Day] == currentDay;
        
        areaIdx = strcmpi(string({allData.TargetArea_ROI}), string(targetArea));
        sessList = allData(targetDayIdx & areaIdx);
        if isempty(sessList), continue; end
        
        dayData = cell(1, numCols);
        for s = 1:length(sessList)
            cData = sessList(s).ConditionData;
            cNames = fieldnames(cData);
            baseIdx = find(contains(lower(cNames), 'baseline') | contains(lower(cNames), 'norm'), 1);
            if isempty(baseIdx), continue; end
            baseName = cNames{baseIdx};
            
            % Baseline Odd/Even
            fullBase = cData.(baseName).LapActivity(:, :, :);
            dayData{1} = [dayData{1}; squeeze(mean(fullBase(:, 1:2:end, :), 2, 'omitnan'))];
            dayData{2} = [dayData{2}; squeeze(mean(fullBase(:, 2:2:end, :), 2, 'omitnan'))];
            
            % Manipulations
            for m = 1:length(manipConds)
                mName = manipConds{m};
                if isfield(cData, mName)
                    dayData{m+2} = [dayData{m+2}; squeeze(mean(cData.(mName).LapActivity(:, :, :), 2, 'omitnan'))];
                end
            end
        end
        
        % --- SORTING ANCHOR (Baseline Odd) ---
        if isempty(dayData{1}), continue; end
        baseOdd = dayData{1};
        % Self-normalize just the reference for sorting purposes
        refNorm = (baseOdd - min(baseOdd,[],2,'omitnan')) ./ (max(baseOdd,[],2,'omitnan') - min(baseOdd,[],2,'omitnan'));
        [~, peakIdx] = max(refNorm, [], 2);
        [~, sortIdx] = sort(peakIdx);
        
        % 
        maxSomas = 500;
        totalSomas = length(sortIdx);
        if totalSomas > maxSomas
            downsampleIdx = round(linspace(1, totalSomas, maxSomas));
            sortIdx = sortIdx(downsampleIdx);
        end
        
        % --- PLOTTING ---
        for col = 1:numCols
            ax = nexttile();
            currRawData = dayData{col};
            currLabel = lower(plotLabels{col});
            
            if contains(currLabel, 'omit')
                activeColor = omitColor;
            elseif contains(currLabel, 'swap')
                activeColor = swapColor;
            else
                activeColor = baseColor;
            end
            
            if ~isempty(currRawData)
                % SELF-NORMALIZATION: Each plot scaled to its own 0-1 range
                cMin = min(currRawData, [], 2, 'omitnan');
                cMax = max(currRawData, [], 2, 'omitnan');
                cRange = cMax - cMin;
                cRange(cRange == 0) = 1;
                normData = (currRawData - cMin) ./ cRange;
                
                % Slice data with our capped/downsampled sort index
                imagesc(normData(sortIdx, :));
                colormap(ax, flipud(gray)); clim([0 1]);
                set(ax, 'CLim', [0.25 0.75], 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 11);
                
                % --- REMOVE Y AXIS TICKS ---
                yticks(ax, []);
                
                % --- Scale Bar (ONLY ON PLOT 1, Physically Proportional) ---
                numN = size(sortIdx, 1); 
                if col == 1
                    hold(ax, 'on');
                    
                    % Calculate visual rows 100 true somas occupy in this plot
                    if totalSomas > maxSomas
                        visualBarLength = 500 * (maxSomas / totalSomas);
                    else
                        visualBarLength = 500; 
                    end
                    
                    % Anchored high near 75% of visual height
                    y_start = numN * 0.75; 
                    y_end   = y_start + visualBarLength; 
                    y_mid   = (y_start + y_end) / 2;  
                    
                    x_line = -4; 
                    x_text = -12; 
                    
                    if totalSomas >= 100
                        line(ax, [x_line x_line], [y_start y_end], ...
                            'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
                        
                        text(ax, x_text, y_mid, '500 somas', ...
                            'Rotation', 90, ...
                            'HorizontalAlignment', 'center', ... 
                            'VerticalAlignment', 'middle', ...   
                            'FontSize', 12, ...
                            'FontWeight', 'normal', ...
                            'Clipping', 'off');
                    end
                end
                
                % Title
                if d == 1
                    tH = title(plotLabels{col}, 'Interpreter', 'none');
                    tH.Color = activeColor;
                end
                
                % Landmark Ticks and Lines
                landmarks = [40, 80, 120, 160];
                xticks(landmarks);
                xticklabels({'40', '80', '120', '160'});
                
                for l = 1:4
                    isHighlighted = contains(currLabel, num2str(l));
                    lColor = [0 0 0]; lWidth = 0.5;
                    if isHighlighted
                        lColor = activeColor; lWidth = 1.8;
                    end
                    xline(landmarks(l), '--', 'Color', lColor, 'LineWidth', lWidth);
                end
            end
            
            % if col == 1
            %     ylabel(sprintf('Day %d', currentDay), 'FontWeight', 'bold');
            % end
            
            % --- INTEGRATED CUSTOM COLORBAR BLOCK ---
            if d == numDays && col == numCols
                cb = colorbar(ax);
                
                % 1. Force ticks to the right
                cb.Location = 'eastoutside'; 
                
                % 2. Manual position override [left bottom width height]
                pos = get(ax, 'Position');
                cb.Position = [pos(1)+pos(3)+0.03, pos(2)+pos(4)*0.3, 0.012, pos(4)*0.4];
                
                % 3. Ticks and Styling
                cb.Ticks = [0.25 0.75];
                cb.TickLabels = {'0.25', '0.75'};
                cb.TickDirection = 'out';
                cb.Box = 'off';
                
                % 4. FORCE LABEL INSIDE
                cb.Label.String = 'Activity (norm.)';
                cb.Label.FontSize = 10;
                cb.Label.Rotation = 90;
                
                % MANUALLY NUDGE THE LABEL POSITION
                cb.Label.Units = 'normalized';
                cb.Label.Position = [-0.8, 0.5, 0]; 
                
                % Ensure alignment is centered so it doesn't drift
                cb.Label.HorizontalAlignment = 'center';
                cb.Label.VerticalAlignment = 'bottom';
            end
        end
    end
    
    xlabel(t, 'Position (cm)', 'FontSize', 12);
    
    if nargin > 3 && ~isempty(savePath)
        [fDir, fName, ~] = fileparts(savePath);
        baseCleanPath = fullfile(fDir, fName);
        
    
        exportgraphics(fig, [baseCleanPath, '.png'], 'Resolution', 600);
        

        exportgraphics(fig, [baseCleanPath, '.svg'], 'ContentType', 'vector');
    end
end
% function fig = plotConditionTuning(allData, targetArea, days, savePath)
%     % plotConditionTuning: Individual plot normalization, fixed Baseline-Odd sorting.
% 
%     % 1. Identify manipulation conditions
%     allConds = {};
%     for i = 1:length(allData)
%         if isfield(allData(i), 'ConditionData') && ~isempty(allData(i).ConditionData)
%             names = fieldnames(allData(i).ConditionData);
%             manips = names(~contains(lower(names), 'baseline') & ~contains(lower(names), 'norm'));
%             allConds = [allConds, manips'];
%         end
%     end
%     manipConds = unique(allConds, 'stable');
%     plotLabels = [{'Baseline Odd', 'Baseline Even'}, manipConds];
%     numCols = length(plotLabels);
%     numDays = length(days);
% 
%     % Colors
%     omitColor = [0.9 0.2 0.2]; % Red
%     swapColor = [0.2 0.4 0.9]; % Blue
%     baseColor = [0 0 0];       % Black
% 
%     fig = figure('Color', 'w', 'Position', [50 100 300*numCols + 120 400*numDays]);
%     t = tiledlayout(numDays, numCols, 'TileSpacing', 'none', 'Padding', 'compact');
% 
%     for d = 1:numDays
%         currentDay = days(d);
%         targetDayIdx = [allData.Day] == currentDay;
% 
%         areaIdx = strcmpi(string({allData.TargetArea_ROI}), string(targetArea));
%         sessList = allData(targetDayIdx & areaIdx);
%         if isempty(sessList), continue; end
% 
%         dayData = cell(1, numCols);
%         for s = 1:length(sessList)
%             cData = sessList(s).ConditionData;
%             cNames = fieldnames(cData);
%             baseIdx = find(contains(lower(cNames), 'baseline') | contains(lower(cNames), 'norm'), 1);
%             if isempty(baseIdx), continue; end
%             baseName = cNames{baseIdx};
% 
%             % ROI selection
%             % if isfield(sessList(s), 'highlyCorrBoutons') && ~isempty(sessList(s).highlyCorrBoutons)
%             %     idx = sessList(s).highlyCorrBoutons;
%             % else
%             %     idx = 1:sessList(s).NumCells;
%             % end
% 
%             % Baseline Odd/Even
%             fullBase = cData.(baseName).LapActivity(:, :, :);
%             dayData{1} = [dayData{1}; squeeze(mean(fullBase(:, 1:2:end, :), 2, 'omitnan'))];
%             dayData{2} = [dayData{2}; squeeze(mean(fullBase(:, 2:2:end, :), 2, 'omitnan'))];
% 
%             % Manipulations
%             for m = 1:length(manipConds)
%                 mName = manipConds{m};
%                 if isfield(cData, mName)
%                     dayData{m+2} = [dayData{m+2}; squeeze(mean(cData.(mName).LapActivity(:, :, :), 2, 'omitnan'))];
%                 end
%             end
%         end
% 
%         % --- SORTING ANCHOR (Baseline Odd) ---
%         if isempty(dayData{1}), continue; end
%         baseOdd = dayData{1};
%         % Self-normalize just the reference for sorting purposes
%         refNorm = (baseOdd - min(baseOdd,[],2,'omitnan')) ./ (max(baseOdd,[],2,'omitnan') - min(baseOdd,[],2,'omitnan'));
%         [~, peakIdx] = max(refNorm, [], 2);
%         [~, sortIdx] = sort(peakIdx);
% 
%         % --- PLOTTING ---
%         for col = 1:numCols
%             ax = nexttile();
%             currRawData = dayData{col};
%             currLabel = lower(plotLabels{col});
% 
%             if contains(currLabel, 'omit')
%                 activeColor = omitColor;
%             elseif contains(currLabel, 'swap')
%                 activeColor = swapColor;
%             else
%                 activeColor = baseColor;
%             end
% 
%             if ~isempty(currRawData)
%                 % SELF-NORMALIZATION: Each plot scaled to its own 0-1 range
%                 cMin = min(currRawData, [], 2, 'omitnan');
%                 cMax = max(currRawData, [], 2, 'omitnan');
%                 cRange = cMax - cMin;
%                 cRange(cRange == 0) = 1;
%                 normData = (currRawData - cMin) ./ cRange;
% 
%                 imagesc(normData(sortIdx, :));
%                 colormap(ax, flipud(gray)); clim([0 1]);
%                 set(ax, 'CLim', [0.25 0.75], 'YDir', 'normal', 'TickDir', 'out', 'FontSize', 10);
% 
% 
%                 % Title
%                 if d == 1
%                     tH = title(plotLabels{col}, 'Interpreter', 'none');
%                     tH.Color = activeColor;
%                 end
% 
%                 % Landmark Ticks and Lines
%                 landmarks = [40, 80, 120, 160];
%                 xticks(landmarks);
%                 xticklabels({'40', '80', '120', '160'});
% 
%                 for l = 1:4
%                     isHighlighted = contains(currLabel, num2str(l));
%                     lColor = [0 0 0]; lWidth = 0.5;
%                     if isHighlighted
%                         lColor = activeColor; lWidth = 1.8;
%                     end
%                     xline(landmarks(l), '--', 'Color', lColor, 'LineWidth', lWidth);
%                 end
%             end
% 
%             if col == 1
%                 ylabel(sprintf('Day %d\n(n=%d Somas)', currentDay, size(dayData{1},1)), 'FontWeight', 'bold');
%             else
%                 set(ax, 'YTickLabel', []);
%             end
% 
%             if col == numCols
%                 cb = colorbar;
%                 ylabel(cb, '\DeltaF/F (norm.)');
%             end
%         end
%     end
% 
%     xlabel(t, 'Position (cm)', 'FontSize', 12, 'FontWeight', 'bold');
% 
%     if nargin > 3 && ~isempty(savePath)
%         exportgraphics(fig, savePath, 'Resolution', 300);
%     end
% end
% 
% % interagte this 
% numN = size(normEvenSorted, 1);
%             if numN >= 100
%                 hold(ax, 'on');
% 
%                 % Anchor the bottom of the bar to 20% of the plot height
%                 % This ensures a consistent visual "starting height" across panels
%                 y_start = numN * 0.2; 
%                 y_end   = y_start + 100; % Absolute 100 units long
%                 y_mid   = y_start + 50;  % The exact center of the bar
% 
%                 % Horizontal spacing (snug to the plot)
%                 x_line = -6; 
%                 x_text = -14; 
% 
%                 % 1. Draw the Line
%                 line(ax, [x_line x_line], [y_start y_end], ...
%                     'Color', 'k', 'LineWidth', 2, 'Clipping', 'off');
% 
%                 % 2. Draw the Text (Centered on the line)
%                 text(ax, x_text, y_mid, '100 somas', ...
%                     'Rotation', 90, ...
%                     'HorizontalAlignment', 'center', ... % Centers the text vertically
%                     'VerticalAlignment', 'middle', ...   % Centers the text horizontally relative to x_text
%                     'FontSize', 11, ...
%                     'FontWeight', 'normal', ...
%                     'Clipping', 'off');
%             end