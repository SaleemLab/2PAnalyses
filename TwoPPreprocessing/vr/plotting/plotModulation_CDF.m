function figHandle = plotModulation_CDF(allData)
% plotModulation_CDF Plots superimposed CDFs of Modulation Index for specified days,
% using hardcoded filters for ROI stability and ratio variance.
    
    % 1. Mouse Exclusion
    excludeMice = {'M24043', 'M24046', 'M24048', 'M24049'}; 
    
    % 2. Days and Type to Plot
    daysToPlot = [1 2 3 4 5];  % NOTE: Changed from [1 3 5] to match the heatmap function
    typeToPlot = 'Boutons';
    
    % 3. Ratio Threshold Filter (Ratio > Threshold)
    ratioThreshold = 10; 
    useRatioFilter = true;
    
    % 4. Hardcoded Save Path (For consistency)
    hardcodedSavePath = 'Z:\ibn-vision\USERS\Sonali\Figures\Heatmaps_SpatialTuningCurves\CDF-HalvesCorrAndRatioVarianceToTuningVariance.png';
    % ----------------------------------------------------------
    
    % The input parser is removed since all filtering parameters are now hardcoded.
    
    uMice = setdiff(unique({allData.MouseID}), excludeMice, 'stable');
    nMice = length(uMice);
    colors = lines(max(daysToPlot)); % Consistent colors for days
    
    % Create figure and tiled layout
    figHandle = figure('Position', [100 100 350*nMice 300]);
    t = tiledlayout(figHandle, 1, nMice, 'TileSpacing', 'compact', 'Padding', 'compact');
    axList = gobjects(0);
    
    for m = 1:nMice
        ax = nexttile; axList(end+1) = ax;
        hold(ax, 'on');
        hasData = false;
        
        for d = daysToPlot
            % 1. Find the specific session based on hardcoded filters
            thisSession = allData(strcmp({allData.MouseID}, uMice{m}) & ...
                                  [allData.Day] == d & ...
                                  strcmpi({allData.Type}, typeToPlot));
                                  
            if isempty(thisSession), continue; end
            if length(thisSession) > 1, thisSession = thisSession(1); end
            
            % --- ROI FILTERING LOGIC ---
            modulationIndex = thisSession.Modulation;
            
            % Start with all ROIs
            roisToKeepIdx = 1:length(modulationIndex);
            
            % FILTERING 1: STABILITY INDEX (Always Active if field exists)
            if isfield(thisSession, 'lapCorr_HalvesStableIdx') && ~isempty(thisSession.lapCorr_HalvesStableIdx)
                stableRoiIndices = thisSession.lapCorr_HalvesStableIdx;
                roisToKeepIdx = intersect(roisToKeepIdx, stableRoiIndices);
            end
            
            % FILTERING 2: RATIO THRESHOLD (Hardcoded to > 40)
            if useRatioFilter
                if isfield(thisSession, 'ratioVarToTuningVar') && ~isempty(thisSession.ratioVarToTuningVar)
                    ratioValues = thisSession.ratioVarToTuningVar;
                    
                    % Comparison: Find indices where the ratio is GREATER THAN the hardcoded threshold (40)
                    ratioIdx = find(ratioValues > ratioThreshold); 
                    
                    roisToKeepIdx = intersect(roisToKeepIdx, ratioIdx);
                    
                    fprintf('  -> CDF Mouse %s Day %d: %d ROIs kept after filtering.\n', ...
                        uMice{m}, d, length(roisToKeepIdx));
                end
            end
            % -------------------------

            if isempty(roisToKeepIdx)
                warning('CDF: No ROIs satisfy ALL active filters for %s Day %d. Skipping.', uMice{m}, d);
                continue;
            end
            
            % --- APPLY FILTER AND PLOT CDF ---
            
            % Select the Modulation Index for the filtered ROIs
            mi = modulationIndex(roisToKeepIdx);
            mi = mi(~isnan(mi)); % Remove NaNs
            
            if ~isempty(mi)
                [f, x] = ecdf(mi);
                plot(ax, x, f, 'Color', colors(d,:), 'LineWidth', 2, ...
                     'DisplayName', sprintf('Day %d', d));
                hasData = true;
            end
        end
        
        % Formatting
        title(ax, uMice{m}, 'Interpreter', 'none');
        set(ax, 'TickDir', 'out', 'box', 'off');
        if m==1, ylabel(ax, 'Cumul. Probability'); end
        if m==round(nMice/2), xlabel(ax, 'Modulation Index'); end
        if hasData
             xline(ax, 1, 'k:', 'HandleVisibility', 'off');
             yline(ax, 0.5, 'k:', 'HandleVisibility', 'off');
        end
    end
    
    % Add legend to last plot
    if isgraphics(ax), legend(ax, 'Location', 'best'); end
    linkaxes(axList(isgraphics(axList)), 'x');
    
    % Global Title reflects the filtering
    filterInfo = sprintf(' (Ratio > %g Filtered)', ratioThreshold);
    title(t, sprintf('Modulation CDF (%s) - Days %s %s', typeToPlot, num2str(daysToPlot), filterInfo), ...
        'FontSize', 14);
        
    % Save using hardcoded path
    saveas(figHandle, hardcodedSavePath);
end