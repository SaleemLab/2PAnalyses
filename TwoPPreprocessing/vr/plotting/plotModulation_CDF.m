function figHandle = plotModulation_CDF(allData)
% plotModulation_CDF Plots CDFs of Range (Max-Min) for Odd-Even Stable ROIs.
    
    % 1. Settings
    excludeMice = {'M24043', 'M24046', 'M24048', 'M24049'}; 
    daysToPlot = [1 3 5]; 
    typeToPlot = 'Boutons';
    hardcodedSavePath = 'Z:\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\CDF_OddEven_Stability.png';
    
    uMice = setdiff(unique({allData.MouseID}), excludeMice, 'stable');
    nMice = length(uMice);
    colors = lines(max(daysToPlot)); 
    
    figHandle = figure('Position', [100 100 350*nMice 350], 'Color', 'w');
    t = tiledlayout(figHandle, 1, nMice, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for m = 1:nMice
        ax = nexttile; hold(ax, 'on');
        
        for d = daysToPlot
            % Match Mouse, Day, and Type
            thisSession = allData(strcmp({allData.MouseID}, uMice{m}) & ...
                                  [allData.Day] == d & ...
                                  strcmpi({allData.Type}, typeToPlot));
                                  
            if isempty(thisSession), continue; end
            if length(thisSession) > 1, thisSession = thisSession(1); end
            
            % --- CORRECTED FIELD NAME: MeanTuning ---
            tc = thisSession.MeanTuning; 
            
            % Calculate Range (Max-Min) for z-scored dFF
            % tc is [ROIs x Bins], so we operate on dimension 2
            rangeIndex = max(tc, [], 2) - min(tc, [], 2);
            
            % --- ODD-EVEN STABILITY FILTER ---
            if isfield(thisSession, 'lapCorr_OddEvenStableIdx') && ~isempty(thisSession.lapCorr_OddEvenStableIdx)
                stableIdx = thisSession.lapCorr_OddEvenStableIdx;
                
                % Handle logical vs numeric indices automatically
                if islogical(stableIdx), stableIdx = find(stableIdx); end
                
                % Select only the stable ROIs
                mi = rangeIndex(stableIdx);
                mi = mi(~isnan(mi) & ~isinf(mi)); 
                
                if ~isempty(mi)
                    [f, x] = ecdf(mi);
                    % Legend shows Day and the number of Stable ROIs (n)
                    plot(ax, x, f, 'Color', colors(d,:), 'LineWidth', 2, ...
                         'DisplayName', sprintf('Day %d (n=%d)', d, length(mi)));
                end
            end
        end
        
        % Formatting
        title(ax, uMice{m}, 'Interpreter', 'none');
        set(ax, 'TickDir', 'out', 'box', 'off');
        xlim(ax, [0 5]); 
        if m==1, ylabel(ax, 'Cumul. Probability'); end
        if m==round(nMice/2), xlabel(ax, 'Spatial Range (Max-Min Z-score)'); end
        
        % Reference line for a clear spatial field
        xline(ax, 5, 'k:', 'HandleVisibility', 'off'); 
        legend(ax, 'Location', 'southeast', 'FontSize', 8);
    end
    
    title(t, sprintf('Modulation of Odd-Even Stable %s', typeToPlot), 'FontSize', 14);
    
    % Save logic
    if ~exist(fileparts(hardcodedSavePath), 'dir'), mkdir(fileparts(hardcodedSavePath)); end
    saveas(figHandle, hardcodedSavePath);
end