function plotDistributionsROISelectionCriteria(allData)
    % Define the target directory
    saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\DistibutionsAllCritera';
    
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    metrics = {'realRangePercentileRank', 'realPeakPercentileRank', ...
               'ratioVarToTuningVar', 'ratioVarToTuningRange', ...
               'lapCorr_HalvesRho', 'lapCorr_OddEvenRho'};
    uniqueMice = unique({allData.MouseID});
    
    fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
   
    for m = 1:length(metrics)
        metricName = metrics{m};
        pdfFileName = fullfile(saveDir, sprintf('Distributions_%s.png', metricName));
        if exist(pdfFileName, 'file'), delete(pdfFileName); end
        
        for i = 1:length(uniqueMice)
            mouseID = uniqueMice{i};
            mouseIdx = strcmp({allData.MouseID}, mouseID);
            mouseData = allData(mouseIdx);
            numSessions = length(mouseData);
            
            if numSessions == 0, continue; end
            clf(fig);
            
            numCols = 2;
            numRows = ceil(numSessions / numCols);
            t = tiledlayout(numRows, numCols, 'Padding', 'loose', 'TileSpacing', 'compact');
            sgtitle(t, sprintf('Mouse: %s | Metric: %s', mouseID, metricName), 'FontSize', 14, 'FontWeight', 'bold');
            
            for s = 1:numSessions
                ax = nexttile; % Capture the specific axes handle for this tile
                hold(ax, 'on');
                
                % FORCED BOX OFF: This uses the specific handle to ensure it hits the right axes
                ax.Box = 'off'; 
                
                try
                    if isfield(mouseData(s), metricName) && ~isempty(mouseData(s).(metricName))
                        currentData = mouseData(s).(metricName);
                        histogram(ax, currentData); % Explicitly plot into 'ax'
                        
                        medVal = nanmedian(currentData);
                        xl = xline(ax, medVal, '--r', 'LineWidth', 1.5);
                        xl.Label = sprintf('Med: %.2f', medVal);
                        
                        if contains(metricName, 'Percentile')
                            xlim(ax, [0 100]);
                        elseif contains(metricName, 'Rho')
                            xlim(ax, [-1 1]);
                        end
                    else
                        text(ax, 0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
                    end
                catch ME
                    text(ax, 0.5, 0.5, 'Plotting Error', 'HorizontalAlignment', 'center', 'Color', 'r');
                end
                
                subTitleStr = sprintf('Sess: %s | Type: %s\nROIs: %d | Laps: %d', ...
                    string(mouseData(s).Session), ...
                    string(mouseData(s).TypeImaged_ROI), ...
                    mouseData(s).NumCells, ...
                    mouseData(s).NumLaps);
                title(ax, subTitleStr, 'FontSize', 9);
            end
            
            exportgraphics(fig, pdfFileName, 'Append', true);
        end
    end
    close(fig);
end