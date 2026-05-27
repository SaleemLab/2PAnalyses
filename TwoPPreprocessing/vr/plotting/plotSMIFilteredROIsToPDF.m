function plotSMIFilteredROIsToPDF(sessionMatrix)
    % This version creates one PDF per session with MouseID in the filename
    
    numSessions = size(sessionMatrix, 1);
    
    for s = 1:numSessions
        % Stability fix: Isolate row to avoid comma-separated list errors
        currRow = sessionMatrix(s, :);
        
        if istable(sessionMatrix)
            cond = currRow.ConditionData{1}.Baseline;
            fveValues = currRow.FVE{1};
            sessionID = char(currRow.Session{1});
            mouseID = char(currRow.MouseID{1}); 
        else
            % FIX: Add (1) indexing to currRow to ensure a single value for struct access
            cond = currRow(1).ConditionData.Baseline;
            fveValues = currRow(1).FVE;
            sessionID = char(currRow(1).Session);
            mouseID = char(currRow(1).MouseID);
        end
        
        % Check for valid SMI results
        if ~isfield(cond, 'SMI'), continue; end
        validIdx = find(~isnan(cond.SMI)); 
        
        if isempty(validIdx), continue; end
        
        % Define path and filename
        filepath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\SMI\ROIS-Heatmaps-SMI';
        
        % Ensure the directory exists
        if ~exist(filepath, 'dir'), mkdir(filepath); end
        
        Name = sprintf('%s_%s_ROI_Analysis.pdf', mouseID, sessionID);
        pdfName = fullfile(filepath, Name); 
        
        if exist(pdfName, 'file'), delete(pdfName); end
        
        lapActivity = cond.LapActivity; 
        numLaps = size(lapActivity, 2);
        bins = 1:size(lapActivity, 3);
        
        fprintf('Generating PDF for Mouse %s, Session %s...\n', mouseID, sessionID);
        
        for j = 1:length(validIdx)
            roiIdx = validIdx(j);
            
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 600]);
            tlo = tiledlayout(2, 1, 'Padding', 'compact');
            
            title(tlo, sprintf('Mouse: %s | Sess: %s | ROI: %d | SMI: %.3f | FVE: %.3f', ...
                mouseID, sessionID, roiIdx, cond.SMI(roiIdx), fveValues(roiIdx)), 'FontWeight', 'bold');
            
            % 1. Heatmap
            nexttile;
            roiHeatmap = squeeze(lapActivity(roiIdx, :, :));
            imagesc(bins, 1:numLaps, roiHeatmap);
            colormap(jet); colorbar; 
            ylabel('Laps'); title('Spatial Activity Heatmap');
            
            % 2. Line Plot
            nexttile;
            meanTuning = mean(roiHeatmap, 1, 'omitnan');
            plot(bins, meanTuning, 'k', 'LineWidth', 2);
            xlabel('Spatial Bins'); ylabel('Avg Activity'); 
            grid on; xlim([min(bins) max(bins)]);
            
            % Save to the session-specific PDF
            exportgraphics(fig, pdfName, 'Append', true, 'ContentType', 'vector');
            close(fig); 
        end
    end
end