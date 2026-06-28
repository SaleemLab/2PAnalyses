function plotGausTempFit_allROIs(response_1D, pdfFullPath, conditionLabel, useField)
% PLOTSINGLECONDITIONFITS_SPEEDACTIVITY_ALLROIS
%
% Iterates through EVERY single ROI in your dataset, creating a clean 
% 1D speed-tuning visualization (Stationary vs Moving + Fit) per page,
% and appends them sequentially into a single multi-page master PDF.
%
% COLUMNS (Per Page):
%   Left Slot:   Disconnected Stationary Firing Rate at X = 0
%   Right Slot:  Locomotion Tuning Domain (Log-scaled velocity + Fit Line)
%
% USAGE:
%   plotSingleConditionFits_SpeedActivity_allROIs(response, 'Z:\path\to\output.pdf', 'Darkness');

    if nargin < 3 || isempty(conditionLabel), conditionLabel = 'Condition'; end
    if nargin < 4, useField = 'dFFNeuropilCorrected'; end
    
    [targetDir, ~, ~] = fileparts(pdfFullPath);
    if ~isempty(targetDir) && ~exist(targetDir, 'dir')
        mkdir(targetDir);
    end
    
    if exist(pdfFullPath, 'file')
        delete(pdfFullPath);
    end
    
    % Extract base array sizes directly from the 1D profile dimension scale
    totalROIs = numel(response_1D.tuningCurve.(useField).statMean);
    fprintf('=== Multi-Page 1D PDF Engine Initialized ===\nCompiling %d total units for %s. Master file: %s\n\n', totalROIs, conditionLabel, pdfFullPath);
    
    gaussFun = @(params, xdata) params(1) + params(2) .* exp(-(((xdata - params(3)).^2) / (2 * (params(4).^2))));
    
    colorLine = [0.15, 0.15, 0.15];          
    sigGreenBackground  = [0.85, 0.95, 0.85]; % Soft green background for p <= 0.01
    defaultWhiteBackground = [1, 1, 1];        % Plain white background
    
    for targetROI = 1:totalROIs
        fprintf('Appending page for Unit %d/%d...\n', targetROI, totalROIs);
        
        % Centered dashboard layout shape optimized for a standalone clean tuning profile
        figHandle = figure('Name', sprintf('ROI %d PDF Canvas', targetROI), ...
                           'Position', [50, 50, 750, 420], 'Color', 'w', 'Visible', 'off');
        
        vOffset = 0.15;
        rowHeights = 0.70;
        
        % --- EXTRACT 1D SPEED TUNING PROPERTIES ---
        edges   = response_1D.tuningCurve.speedBins;
        centers = edges(1:end-1) + diff(edges)/2;
        y_stat  = response_1D.tuningCurve.(useField).statMean(targetROI);
        y_statE = response_1D.tuningCurve.(useField).statSEM(targetROI);
        y_move  = response_1D.tuningCurve.(useField).moveMean(targetROI, :);
        y_moveE = response_1D.tuningCurve.(useField).moveSEM(targetROI, :);
        
        clsField = 'classification';

        % 
%         if isfield(response_1D.tuningCurve.(useField), 'classification_80_20')
%             clsField = 'classification_80_20';
%         else
%             clsField = 'classification';
%         end
        
        params  = response_1D.tuningCurve.(useField).(clsField).fitParams(targetROI, :);
        r2      = response_1D.tuningCurve.(useField).(clsField).R2(targetROI);
        typeStr = response_1D.tuningCurve.(useField).(clsField).tuningType{targetROI};
        
        pVal = NaN;
        if isfield(response_1D.tuningCurve.(useField), 'pValMoving')
            pVal = response_1D.tuningCurve.(useField).pValMoving(targetROI);
        elseif isfield(response_1D.tuningCurve.(useField), 'pValFull')
            pVal = response_1D.tuningCurve.(useField).pValFull(targetROI);
        end
        if isnan(pVal), pStr = 'NaN'; else, pStr = sprintf('%.4f', pVal); end
        
        valid = ~isnan(y_move); centers = centers(valid); y_move = y_move(valid); y_moveE = y_moveE(valid);
        
        all1DY = [y_stat, y_move];
        yLimits1D = [min(all1DY)*0.85, max(all1DY)*1.15];
        if yLimits1D(1) == yLimits1D(2), yLimits1D = [yLimits1D(1)-0.1, yLimits1D(1)+0.1]; end
        
        % --- RENDERING PANELS ---
        axStat = axes('Position', [0.10, vOffset, 0.05, rowHeights]); hold(axStat, 'on');
        axMove = axes('Position', [0.16, vOffset, 0.74, rowHeights]); hold(axMove, 'on');
        
        % Stationary Dot at X = 0
        errorbar(axStat, 0, y_stat, y_statE, 'ok', 'MarkerSize', 6, 'LineWidth', 1.2, 'MarkerEdgeColor', colorLine);
        xlim(axStat, [-0.5, 0.5]); axStat.XTick = 0; axStat.XTickLabel = {'0'}; ylim(axStat, yLimits1D);
        set(axStat, 'TickDir', 'out', 'Box', 'off'); ylabel(axStat, '\DeltaF/F (Neu)', 'FontWeight', 'bold');
        
        % Continuous Velocity Domain
        errorbar(axMove, centers, y_move, y_moveE, 'ok', 'MarkerFaceColor', colorLine, 'MarkerEdgeColor', colorLine, 'MarkerSize', 5);
        fineX = logspace(log10(min(centers)), log10(max(centers)), 200);
        
        if all(isfinite(params))
            plot(axMove, fineX, gaussFun(params, fineX), 'Color', [0, 0.44, 0.74], 'LineWidth', 2.5);
        end
        set(axMove, 'XScale', 'log', 'YTickLabel', '', 'Box', 'off', 'TickDir', 'out'); axMove.YAxis.Visible = 'off';
        
        axMove.XTick = [2, 5, 10, 20, 30, 40]; 
        axMove.XTickLabel = {'2', '5', '10', '20', '30', '40'}; 
        xlim(axMove, [min(centers)*0.85, max(centers)*1.15]); ylim(axMove, yLimits1D); grid(axMove, 'on');
        xlabel(axMove, 'Speed (cm/s)', 'FontWeight', 'bold');
        
        % Conditional header mapping block
        titleText = sprintf('\\bf%s Speed Profile \\rm(%s)\\newline\\bfR^2 = %.3f | p_{shuffle} = %s', conditionLabel, upper(typeStr), r2, pStr);
        tH = title(axMove, titleText, 'Interpreter', 'tex', 'FontSize', 11);
        if ~isnan(pVal) && pVal <= 0.01
            set(tH, 'BackgroundColor', sigGreenBackground, 'EdgeColor', [0.4, 0.7, 0.4], 'Margin', 3);
        else
            set(tH, 'BackgroundColor', defaultWhiteBackground, 'EdgeColor', 'none');
        end
        
        % Append current slice vector frame into file structure
        if targetROI == 1
            exportgraphics(figHandle, pdfFullPath, 'ContentType', 'vector');
        else
            exportgraphics(figHandle, pdfFullPath, 'ContentType', 'vector', 'Append', true);
        end
        close(figHandle);
    end
    fprintf('\n=== 1D Compilation Completed! ===\nManual updated with green significance tiles at:\n%s\n', pdfFullPath);
end