function plotRunningTuningCurves(respGray, respDark, pdfPath)
% pdfPath: Full path ending in .pdf
% respDark: Can be empty [] if only plotting Gray Screen

% 1. Setup Path and Extension Safety
[folder, name, ext] = fileparts(pdfPath);
if ~strcmpi(ext, '.pdf')
    pdfPath = fullfile(folder, [name '.pdf']);
end

% Check if Darkness data was provided
hasDark = ~isempty(respDark) && isfield(respDark, 'tuningCurve');

% 2. Identify Data Fields
fieldsGray = fieldnames(respGray.tuningCurve);
targetFields = {'dFF', 'spks', 'dFFNeuropilCorrected', 'spikes', 'spk'};

if hasDark
    fieldsDark = fieldnames(respDark.tuningCurve);
    plotFields = intersect(intersect(fieldsGray, fieldsDark), targetFields, 'stable');
else
    plotFields = intersect(fieldsGray, targetFields, 'stable');
end

numROIs = size(respGray.tuningCurve.(plotFields{1}).statMean, 1);
numFields = length(plotFields);

% 3. Create Figure
fig = figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.9], 'Visible', 'off', 'Color', 'w');

fprintf('Generating PDF for %d ROIs (Comparison: %d)...\n', numROIs, hasDark);

for iROI = 1:numROIs
    clf(fig);
    tlo = tiledlayout(fig, numFields, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    % Main Title logic
    mainTitle = sprintf('ROI %d: Gray Screen', iROI);
    if hasDark, mainTitle = [mainTitle ' (Solid) vs Darkness (Dashed)']; end
    title(tlo, mainTitle, 'FontWeight', 'bold');

    for f = 1:numFields
        fname = plotFields{f};
        ax = nexttile(tlo);
        hold(ax, 'on');
        
        % --- Always Plot Gray Screen ---
        renderSuperLayer(ax, respGray.tuningCurve, iROI, fname, 'k', '-');
        
        % --- Plot Darkness ONLY if present ---
        if hasDark
            renderSuperLayer(ax, respDark.tuningCurve, iROI, fname, [0.4 0.4 0.4], '--');
        end
        
        % Formatting
        ylabel(ax, fname, 'Interpreter', 'none');
        set(ax, 'XScale', 'log'); grid on;
        xticks(ax, [0.5, 1, 5, 10, 20, 40, 60, 100]);
        xticklabels(ax, {'0','1','5','10','20','40','60','100'});
        xlim(ax, [0.4, 110]);
        
        % --- Dynamic Significance Title ---
        sigG = respGray.tuningCurve.(fname).isSignificant_999(iROI);
        colG = '\color[rgb]{0, 0.6, 0}'; if ~sigG, colG = '\color[rgb]{0.8, 0, 0}'; end
        colB = '\color{black}';
        
        if hasDark
            sigD = respDark.tuningCurve.(fname).isSignificant_999(iROI);
            colD = '\color[rgb]{0, 0.6, 0}'; if ~sigD, colD = '\color[rgb]{0.8, 0, 0}'; end
            sigText = sprintf('%s%s Gray Sig%s,  %s%s Dark Sig%s', colB, colG, colB, colB, colD, colB);
        else
            sigText = sprintf('%s%s Gray Sig%s', colB, colG, colB);
        end
        title(ax, sigText, 'FontSize', 9, 'FontAngle', 'italic');
    end
    
    xlabel(tlo, 'Running Speed (cm/s)', 'FontWeight', 'bold');

    % --- Export ---
    if iROI == 1
        exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    else
        exportgraphics(fig, pdfPath, 'Append', true, 'ContentType', 'vector');
    end
    
    if mod(iROI, 50) == 0, fprintf('ROI %d/%d appended.\n', iROI, numROIs); end
end
close(fig);
end

function renderSuperLayer(ax, tc, roi, fn, col, lstyle)
    x_vals = [0.5, tc.speedBins(1:end-1) + diff(tc.speedBins)/2];
    y_mean = [tc.(fn).statMean(roi); tc.(fn).moveMean(roi, :)'];
    y_err  = [tc.(fn).statSEM(roi);  tc.(fn).moveSEM(roi, :)'];
    
    errorbar(ax, x_vals, y_mean, y_err, 'Color', col, 'LineStyle', lstyle, ...
        'LineWidth', 1.2, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', col);
    
    % Only add occupancy text for the top signal (dFF) to keep it clean
    if strcmp(fn, 'dFF')
        occ = [tc.occupancy.stationary, tc.occupancy.moving];
        for i = 1:length(x_vals)
            text(ax, x_vals(i), y_mean(i), sprintf('%.0fs', occ(i)), ...
                'VerticalAlignment', 'bottom', 'FontSize', 6, 'Color', col);
        end
    end
end