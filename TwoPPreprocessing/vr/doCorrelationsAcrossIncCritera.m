function corrs = doCorrelationsAcrossIncCritera(sessionFileInfo, signalToUse, doPlot)
% DOCORRELATIONSACROSSINCCRITERA Calculates and plots a standard cross-correlation 
% matrix across ROI metrics.

%% 
if nargin < 2 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 3 || isempty(doPlot);  doPlot = 1; end

ROIDataFile = sessionFileInfo.otherSessFilePaths.sessionROIData;
corrs = struct();

if exist(ROIDataFile, 'file') ~= 2
    warning('Target ROI data file does not exist: %s', ROIDataFile);
    return;
end

%%
loadedMetrics = load(ROIDataFile);

medianEV         = loadedMetrics.crossValExpVar.medianExpVar; 
varToTuningVar   = loadedMetrics.tuningCurveVariance.ratioVarToTuningVar;
varToTuningRange = loadedMetrics.tuningCurveVariance.ratioVarToTuningRange; 
oddEvenRho       = loadedMetrics.lapCorr_OddEven.rho;
halvesRho        = loadedMetrics.lapCorr_Halves.rho;

roisToKeepMask = false(size(varToTuningVar));
if isfield(loadedMetrics.highlyCorrBoutons, 'roisToKeep')
    roisToKeepMask(loadedMetrics.highlyCorrBoutons.roisToKeep) = true;
end

%% 
% exclude the shuffles and unique boutons critera 
metricNames = {'Median_EV', 'Var_to_TuningVar', 'Var_to_TuningRange', ...
               'OddEven_Rho', 'Halves_Rho'};

metricTable = table(medianEV(:), varToTuningVar(:), varToTuningRange(:), ...
                    oddEvenRho(:), halvesRho(:), ...
                    'VariableNames', metricNames);

%% cross correlation (spearman pairwise) 
[corrMatrix, pMatrix] = corr(table2array(metricTable), 'Type', 'Spearman', 'Rows', 'pairwise');

corrs.metricTable = metricTable;
corrs.corrMatrix = corrMatrix;
corrs.pMatrix = pMatrix;
corrs.metricNames = metricNames;

%% heatmap
if doPlot
    figure('Color', 'w', 'Position', [150, 150, 750, 650]);
    
    imagesc(corrMatrix);
    colorbar;
    clim([-1 1]);
    
    % High-contrast custom diverging map (Blue [-1] -> White [0] -> Red
    % [+1]) ; gemini 
    try
        colormap(gca, bluered(256));
    catch
        r = [zeros(1,128), linspace(0,1,128)];
        g = [linspace(0,1,128), linspace(1,0,128)];
        b = [linspace(1,0,128), zeros(1,128)];
        colormap(gca, [r; g; b]'/1.1); 
    end
    

    numMetrics = numel(metricNames);
    set(gca, 'XTick', 1:numMetrics, 'XTickLabel', strrep(metricNames, '_', ' '), 'XTickLabelRotation', 45, ...
             'YTick', 1:numMetrics, 'YTickLabel', strrep(metricNames, '_', ' '), ...
             'TickLabelInterpreter', 'none', 'FontSize', 10);
         
    title('Cross-correlation matrix (Spearman rho)', 'FontSize', 12);
    axis square;
    set(gca, 'GridColor', 'w', 'GridAlpha', 0.4, 'Layer', 'top');
    
    % 
    for i = 1:numMetrics
        for j = 1:numMetrics
            val = corrMatrix(i,j);
            text(j, i, sprintf('%.2f', val), ...
                'HorizontalAlignment', 'center', ...
                'Color', inlineIf(abs(val) > 0.55, 'w', 'k'), ...
                'FontSize', 10, ...
                'FontWeight', inlineIf(abs(val) > 0.7, 'bold', 'normal'));
        end
    end
end
end

%%helper 
function val = inlineIf(condition, trueVal, falseVal)
    if condition; val = trueVal; else; val = falseVal; end
end