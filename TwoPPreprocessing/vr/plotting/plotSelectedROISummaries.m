function plotSelectedROISummaries(sessionFileInfo, response, selectedROIs, applySmoothing, signalToUse)
% 

if nargin < 4, applySmoothing = true; end
if nargin < 5, signalToUse = 'dFF'; end

%%
outputFilePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
metricsLoaded = false;
if exist(outputFilePath, 'file') == 2
    loadedMetrics = load(outputFilePath); 
    peakSigIndex = loadedMetrics.nullDist_PeakTuningMetric.isSignificantByPeakShuffling.(signalToUse);
    rangeSigIndex = loadedMetrics.nullDist_RangeTuningMetric.isSignificantByRange.(signalToUse);
    varToTuningVar = loadedMetrics.tuningCurveVariance.ratioVarToTuningVar;
    varToTuningRange = loadedMetrics.tuningCurveVariance.ratioVarToTuningRange;
    halvesCorrIdxRho = loadedMetrics.lapCorr_Halves.rho;
    oddEvenCorrIdxRho = loadedMetrics.lapCorr_OddEven.rho; 
    
    roisToKeepMask = false(size(varToTuningVar));
    if isfield(loadedMetrics.highlyCorrBoutons, 'roisToKeep')
        roisToKeepMask(loadedMetrics.highlyCorrBoutons.roisToKeep) = true;
    end
    metricsLoaded = true;
end

%% 
saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\DistibutionsAllCritera';
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

savePath = fullfile(saveDir, sprintf('%s_%s_selectedROIs_summary.png', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name));

lapActivityFull = response.lapPositionActivity.(signalToUse);
numROIs = length(selectedROIs);

rowHeight = 400; % Increased row height for better resolution
totalHeight = rowHeight * numROIs + 150; 
figWidth = 1200; 

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 figWidth totalHeight]);
t = tiledlayout(numROIs, 2, 'TileSpacing', 'loose', 'Padding', 'compact');


sgtitle(t, sprintf('Selected ROI Summary: %s | %s', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name), ...
    'FontWeight', 'bold', 'FontSize', 22);


T_VV = 20; T_VR = 1; T_HALVES = 0.4;

%% 
for i = 1:numROIs
    neuronIdx = selectedROIs(i);
    roiActivity = squeeze(lapActivityFull(neuronIdx, :, :));
    
    if applySmoothing
        w = gausswin(5); w = w / sum(w);
        for iL = 1:size(roiActivity, 1)
            trace = roiActivity(iL, :);
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            trace = filtfilt(w, 1, trace); trace(nanMask) = NaN;
            roiActivity(iL, :) = trace;
        end
    end
    
    meanActivity = mean(roiActivity, 1, 'omitnan');
    semActivity = std(roiActivity, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiActivity), 1));
    normLapActivity = normalize(roiActivity, 2, 'range');

    % 
    nexttile; hold on;
    x = 1:size(meanActivity, 2);
    fill([x fliplr(x)], [meanActivity + semActivity, fliplr(meanActivity - semActivity)], ...
         [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, meanActivity, 'k', 'LineWidth', 2.5);
    xline([50 70 90 110], 'k--', 'LineWidth', 1.5); 
    
    % 
    xticks([0 50 70 90 110 140]);
    set(gca, 'FontSize', 14); 
    ylabel('Mean dFF', 'FontSize', 16);
    xlabel('Position (cm)', 'FontSize', 16);
    
    if metricsLoaded
        passes = (peakSigIndex(neuronIdx) || rangeSigIndex(neuronIdx)) && ...
                 (varToTuningVar(neuronIdx) <= T_VV) && (varToTuningRange(neuronIdx) <= T_VR) && ...
                 roisToKeepMask(neuronIdx) && (halvesCorrIdxRho(neuronIdx) >= T_HALVES);
        tCol = [0 0.5 0]; if ~passes, tCol = [0.8 0 0]; end
        
        critTxt = sprintf('SigP|R: %d|%d, V/V: %.1f, V/R: %.1f, Half: %.2f, OdEv: %.2f, Matched:%d', ...
            peakSigIndex(neuronIdx), rangeSigIndex(neuronIdx), varToTuningVar(neuronIdx), ...
            varToTuningRange(neuronIdx), halvesCorrIdxRho(neuronIdx), oddEvenCorrIdxRho(neuronIdx), ...
            roisToKeepMask(neuronIdx));
            
        % ROI Title and Subtitle (Increased Sizes)
        title(sprintf('ROI %d', neuronIdx), 'Color', tCol, 'FontSize', 18);
        subtitle(critTxt, 'FontSize', 12, 'FontAngle', 'italic', 'FontWeight', 'bold');
    end

    %
    nexttile;
    imagesc(normLapActivity);
    caxis([0 1]); colormap(flipud(gray));
    xline([50 70 90 110], 'k--', 'LineWidth', 1.5);
    
    xticks([0 50 70 90 110 140]);
    set(gca, 'FontSize', 14);
    title('Lap-by-position activity', 'FontSize', 16);
    ylabel('Lap #', 'FontSize', 16);
    xlabel('Position (cm)', 'FontSize', 16);
end

%%
exportgraphics(fig, savePath, 'Resolution', 300);
close(fig);
fprintf('Figure saved with increased text size to: %s\n', savePath);
end
