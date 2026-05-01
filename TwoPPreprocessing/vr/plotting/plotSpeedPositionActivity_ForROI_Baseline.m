function plotSpeedPositionActivity_ForROI_Baseline(sessionFileInfo, response, targetROI, applySmoothing, smoothSigma)
% plotSpeedPositionActivity_ForROI_Baseline - Heatmap for baseline trials only.
    %% set default options
    if nargin < 4, applySmoothing = true; end
    if nargin < 5, smoothSigma = [1.1, 1.5]; end
    
    %% 
    % Identify baseline indices
    conds = fieldnames(response.trialIndicesByCondition);
    baseIdx = find(contains(lower(conds), 'baseline'), 1);
    if isempty(baseIdx), error('No "baseline" condition found.'); end
    baseLaps = response.trialIndicesByCondition.(conds{baseIdx});
    
    % Reconstruct the tuning surface using ONLY baseline laps
    % Assuming response.speedPositionActivity stores raw lap-wise matrices:
    % (Lap x Position x ROI)
    rawLapData = response.speedPositionActivity.Matrix(baseLaps, :, :, targetROI);
    
    % Average across the baseline lap dimension to create the surface
    tuningSurface = squeeze(mean(rawLapData, 1, 'omitnan'));
    
    speedCenters  = response.speedPositionActivity.speedBinCenters;
    numPosBins    = size(tuningSurface, 2);

    %% 
    if applySmoothing
        mask = ~isnan(tuningSurface);
        dataZeroed = tuningSurface;
        dataZeroed(~mask) = 0; 
        blurredData = imgaussfilt(dataZeroed, smoothSigma, 'Padding', 'replicate');
        blurredMask = imgaussfilt(double(mask), smoothSigma, 'Padding', 'replicate');
        tuningSurface = blurredData ./ blurredMask;
        tuningSurface(isnan(tuningSurface)) = 0;
    else
        tuningSurface(isnan(tuningSurface)) = 0;
    end

    %%
    figHandle = figure('Name', sprintf('Baseline Soma %d', targetROI), ...
        'Position', [100 100 650 550], 'Color', 'w');
    
    leftMargin = 0.13;
    plotWidth  = 0.70; 
    
    % Collapsed spatial profile
    ax1 = subplot('Position', [leftMargin, 0.72, plotWidth, 0.20]); 
    collapsedSpatial = mean(tuningSurface, 1, 'omitnan');
    plot(1:numPosBins, collapsedSpatial, 'k', 'LineWidth', 1.5);
    hold on;
    xline([40, 80, 120, 160, 200], '--', 'Color', [0.7 0.7 0.7], 'Alpha', 0.5);
    
    ylabel('\DeltaF/F');
    set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
    title(sprintf('Bouton %d: Baseline Speed-Position', targetROI));

    % Speed-position heatmap
    ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.55]); 
    imagesc(1:numPosBins, speedCenters, tuningSurface);
    set(ax2, 'YDir', 'normal'); 
    colormap(parula);
    
    activeData = tuningSurface(tuningSurface > 0);
    maxVal = prctile(activeData, 98.5);
    if isempty(maxVal) || maxVal == 0, maxVal = 1; end
    set(ax2, 'CLim', [0, maxVal]);
    
    c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.03, 0.55]);
    c.Label.String = '\DeltaF/F';

    %% 
    hold on;
    xline([40, 80, 120, 160, 200], '--w', 'Alpha', 0.4, 'LineWidth', 1.5);
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    xlim([1, numPosBins]);
    
    set(ax2, 'YScale', 'log', 'YMinorTick', 'on', 'TickDir', 'out');
    yticks([2, 5, 10, 20, 30]); 
    ylim([min(speedCenters), max(speedCenters)]);
    
    xlabel('Position (cm)');
    ylabel('Running speed (cm/s)');
    set(ax2, 'Box', 'off', 'FontSize', 11);
    
    linkaxes([ax1, ax2], 'x');

    %% 5. Saving
    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionMaps', 'BaselineOnly');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
    
    saveName = sprintf('%s_%s_%03d_baseline_speedpos.png', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, targetROI);
    exportgraphics(figHandle, fullfile(figSaveDir, saveName), 'Resolution', 300);
end