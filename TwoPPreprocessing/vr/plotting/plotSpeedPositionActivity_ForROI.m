function plotSpeedPositionActivity_ForROI(sessionFileInfo, response, targetROI, applySmoothing, smoothSigma)
% plotspeedpositionactivity_forroi - plots a 2d heatmap with a shared-axis 1d profile.
% ensures perfect vertical alignment by accounting for the colorbar width.

    %% set default options
    if nargin < 4, applySmoothing = true; end
    if nargin < 5, smoothSigma = [1.1, 1.5]; end
    
    %% 2. extract data
    tuningSurface = response.speedPositionActivity.matrix(:, :, targetROI);
    speedCenters  = response.speedPositionActivity.speedBinCenters;
    numPosBins    = size(tuningSurface, 2);

    %% smoothing logic 
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

    %% create figure
    figHandle = figure('Name', sprintf('Soma %d', targetROI), ...
        'Position', [100 100 650 550], 'Color', 'w');
    
    % define a shared width and left margin to ensure alignment
    leftMargin = 0.13;
    plotWidth  = 0.70; 

    % collapsed spatial profile ---
    ax1 = subplot('Position', [leftMargin, 0.72, plotWidth, 0.20]); 
    collapsedSpatial = mean(tuningSurface, 1, 'omitnan');
    plot(1:numPosBins, collapsedSpatial, 'k', 'LineWidth', 1.5);
    hold on;
    xline([40, 80, 120, 160, 200], '--', 'Color', [0.7 0.7 0.7], 'Alpha', 0.5);
    
    ylabel('\DeltaF/F');
    set(ax1, 'Box', 'off', 'TickDir', 'out', 'XTickLabel', [], 'XLim', [1, numPosBins]);
    title(sprintf('Soma %d: speed-position activity', targetROI));

    %speed-position heatmap
    ax2 = subplot('Position', [leftMargin, 0.12, plotWidth, 0.55]); 
    imagesc(1:numPosBins, speedCenters, tuningSurface);
    set(ax2, 'YDir', 'normal'); 
    
    colormap(parula);
    
    % contrast (98.5th percentile)
    activeData = tuningSurface(tuningSurface > 0);
    if ~isempty(activeData)
        maxVal = prctile(activeData, 98.5);
    else
        maxVal = 1;
    end
    set(ax2, 'CLim', [0, maxVal]);

    % add colorbar WITHOUT resizing the axis
    c = colorbar(ax2, 'Position', [leftMargin + plotWidth + 0.02, 0.12, 0.03, 0.55]);
    c.Label.String = '\DeltaF/F [NeuC]';

    %% reference lines and axes
    hold on;
    xline([40, 80, 120, 160, 200], '--w', 'Alpha', 0.4, 'LineWidth', 1.5);
    
    % x-axis formatting
    xticks([1 40 80 120 160 200]);
    xticklabels({'1', '40', '80', '120', '160', '200'});
    xlim([1, numPosBins]);

    % y-axis formatting (log scale)
    set(ax2, 'YScale', 'log'); 
    set(ax2, 'YMinorTick', 'on', 'TickDir', 'out');
    yticks([2, 5, 10, 20, 30]); 
    yticklabels({'2', '5', '10', '20', '30'});
    ylim([min(speedCenters), max(speedCenters)]);

    % labels and styling
    xlabel('Position (cm)');
    ylabel('Running speed (cm/s)');
    set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);

    %% link the x-axes for interactive use
    linkaxes([ax1, ax2], 'x');

    %% 7. directory setup and saving
    figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures', 'SpeedPositionMaps');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
    
    smoothLabel = 'raw'; if applySmoothing, smoothLabel = 'smoothed'; end
    saveName = sprintf('%s_%s_bouton%03d_%s_speedpos.png', ...
        sessionFileInfo.animal_name, sessionFileInfo.session_name, targetROI, smoothLabel);
    
    exportgraphics(figHandle, fullfile(figSaveDir, saveName), 'Resolution', 300);
    fprintf('figure saved for roi %d: %s\n', targetROI, saveName);
end