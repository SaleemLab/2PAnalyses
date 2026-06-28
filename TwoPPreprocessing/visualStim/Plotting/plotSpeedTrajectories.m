function figHandle = plotSpeedTrajectories(response)
    if ~isfield(response, 'lapPositionRunningSpeed')
        error('Input structure must contain the field "lapPositionRunningSpeed".');
    end

    binnedSpeed = response.lapPositionRunningSpeed;
    [nLaps, numBins] = size(binnedSpeed);

    maxSpeedInData = max(binnedSpeed(:), [], 'omitnan');
    if isempty(maxSpeedInData) || isnan(maxSpeedInData) || maxSpeedInData < 30
        maxSpeedInData = 30;
    end
    
    paddingCeiling = maxSpeedInData + 2; 

    figHandle = figure('Color', 'w', 'Position', [200, 200, 500, 380]);
    hold on;

    x_box = [1, numBins, numBins, 1];
    patch(x_box, [1, 1, 15, 15], [0.68, 0.92, 0.98], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    patch(x_box, [15, 15, 23, 23], [0.88, 0.88, 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    patch(x_box, [23, 23, paddingCeiling, paddingCeiling], [0.96, 0.64, 0.76], 'EdgeColor', 'none', 'FaceAlpha', 0.45);

    x_positions = 1:numBins;
    for iLap = 1:nLaps
        lapTrace = binnedSpeed(iLap, :);
        if all(isnan(lapTrace)), continue; end
        plot(x_positions, lapTrace, 'Color', [0.1, 0.1, 0.1, 0.3], 'LineWidth', 1);
    end

    meanSpeedProfile = mean(binnedSpeed, 1, 'omitnan');
    plot(x_positions, meanSpeedProfile, 'Color', 'k', 'LineWidth', 3);

    xlabel('Position (cm)', 'FontName', 'Arial', 'FontSize', 12);
    ylabel('Speed (cm/s)', 'FontName', 'Arial', 'FontSize', 12);

    xlim([1, numBins]);
    if numBins > 150
        xticks([1 40 80 120 160 200]);
        xticklabels({'1', '40', '80', '120', '160', '200'});
    else
        xticks([1 40 80 120 140]);
        xticklabels({'1', '40', '80', '120', '140'});
    end

    ylim([1, paddingCeiling]);
    yticks([1, 15, 30, 40]);

    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca, true); end
    if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end

    box off;
    hold off;
end