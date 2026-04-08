function plotPeakSignificance(response, ROI_idx, signalName)
    % Extracts data for a specific ROI to visualize the significance test
    realPeak = max(mean(response.lapPositionActivity.(signalName)(ROI_idx, :, :), 2, 'omitnan'));
    shuffDist = squeeze(max(response.lapPositionActivity_ShuffleMatrix.(signalName)(ROI_idx, :, :), [], 2));
    threshold = prctile(shuffDist, 95);
    
    figure('Position', [100, 100, 900, 400], 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact');

    % --- Left: The Null Distribution ---
    nexttile;
    histogram(shuffDist, 30, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none'); hold on;
    xline(threshold, '--r', 'LineWidth', 2, 'Label', '95th Percentile');
    plot(realPeak, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    text(realPeak, 2, ' Real Peak', 'Color', 'g', 'FontWeight', 'bold');
    
    title(['ROI ' num2str(ROI_idx) ': Peak Distribution']);
    xlabel('Peak Activity Magnitude'); ylabel('Shuffle Count');
    grid on;

    % --- Right: Tuning Curve vs. Shuffles ---
    nexttile;
    % Plot first 10 shuffles for visual reference
    shuffCurves = squeeze(response.lapPositionActivity_ShuffleMatrix.(signalName)(ROI_idx, :, 1:10));
    plot(shuffCurves, 'Color', [0.8 0.8 0.8, 0.5]); hold on;
    
    % Plot real tuning curve
    realCurve = squeeze(mean(response.lapPositionActivity.(signalName)(ROI_idx, :, :), 2, 'omitnan'));
    plot(realCurve, 'g', 'LineWidth', 2.5);
    
    title('Real Map vs. Shuffled Maps');
    xlabel('Position Bins'); ylabel('Activity');
    legend('Shuffled', '', '', '', '', '', '', '', '', '', 'Real Tuning', 'Location', 'best');
end