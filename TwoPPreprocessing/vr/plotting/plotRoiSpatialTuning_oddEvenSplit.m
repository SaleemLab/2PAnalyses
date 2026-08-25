function fig = plotRoiSpatialTuning(RSPData, sessionIdx, roiIdx, targetName)
    lap_pos_activity = RSPData(sessionIdx).ConditionData.Baseline.LapActivity;
    roi_matrix = squeeze(lap_pos_activity(roiIdx, :, :));
    valid_laps = ~all(isnan(roi_matrix), 2);
    roi_matrix = roi_matrix(valid_laps, :);
    num_laps = size(roi_matrix, 1);
    num_positions = size(roi_matrix, 2);
    position_vector = 1:num_positions;
    w = gausswin(15);
    w = w / sum(w);
    for iLap = 1:num_laps
        trace = roi_matrix(iLap, :);
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace);
        trace(nanMask) = 0;
        smoothed = filtfilt(w, 1, trace);
        smoothed(nanMask) = NaN;
        roi_matrix(iLap, :) = smoothed;
    end

    % --- Odd/even lap split (BEFORE normalization) ---
    odd_laps  = 1:2:num_laps;
    even_laps = 2:2:num_laps;

    % Derive normalization range from ODD laps only, apply to the whole matrix
    odd_min = min(roi_matrix(odd_laps, :), [], 'all', 'omitnan');
    odd_max = max(roi_matrix(odd_laps, :), [], 'all', 'omitnan');
    roi_matrix = (roi_matrix - odd_min) ./ (odd_max - odd_min);

    roi_matrix_odd  = roi_matrix(odd_laps, :);
    roi_matrix_even = roi_matrix(even_laps, :);

    mean_trace_odd  = mean(roi_matrix_odd, 1, 'omitnan');
    sem_trace_odd   = std(roi_matrix_odd, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roi_matrix_odd), 1));

    mean_trace_even = mean(roi_matrix_even, 1, 'omitnan');
    sem_trace_even  = std(roi_matrix_even, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roi_matrix_even), 1));

    fig = figure('Color', 'w', 'Position', [150, 150, 500, 580], 'Name', targetName);
    ax_heatmap = axes('Position', [0.15, 0.42, 0.65, 0.50]);
    imagesc(position_vector, 1:num_laps, roi_matrix);
    colormap(ax_heatmap, flipud(gray));
    clim([0.25 0.75]);
    cb = colorbar('Location', 'eastoutside');
    cb.Ticks = [0.25 0.50 0.75];
    cb.TickLabels = {'0.25', '0.50', '0.75'};
    cb.TickDirection = 'out';
    cb.Box = 'off';
    cb.FontName = 'Arial';
    cb.FontSize = 10;
    cb.Label.String = 'Activity (normalised)';
    cb.Label.FontName = 'Arial';
    cb.Label.FontSize = 12;
    cb.Label.Rotation = 90;
    cb.Label.VerticalAlignment = 'bottom';
    hold on;
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    hold off;
    ylabel('Lap #', 'FontSize', 12);
    title(sprintf('R^2 = %.4f  |  SMI = %.4f', RSPData(sessionIdx).cvExpVar.meanExpVar(roiIdx), RSPData(sessionIdx).SMI.SMI(roiIdx)), 'FontSize', 12, 'FontWeight', 'normal');
    % Enforce only the first and last lap numbers on the Y-axis
    set(ax_heatmap, 'XTick', [40 80 120 160 200], 'XTickLabel', [], ...
        'YTick', [1, num_laps], 'YTickLabel', { '1', num2str(num_laps) }, ...
        'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');

    ax_trace = axes('Position', [0.15, 0.10, 0.65, 0.25]);
    hold on;

    % --- Odd laps shading + mean ---
    x_patch_odd = [position_vector, fliplr(position_vector)];
    y_patch_odd = [(mean_trace_odd + sem_trace_odd), fliplr(mean_trace_odd - sem_trace_odd)];
    nan_mask_odd = isnan(x_patch_odd) | isnan(y_patch_odd);
    x_patch_odd(nan_mask_odd) = [];
    y_patch_odd(nan_mask_odd) = [];
    fill(x_patch_odd, y_patch_odd, [0.85 0.33 0.10], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    h_odd = plot(position_vector, mean_trace_odd, 'Color', [0.85 0.33 0.10], 'LineWidth', 2);

    % --- Even laps shading + mean ---
    x_patch_even = [position_vector, fliplr(position_vector)];
    y_patch_even = [(mean_trace_even + sem_trace_even), fliplr(mean_trace_even - sem_trace_even)];
    nan_mask_even = isnan(x_patch_even) | isnan(y_patch_even);
    x_patch_even(nan_mask_even) = [];
    y_patch_even(nan_mask_even) = [];
    fill(x_patch_even, y_patch_even, [0 0.45 0.74], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    h_even = plot(position_vector, mean_trace_even, 'Color', [0 0.45 0.74], 'LineWidth', 2);

    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    xlabel('Position (cm)', 'FontSize', 12);
    ylabel('Mean \DeltaF/F', 'FontSize', 10);
    legend([h_odd, h_even], {'Odd laps', 'Even laps'}, 'Location', 'best', 'Box', 'off', 'FontSize', 9);
    set(ax_trace, 'XTick', [40 80 120 160 200], 'XTickLabels', {'40', '80', '120', '160', '200'}, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);
    linkaxes([ax_heatmap, ax_trace], 'x');
    xlim(ax_trace, [1, num_positions]);

    all_lo = min([mean_trace_odd - sem_trace_odd, mean_trace_even - sem_trace_even]);
    all_hi = max([mean_trace_odd + sem_trace_odd, mean_trace_even + sem_trace_even]);
    ylim(ax_trace, [all_lo * 0.85, all_hi * 1.25]);

    pos_heat = get(ax_heatmap, 'Position');
    pos_trace = get(ax_trace, 'Position');
    set(ax_trace, 'Position', [pos_trace(1), pos_trace(2), pos_heat(3), pos_trace(4)]);
    offsetAxes(ax_trace);
    hold off;
end