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
    roi_matrix = normalize(roi_matrix, 2, 'range');
    mean_trace = mean(roi_matrix, 1, 'omitnan');
    std_trace  = std(roi_matrix, 0, 1, 'omitnan');
    sem_trace  = std_trace ./ sqrt(sum(~isnan(roi_matrix), 1));
    
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
    
    x_patch = [position_vector, fliplr(position_vector)];
    y_patch = [(mean_trace + sem_trace), fliplr(mean_trace - sem_trace)];
    
    nan_mask = isnan(x_patch) | isnan(y_patch);
    x_patch(nan_mask) = [];
    y_patch(nan_mask) = [];
    
    fill(x_patch, y_patch, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(position_vector, mean_trace, 'Color', 'k', 'LineWidth', 2);
    
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    
    xlabel('Position (cm)', 'FontSize', 12);
    ylabel('Mean \DeltaF/F', 'FontSize', 10);
    set(ax_trace, 'XTick', [40 80 120 160 200], 'XTickLabels', {'40', '80', '120', '160', '200'}, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);
    
    linkaxes([ax_heatmap, ax_trace], 'x');
    xlim(ax_trace, [1, num_positions]);
    ylim(ax_trace, [min(mean_trace - sem_trace) * 0.85, max(mean_trace + sem_trace) * 1.25]); %1.35
    
    pos_heat = get(ax_heatmap, 'Position');
    pos_trace = get(ax_trace, 'Position');
    set(ax_trace, 'Position', [pos_trace(1), pos_trace(2), pos_heat(3), pos_trace(4)]);
%     defaultAxesProperties(ax_trace);
    offsetAxes(ax_trace);
    
    hold off;
end



% %% raw version 
% function fig = plotRoiSpatialTuning(RSPData, sessionIdx, roiIdx, targetName)
%     lap_pos_activity = RSPData(sessionIdx).ConditionData.Baseline.LapActivity;
%     roi_matrix = squeeze(lap_pos_activity(roiIdx, :, :));
%     valid_laps = ~all(isnan(roi_matrix), 2);
%     roi_matrix = roi_matrix(valid_laps, :);
%     num_laps = size(roi_matrix, 1);
%     num_positions = size(roi_matrix, 2);
%     position_vector = 1:num_positions;
%     w = gausswin(15);
%     w = w / sum(w);
% 
%     for iLap = 1:num_laps
%         trace = roi_matrix(iLap, :);
%         if all(isnan(trace)), continue; end
%         nanMask = isnan(trace);
%         trace(nanMask) = 0;
%         smoothed = filtfilt(w, 1, trace);
%         smoothed(nanMask) = NaN;
%         roi_matrix(iLap, :) = smoothed;
%     end
% 
%     % no normalization applied - raw/smoothed activity used directly
% 
%     mean_trace = mean(roi_matrix, 1, 'omitnan');
%     std_trace  = std(roi_matrix, 0, 1, 'omitnan');
%     sem_trace  = std_trace ./ sqrt(sum(~isnan(roi_matrix), 1));
% 
%     fig = figure('Color', 'w', 'Position', [150, 150, 500, 580], 'Name', targetName);
%     ax_heatmap = axes('Position', [0.15, 0.42, 0.65, 0.50]);
%     imagesc(position_vector, 1:num_laps, roi_matrix);
%     colormap(ax_heatmap, flipud(gray));
% 
%     % Set color limits based on actual data range instead of fixed 0.25-0.75
%     clim_lo = prctile(roi_matrix(:), 1);
%     clim_hi = prctile(roi_matrix(:), 99);
%     clim([clim_lo, clim_hi]);
% 
%     cb = colorbar('Location', 'eastoutside');
%     cb.TickDirection = 'out';
%     cb.Box = 'off';
%     cb.FontName = 'Arial';
%     cb.FontSize = 10;
%     cb.Label.String = 'Activity (\DeltaF/F)';
%     cb.Label.FontName = 'Arial';
%     cb.Label.FontSize = 12;
%     cb.Label.Rotation = 90;
%     cb.Label.VerticalAlignment = 'bottom';
% 
%     hold on;
%     xline(40, 'k--', 'LineWidth', 1.5);
%     xline(80, 'k--', 'LineWidth', 1.5);
%     xline(120, 'k--', 'LineWidth', 1.5);
%     xline(160, 'k--', 'LineWidth', 1.5);
%     hold off;
% 
%     ylabel('Lap #', 'FontSize', 12);
%     title(sprintf('R^2 = %.4f  |  SMI = %.4f', RSPData(sessionIdx).cvExpVar.meanExpVar(roiIdx), RSPData(sessionIdx).SMI.SMI(roiIdx)), 'FontSize', 12, 'FontWeight', 'normal');
% 
%     set(ax_heatmap, 'XTick', [40 80 120 160 200], 'XTickLabel', [], ...
%         'YTick', [1, num_laps], 'YTickLabel', { '1', num2str(num_laps) }, ...
%         'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
% 
%     ax_trace = axes('Position', [0.15, 0.10, 0.65, 0.25]);
%     hold on;
%     x_patch = [position_vector, fliplr(position_vector)];
%     y_patch = [(mean_trace + sem_trace), fliplr(mean_trace - sem_trace)];
%     nan_mask = isnan(x_patch) | isnan(y_patch);
%     x_patch(nan_mask) = [];
%     y_patch(nan_mask) = [];
%     fill(x_patch, y_patch, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     plot(position_vector, mean_trace, 'Color', 'k', 'LineWidth', 2);
%     xline(40, 'k--', 'LineWidth', 1.5);
%     xline(80, 'k--', 'LineWidth', 1.5);
%     xline(120, 'k--', 'LineWidth', 1.5);
%     xline(160, 'k--', 'LineWidth', 1.5);
%     xlabel('Position (cm)', 'FontSize', 12);
%     ylabel('Mean \DeltaF/F', 'FontSize', 10);
%     set(ax_trace, 'XTick', [40 80 120 160 200], 'XTickLabels', {'40', '80', '120', '160', '200'}, 'TickDir', 'out', 'box', 'off', 'FontSize', 12);
%     linkaxes([ax_heatmap, ax_trace], 'x');
%     xlim(ax_trace, [1, num_positions]);
%     ylim(ax_trace, [min(mean_trace - sem_trace) * 0.85, max(mean_trace + sem_trace) * 1.25]);
% 
%     pos_heat = get(ax_heatmap, 'Position');
%     pos_trace = get(ax_trace, 'Position');
%     set(ax_trace, 'Position', [pos_trace(1), pos_trace(2), pos_heat(3), pos_trace(4)]);
%     offsetAxes(ax_trace);
%     hold off;
% end
