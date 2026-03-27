function plotTuningInAxes(ax, oddData, evenData)
% Helper to standardize heatmap plotting across all functions

    % 1. Force this axis to be the current one. 
    % This fixes almost all "Incorrect handle" errors for subsequent commands.
    axes(ax); 

    normOdd = normalize(oddData, 2, 'range');
    normEven = normalize(evenData, 2, 'range');
    
    % Sort based on Odd peaks
    [~, peakIdx] = max(normOdd, [], 2, 'omitnan');
    [~, sortIdx] = sort(peakIdx);
    
    % Plot Even sorted by Odd into the current axes (ax)
    imagesc(normEven(sortIdx, :));
    
    % Standard Formatting using universally compatible commands
    colormap(flipud(gray)); 
    caxis([0 1]); % Implicit call is often safer than explicit caxis(ax, [0 1])

    % For properties, explicit handles still work best
    set(ax, 'TickDir', 'out', 'box', 'off', 'FontSize', 10, 'YDir', 'normal');
    
    % Landmarks
    xline(40, 'k--', 'LineWidth', 1.5);
    xline(80, 'k--', 'LineWidth', 1.5);
    xline(120, 'k--', 'LineWidth', 1.5);
    xline(160, 'k--', 'LineWidth', 1.5);
    
    % Ticks
    set(ax, 'XTick', [1 40 80 120 160 200], ...
            'XTickLabel', {'1','40','80','120','160','200'});
    xlabel('Position (cm)');
    colorbar;
end