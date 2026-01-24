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
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    
    % Ticks
    set(ax, 'XTick', [0 50 70 90 110 140], ...
            'XTickLabel', {'0','50','70','90','110','140'});
    xlabel('Position (cm)');
    colorbar;
end