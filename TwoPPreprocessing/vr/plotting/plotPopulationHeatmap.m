function plotPopulationHeatmap(tuningCurves, positionAxis, animalName, sessionName)
    % plotPopulationHeatmap Generates a population response heatmap with a specific style.
    %
    % INPUTS:
    %   tuningCurves - 2D matrix (neurons x position bins) of final tuning curves.
    %   positionAxis - 1D vector for the x-axis (position bins).
    %   animalName   - (Optional) String for the figure title.
    %   sessionName  - (Optional) String for the figure title.

    if nargin < 3, animalName = 'Animal'; end
    if nargin < 4, sessionName = 'Session'; end
    
    %% --- 1. Normalize and Sort the Data ---
    numNeurons = size(tuningCurves, 1);
    normAll = normalize(tuningCurves, 2, 'range');
    [~, peakIdx] = max(normAll, [], 2);
    [~, sortIdx] = sort(peakIdx);

    %% --- 2. Create the Plot with the Specified Style ---
    figure();
    imagesc(positionAxis, 1:numNeurons, normAll(sortIdx, :));
    
    caxis([0 1]);
    colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 18, 'YDir', 'normal');
    
    xline(50, 'w--', 'LineWidth', 2.0);
    xline(70, 'w--', 'LineWidth', 2.0);
    xline(90, 'w--', 'LineWidth', 2.0);
    xline(110, 'w--', 'LineWidth', 2.0);
    
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    xlim([0 140]);
    
    xlabel('Position (cm)');
    ylabel('Sorted ROIs');
    title([animalName ' - ' sessionName ' - All laps sorted']);
    
    h = colorbar;
    ylabel(h, 'Activity (normalised)');
end