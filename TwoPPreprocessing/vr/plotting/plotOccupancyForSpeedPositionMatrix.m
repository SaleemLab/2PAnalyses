function plotOccupancyForSpeedPositionMatrix(sessionfileinfo, speedData, speedEdges, speedCenters)
% plotoccupancyvalidation - creates a diagnostic figure to prove speed binning logic.
% helps verify if every speed row has equal data (quantile) or equal width (fixed).

    % 1. get dimensions from the data
    [numLaps, numPosBins] = size(speedData);
    nSpeedBins = length(speedCenters);

    % 2. calculate occupancy matrix
    % logic: count how many frames fall into every speed x position box
    debugOccupancy = zeros(nSpeedBins, numPosBins);
    for b = 1:numPosBins
        binSpeeds = speedData(:, b);
        for s = 1:nSpeedBins
            % find frames where speed falls within the current bin edges
            idx = binSpeeds >= speedEdges(s) & binSpeeds < speedEdges(s+1);
            debugOccupancy(s, b) = sum(idx); 
        end
    end

    % create the visual proof figure
    figHandle = figure('Name', 'Occupancy Normalization Check', 'Position', [150 150 900 450], 'Color', 'w');

    % left subplot: the 2d occupancy map
    % shows exactly where the mouse spent its time on the track
    subplot(1, 4, 1:3);
    imagesc(1:numPosBins, speedCenters, debugOccupancy);
    set(gca, 'YScale', 'log', 'YDir', 'normal');
    colormap(hot);
    c = colorbar('Location', 'southoutside');
    c.Label.String = 'Samples per Bin (Frames)';
    xlabel('Position (cm)'); 
    ylabel('Speed (cm/s)');
    title(sprintf('Animal: %s | Speed-Position Occupancy', sessionfileinfo.animal_name), 'Interpreter', 'none');

    % right subplot: the "proof" bar chart
    % if using quantilebins, these bars should all be exactly the same length
    subplot(1, 4, 4);
    totalSamplesPerSpeedBin = sum(debugOccupancy, 2);
    barh(speedCenters, totalSamplesPerSpeedBin, 'FaceColor', [0.8 0.2 0.2]);
    set(gca, 'YScale', 'log', 'YDir', 'normal');
    xlabel('Total Samples');
    title('Total Data per Speed Row');
    grid on;

    % 4. save diagnostic figure
    figSaveDir = fullfile(sessionfileinfo.Directories.save_folder, 'Figures', 'Diagnostics');
    if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
    
    saveName = sprintf('%s_%s_OccupancyValidation.png', ...
        sessionfileinfo.animal_name, sessionfileinfo.session_name);
    exportgraphics(figHandle, fullfile(figSaveDir, saveName), 'Resolution', 300);
    
    fprintf('occupancy validation saved to: %s\n', saveName);
end