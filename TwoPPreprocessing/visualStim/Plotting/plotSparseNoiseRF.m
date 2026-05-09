function plotSparseNoiseRF(neuronIdx, sparseNoiseRF)
    % plotSparseNoiseRF
    
    % Check if data exists
    if ~isempty(sparseNoiseRF.initMap) && neuronIdx <= length(sparseNoiseRF.initMap) && ...
       ~isempty(sparseNoiseRF.initMap{neuronIdx})
        
        % Extract RF data
        rfRaw = flipud(sparseNoiseRF.initMap{neuronIdx}(:, :, end, 4));
        
        % Plot on current axes
        imagesc(linspace(-70, 20, size(rfRaw, 2)), ...
                linspace(-20, 40, size(rfRaw, 1)), ...
                imgaussfilt(rfRaw, 1));
        
        set(gca, 'YDir', 'normal'); 
        colormap(gca, 'parula'); 
        colorbar;
        
        hold on;
        xline(0, 'k:', 'Alpha', 0.5); 
        yline(0, 'k:', 'Alpha', 0.5);
        hold off;
        
        xlabel('Azimuth (\circ)'); 
        ylabel('Elevation (\circ)');
        title(['Sparse Noise RF (SVD) - ROI ' num2str(neuronIdx)]);
    else
        % If no data
        cla; % Clear current axis
        text(0.5, 0.5, 'No Sparse Noise Data', 'HorizontalAlignment', 'center');
        axis off;
    end
end