function plotSingleROIReceptiveField(response, roiIndices)
% Normalizes each trace by its specific pre-stimulus baseline
% to visualize activity increases relative to stimulus onset.

psth = response.psthData;
stimVs = vertcat(psth.stimValue);

% Extract unique centers from your data
az_centers = sort(unique(stimVs(:,1)), 'ascend'); 
el_centers = sort(unique(stimVs(:,2)), 'ascend'); 
nAz = numel(az_centers);
nEl = numel(el_centers);

% Time and Window definitions
tVec = psth(1).timeVector;
tMask = tVec >= -1 & tVec <= 3.0; 
baseMask = tVec >= -1 & tVec <= 0; % Pre-stimulus baseline window

figure('Color','w','Position',[100 100 600 850]);
tlo = tiledlayout(numel(roiIndices), 1, 'Padding', 'compact', 'TileSpacing', 'none');

for i = 1:numel(roiIndices)
    roiIdx = roiIndices(i);
    nexttile(tlo); hold on;
    
    % --- 1. Calculate the Heatmap (0.5 to 2.0s mean response) ---
    rfMap = zeros(nEl, nAz);
    for k = 1:numel(psth)
        r = find(el_centers == psth(k).stimValue(2));
        c = find(az_centers == psth(k).stimValue(1));
        
        mu = median(squeeze(psth(k).alignedResponses(roiIdx,:,:)), 2, 'omitnan');
        % Subtract pre-stimulus baseline
        F0 = mean(mu(baseMask), 'omitnan');
        zeroed = mu - F0;
        
        rfMap(r, c) = mean(zeroed(tVec >= 0.5 & tVec <= 2.0), 'omitnan');
    end
    
    % --- 2. Plot Grayscale Background ---
    imagesc(az_centers, el_centers, rfMap);
    colormap(gray);
    clim([0 0.5]); % Scaling to match your reference image
    
    % --- 3. Plot Baseline-Subtracted Traces ---
    % Scale so a 0.5 dF/F signal is clearly visible within the row height
    % (Avg distance between your elevations is ~22 deg)
    traceScale = 25; 
    
    for k = 1:numel(psth)
        pos = psth(k).stimValue;
        mu = median(squeeze(psth(k).alignedResponses(roiIdx,:,:)), 2, 'omitnan');
        
        % Normalize trace by subtracting pre-stimulus mean
        F0 = mean(mu(baseMask), 'omitnan');
        normTrace = mu(tMask) - F0;
        
        % Time-to-Azimuth mapping (Starts at Az center, extends right)
        plotX = pos(1) + (tVec(tMask) * 7); 
        
        % Amplitude-to-Elevation mapping (Baseline is exactly at pos(2))
        plotY = pos(2) + (normTrace * traceScale);
        
        % Plotting the normalized trace
        plot(plotX, plotY, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.3);
    end
    
    % --- Formatting ---
    set(gca, 'XDir', 'reverse', 'YDir', 'normal');
    set(gca, 'XTick', az_centers, 'YTick', el_centers);
    ylabel('elevation (°)');
    title(sprintf('Bouton %d', roiIdx), 'Units', 'normalized', 'Position', [0.01 0.9], 'FontWeight', 'bold');
    
    cb = colorbar; 
    ylabel(cb, 'Baseline-Subtracted \DeltaF/F');
    
    if i < numel(roiIndices), set(gca, 'XTickLabel', []); end
end
xlabel('azimuth (°)');
end