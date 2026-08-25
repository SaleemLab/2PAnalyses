% function plotAllNeuronSummariesToPDF_SparseNoiseAndVR(sessionFileInfo, response, applySmoothing)
% % UPDATED: Now includes Receptive Field subplot if sparseNoiseRF exists
% 
% if nargin < 3, applySmoothing = true; end
% 
% %% Output path
% figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
% if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
% pdfPath = fullfile(figSaveDir, sprintf('%s_%s_AllROIs_Summary_withRF_dFF.pdf', ...
%     sessionFileInfo.animal_name, sessionFileInfo.session_name));
% 
% %% Load RF data if available
% rfDataAvailable = false;
% if isfield(sessionFileInfo.otherSessFilePaths, 'sessionROIData') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file')
%     vars = who('-file', sessionFileInfo.otherSessFilePaths.sessionROIData);
%     if ismember('sparseNoiseRF', vars)
%         load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'sparseNoiseRF');
%         rfDataAvailable = true;
%     end
% end
% 
% lapActivityFull = response.lapPositionActivity.dFF;
% nROIs = size(lapActivityFull, 1);
% 
% for neuronIdx = 1:nROIs
%     roiActivity = squeeze(lapActivityFull(neuronIdx, :, :));
%     if all(isnan(roiActivity), 'all'), continue; end
% 
%     % [Smoothing Logic - remains same as your original]
%     if applySmoothing
%         w = gausswin(10); w = w / sum(w);
%         for iLap = 1:size(roiActivity, 1)
%             trace = roiActivity(iLap, :);
%             if all(isnan(trace)), continue; end
%             nanMask = isnan(trace); trace(nanMask) = 0;
%             smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
%             roiActivity(iLap, :) = smoothed;
%         end
%     end
% 
%     meanActivity = mean(roiActivity, 1, 'omitnan');
%     semActivity = std(roiActivity, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiActivity), 1));
%     normLapActivity = normalize(roiActivity, 2, 'range');
% 
%     % Create figure - wider to accommodate 3 plots
%     fig = figure('Visible', 'off', 'Position', [50 100 1600 450]);
% 
%     %  Mean ± SEM Trace
%     subplot(1, 3, 1); hold on;
%     x = 1:size(meanActivity, 2);
%     fill([x fliplr(x)], [meanActivity + semActivity, fliplr(meanActivity - semActivity)], [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     plot(x, meanActivity, 'k', 'LineWidth', 2);
%     xline(40, 'k--'); xline(80, 'k--'); xline(120, 'k--'); xline(160, 'k--');
%     xticks([1 40 80 120 160 200]); xlabel('Position (cm)'); ylabel('dFF');
%     title(sprintf('ROI %d', neuronIdx));
% 
%     % Heatmap ---
%     subplot(1, 3, 2);
%     imagesc(normLapActivity); caxis([0 1]); colormap(subplot(1,3,2), flipud(gray));
%     xline(40, 'k--'); xline(80, 'k--'); xline(120, 'k--'); xline(160, 'k--');
%     xticks([1 40 80 120 160 200]); xlabel('Position (cm)'); ylabel('Lap #');
%     title('Lap-Position-Activity');
%     colorbar; ylabel(colorbar, sprintf('dFF'));
% 
%     % snt
%     subplot(1, 3, 3);
%     if rfDataAvailable && ~isempty(sparseNoiseRF.initMap{neuronIdx})
%         % Extract RF (SVD map from contrast channel)
%         % Map structure: [elevation x azimuth x time/SVD x channel]
%         % Contrast is the 4th channel
%         rfRaw = sparseNoiseRF.initMap{neuronIdx}(:, :, end, 4); 
%         rfSmooth = imgaussfilt(rfRaw, 1);
% 
%         % Coordinate mapping (based on your 12x8 grid info)
%         % Azimuth ~ -70 to 20, Elevation ~ -20 to 40
%         az = linspace(-70, 20, size(rfRaw, 2));
%         el = linspace(-20, 40, size(rfRaw, 1));
% 
%         imagesc(az, el, rfSmooth);
%         set(gca, 'YDir', 'normal');
%         colormap(subplot(1,3,3), 'parula'); % Distinct color for RF
%         xline(0, 'k:'); yline(0, 'k:');
%         xlabel('Azimuth'); ylabel('Elevation');
%         title('SVD Receptive Field');
%         colorbar;
%     else
%         text(0.5, 0.5, 'No RF Data', 'HorizontalAlignment', 'center');
%         axis off;
%     end
% 
%     exportgraphics(fig, pdfPath, 'Append', true);
%     close(fig);
% end
% end

function plotAllNeuronSummariesToPDF_SparseNoiseAndVR(sessionFileInfo, response, applySmoothing)
% UPDATED: Receptive Field elevation corrected (Row 1 = Bottom / -20 deg)
if nargin < 3, applySmoothing = true; end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pdfPath = fullfile(figSaveDir, sprintf('%s_%s_AllROIs_Summary_withRF_dFF_fliptest.pdf', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name));

%% Load RF data if available
rfDataAvailable = false;
if isfield(sessionFileInfo.otherSessFilePaths, 'sessionROIData') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file')
    vars = who('-file', sessionFileInfo.otherSessFilePaths.sessionROIData);
    if ismember('sparseNoiseRF', vars)
        load(sessionFileInfo.otherSessFilePaths.sessionROIData, 'sparseNoiseRF');
        rfDataAvailable = true;
    end
end

lapActivityFull = response.lapPositionActivity.spks;
nROIs = size(lapActivityFull, 1);

for neuronIdx = 1:nROIs
    roiActivity = squeeze(lapActivityFull(neuronIdx, :, :));
    if all(isnan(roiActivity), 'all'), continue; end
    
    % Smoothing Logic
    if applySmoothing
        w = gausswin(10); w = w / sum(w);
        for iLap = 1:size(roiActivity, 1)
            trace = roiActivity(iLap, :);
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace); smoothed(nanMask) = NaN;
            roiActivity(iLap, :) = smoothed;
        end
    end
    
    meanActivity = mean(roiActivity, 1, 'omitnan');
    semActivity = std(roiActivity, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiActivity), 1));
    normLapActivity = normalize(roiActivity, 2, 'range');
    
    % Create figure
    fig = figure('Visible', 'off', 'Position', [50 100 1600 450], 'Color', 'w');
    
    % --- Subplot 1: Mean ± SEM Trace ---
    subplot(1, 3, 1); hold on;
    x = 1:size(meanActivity, 2);
    fill([x fliplr(x)], [meanActivity + semActivity, fliplr(meanActivity - semActivity)], [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x, meanActivity, 'k', 'LineWidth', 2);
    for p = [40 80 120 160], xline(p, 'k--', 'Alpha', 0.3); end
    xticks([1 40 80 120 160 200]); xlabel('Position (cm)'); ylabel('dFF');
    title(sprintf('ROI %d Mean Activity', neuronIdx));
    
    % --- Subplot 2: Heatmap ---
    subplot(1, 3, 2);
    imagesc(normLapActivity); caxis([0 1]); colormap(subplot(1,3,2), flipud(gray));
    for p = [40 80 120 160], xline(p, 'k--', 'Alpha', 0.3); end
    xticks([1 40 80 120 160 200]); xlabel('Position (cm)'); ylabel('Lap #');
    title('Lap-Position-Activity');
    colorbar; ylabel(colorbar, 'dFF (norm)');
    
    % --- Subplot 3: Receptive Field (SNT) ---
    subplot(1, 3, 3);
    if rfDataAvailable && ~isempty(sparseNoiseRF.initMap{neuronIdx})
        % Extract RF (SVD map from contrast channel 4)
        rfRaw = sparseNoiseRF.initMap{neuronIdx}(:, :, end, 4); 
        
        % Logic update: Row 1 of rfRaw is already the bottom of the screen.
        % Do not use flipud. Apply smoothing directly.
        rfSmooth = imgaussfilt(rfRaw, 1);
        
        % Coordinate mapping
        az = linspace(-70, 20, size(rfRaw, 2));
        el = linspace(-20, 40, size(rfRaw, 1)); % -20 at Index 1
        
        imagesc(az, el, rfSmooth);
        
        % This ensures the Y-axis value -20 is at the bottom visually
        set(gca, 'YDir', 'normal'); 
        
        colormap(subplot(1,3,3), 'parula');
        xline(0, 'k:', 'Alpha', 0.5); yline(0, 'k:', 'Alpha', 0.5);
        xlabel('Azimuth (°)'); ylabel('Elevation (°)');
        title('SVD Receptive Field');
        colorbar;
    else
        text(0.5, 0.5, 'No RF Data', 'HorizontalAlignment', 'center');
        axis off;
    end
    
    % Export to PDF
    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end
end