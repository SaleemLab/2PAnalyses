function [r, p, stableIdx] = checkOddEvenCorrelation(sessionFileInfo, response, signalToUse, applySpatialSmoothing, plotFlag)
% Calculates the odd vs. even laps spatial correlation for
% each RO
%
%   r         - [nROIs x 1] vector of Pearson correlation coefficients (r-value)
%   p         - [nROIs x 1] vector of p-values for each correlation
%   stableIdx - [nStableROIs x 1] indices of ROIs with r > 0.5

%% Handle optional inputs
if nargin < 3; signalToUse = 'dFF'; end % changed dff jan 2026 
if nargin < 4; applySpatialSmoothing = true; end 
if nargin < 5; plotFlag = true; end

%% Save figure save path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

filename = fullfile(figSaveDir, ...
    [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' signalToUse '_oddEvenCorr_SortedbyOdd.png']);

%% Get data 
% lapPositionActivity is (ROI x Laps x Position)
lapPositionActivity = response.lapPositionActivity.(signalToUse);

%% Optional spatial smoothning before computing correlations
if applySpatialSmoothing
    fprintf('Applying spatial smoothing...\n');
    w = gausswin(5); % 10-bin Gaussian window
    w = w / sum(w);
    for iCell = 1:size(lapPositionActivity, 1)
        for iLap = 1:size(lapPositionActivity, 2)
            trace = squeeze(lapPositionActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            
            % Handle NaNs for filtering
            nanMask = isnan(trace);
            trace(nanMask) = 0; % Temporarily set NaNs to 0 for filtfilt
            
            smoothed = filtfilt(w, 1, trace);
            
            smoothed(nanMask) = NaN; % Restore NaNs
            lapPositionActivity(iCell, iLap, :) = smoothed;
        end
    end
end

%% Split odd and even laps
oddLaps = lapPositionActivity(:, 1:2:end, :); 
evenLaps = lapPositionActivity(:, 2:2:end, :);

%%  Average tuning curves for each half
% Results (meanOdd, meanEven) are (ROI x Position)
meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));

%% Normalise tuning curves
% This scales each cell's tuning curve to the [0, 1] range.
normOdd = normalize(meanOdd, 2, 'range');
normEven = normalize(meanEven, 2, 'range');

%% Calculate Pearson correlation
% Correlate the (Position x ROI) matrices.
% correlation between cell i (odd) and cell j (even)
[r_matrix, p_matrix] = corr(normOdd', normEven', 'rows', 'pairwise');

% We only care about the diagonal: corr(cell_i_odd, cell_i_even)
r = diag(r_matrix);
p = diag(p_matrix);

%% Identify "stable" cells (0.6 threshold @aman)
stableThresh = 0.4;
stableIdx = find(r > stableThresh);

%%
lapCorr_OddEven.rho = r;
lapCorr_OddEven.p = p;
lapCorr_OddEven.stableThreshold = stableThresh;
lapCorr_OddEven.stableIdx = stableIdx;

%% Save 
% Check if the file path exists and save the variables
if isfield(sessionFileInfo, 'otherSessFilePaths') && exist(sessionFileInfo.otherSessFilePaths.sessionROIData, 'file') == 2
    
    disp(['Saving Odd-Even Lap significance results to: ', sessionFileInfo.otherSessFilePaths.sessionROIData]);
    
    save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
       "lapCorr_OddEven", ...
         '-append')
         
elseif isfield(sessionFileInfo, 'otherSessFilePaths')
    warning('sessionROIData file not found at: %s. Cannot append peak significance data.', ...
        sessionFileInfo.otherSessFilePaths.sessionROIData);
else
    warning('sessionFileInfo.otherSessFilePaths field not found. Cannot save peak significance data.');
end
% %% Save in response (temp) 
% response.roiOddEvenCorr.r = r;
% response.roiOddEvenCorr.p = p; 

%% Plotting
if plotFlag
    fprintf('Plotting %d stable cells (r > %.2f)...\n', length(stableIdx), stableThresh);

    
    
    % Get the tuning curves for stable cells
    normOddStable = normOdd(stableIdx, :);
    normEvenStable = normEven(stableIdx, :);
    
    % Sort them by the peak of the odd-lap tuning curve
    [~, peakIdx] = max(normOddStable, [], 2);
    [~, sortIdx] = sort(peakIdx);
    
    fig1 = figure('Name', 'Odd-Even Lap Correlations');
    
    % Odd Laps
    subplot(121)
    imagesc(normOddStable(sortIdx,:));
    caxis([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    title(sprintf('Odd Laps (n=%d)', size(oddLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROI (Sorted)');
    
    % Even Laps
    subplot(122)
    imagesc(normEvenStable(sortIdx,:));
    caxis([0 1]); colormap(flipud(gray));
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 12, 'YDir', 'normal');
    xline(50, 'k--', 'LineWidth', 1.5);
    xline(70, 'k--', 'LineWidth', 1.5);
    xline(90, 'k--', 'LineWidth', 1.5);
    xline(110, 'k--', 'LineWidth', 1.5);
    xticks([0 50 70 90 110 140]);
    xticklabels({'0', '50', '70', '90', '110', '140'});
    title(sprintf('Even Laps (n=%d)', size(evenLaps, 2)));
    xlabel('Position (cm)');
    ylabel('Stable ROIs (Sorted)');
    

    %% Save
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperPosition', [0 0 11 8.5], ...
             'PaperOrientation', 'landscape');
    print(gcf, filename, '-dpng', '-r300');
    
end

end