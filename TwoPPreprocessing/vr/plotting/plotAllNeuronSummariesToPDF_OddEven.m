function plotAllNeuronSummariesToPDF_OddEven(sessionFileInfo, response, signalToUse, applySmoothing)
%   Plots ROI summaries (odd vs even), sorted by peak location in odd.
%   Line plots show raw dFF; Heatmaps are normalized 0-1.

if nargin < 3; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 4; applySmoothing = true; end

%% Output path
figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'), mkdir(figSaveDir); end
pdfPath = fullfile(figSaveDir, [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_OddEven_ROI_MiniSnakesAndLines.pdf']);

%% Extract activity matrix
lapActivity = response.lapPositionActivity.(signalToUse);

% Optional spatial smoothing
if applySmoothing
    w = gausswin(5); w = w / sum(w);
    for iCell = 1:size(lapActivity, 1)
        for iLap = 1:size(lapActivity, 2)
            trace = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace);
            trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            lapActivity(iCell, iLap, :) = smoothed;
        end
    end
end

%% Split odd and even laps
oddLaps = lapActivity(:, 1:2:end, :);
evenLaps = lapActivity(:, 2:2:end, :);

% Average across laps
meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));

% Normalize ONLY for sorting purposes
normOddForSort = normalize(meanOdd, 2, 'range');

% Sort cells by peak location in odd lap average
[~, peakIdx] = max(normOddForSort, [], 2);
[~, sortIdx] = sort(peakIdx);

%% Plot Per ROI in Sorted Order
fprintf('Generating PDF for %d ROIs (Sorted by Odd Peak)...\n', size(lapActivity,1));

for i = 1:length(sortIdx)
    roiIdx = sortIdx(i);
    
    % Data for this ROI
    roiOddLaps = squeeze(oddLaps(roiIdx, :, :));
    roiEvenLaps = squeeze(evenLaps(roiIdx, :, :));
    
    if all(isnan(roiOddLaps), 'all'), continue; end
    
    % Raw Mean and SEM (not normalized)
    mO = meanOdd(roiIdx, :);
    mE = meanEven(roiIdx, :);
    
    semO = std(roiOddLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiOddLaps), 1));
    semE = std(roiEvenLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(roiEvenLaps), 1));

    % Create Figure
    fig = figure('Visible', 'off', 'Position', [100 100 1400 400]);
    x = 1:size(lapActivity, 3);

    % 
    subplot(1, 3, 1); hold on;
    % Odd
    fill([x fliplr(x)], [mO + semO, fliplr(mO - semO)], [0.2 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(x, mO, 'b', 'LineWidth', 2);
    % Even
    fill([x fliplr(x)], [mE + semE, fliplr(mE - semE)], [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(x, mE, 'r', 'LineWidth', 2);
    
    title(sprintf('ROI %d (dFF)', roiIdx));
    ylabel(signalToUse); 
    legend({'Odd', 'Even'}, 'Location', 'best');
    set(gca, 'FontSize', 10);

    % --- Odd Heatmap (Normalized 0-1) ---
    subplot(1, 3, 2);
    imagesc(normalize(roiOddLaps, 2, 'range'));
    title('Odd Laps');

    % --- Even Heatmap (Normalized 0-1) ---
    subplot(1, 3, 3);
    imagesc(normalize(roiEvenLaps, 2, 'range'));
    title('Even Laps');

    % Formatting
    for s = 1:3
        ax = subplot(1, 3, s);
        colormap(ax, flipud(gray));
        xticks(ax, [1 50 70 90 110 140]);
        xticklabels(ax, {'1', '50', '70', '90', '110', '140'});
        xline(50, 'k--'); xline(70, 'k--'); xline(90, 'k--'); xline(110, 'k--');
        xlabel(ax, 'Position (cm)');
        set(ax, 'TickDir', 'out', 'box', 'off');
        if s > 1
            set(ax, 'CLim', [0 1], 'YDir', 'normal');
            ylabel(ax, 'Lap #');
            colorbar;
        end
    end

    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end

fprintf('Done. Saved to: %s\n', pdfPath);
end