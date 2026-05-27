function plotOddEvenTrialsVR_IncAndExcFlaggedTrials(sessionFileInfo, response, applySmoothing, signalToUse)

if nargin < 3 || isempty(applySmoothing); applySmoothing = true; end
if nargin < 4 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end

figSaveDir = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
pdfPath = fullfile(figSaveDir, sprintf('%s_%s_FlaggedTrialExclusionSummary_%s.pdf', ...
    sessionFileInfo.animal_name, sessionFileInfo.session_name, signalToUse));
if exist(pdfPath, 'file'); delete(pdfPath); end

lapActivityFull = response.lapPositionActivity.(signalToUse);
[nROIs, nLaps, nBins] = size(lapActivityFull);

flaggedLaps = [];
if isfield(response, 'flaggedLaps')
    flaggedLaps = response.flaggedLaps;
end

unflaggedLapsMask = true(1, nLaps);
unflaggedLapsMask(flaggedLaps) = false;
cleanLapsIdx = unflaggedLapsMask;

for neuronIdx = 1:nROIs
    roiActivityRaw = squeeze(lapActivityFull(neuronIdx, :, :));
    if all(isnan(roiActivityRaw), 'all'); continue; end
    
    roiActivity = roiActivityRaw;
    if applySmoothing
        w_space = gausswin(10); w_space = w_space / sum(w_space);
        for iL = 1:nLaps
            trace = roiActivity(iL, :);
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w_space, 1, trace);
            smoothed(nanMask) = NaN;
            roiActivity(iL, :) = smoothed;
        end
    end
    
    oddIdxAll  = 1:2:nLaps;
    evenIdxAll = 2:2:nLaps;
    
    oddLapsAll  = roiActivity(oddIdxAll, :);
    evenLapsAll = roiActivity(evenIdxAll, :);
    
    mOddAll  = mean(oddLapsAll, 1, 'omitnan');
    sOddAll  = std(oddLapsAll, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(oddLapsAll), 1));
    mEvenAll = mean(evenLapsAll, 1, 'omitnan');
    sEvenAll = std(evenLapsAll, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(evenLapsAll), 1));
    
    [~, flaggedInOdd]  = intersect(oddIdxAll, flaggedLaps);
    [~, flaggedInEven] = intersect(evenIdxAll, flaggedLaps);
    
    roiActivityClean = roiActivity(cleanLapsIdx, :);
    nCleanLaps = size(roiActivityClean, 1);
    
    oddIdxClean  = 1:2:nCleanLaps;
    evenIdxClean = 2:2:nCleanLaps;
    
    oddLapsClean  = roiActivityClean(oddIdxClean, :);
    evenLapsClean = roiActivityClean(evenIdxClean, :);
    
    mOddClean  = mean(oddLapsClean, 1, 'omitnan');
    sOddClean  = std(oddLapsClean, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(oddLapsClean), 1));
    mEvenClean = mean(evenLapsClean, 1, 'omitnan');
    sEvenClean = std(evenLapsClean, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(evenLapsClean), 1));
    
    fig = figure('Visible', 'off', 'Position', [50 50 1400 850]);
    xBins = 1:nBins;
    
    subplot(2, 2, 1);
    imagesc(normalize(oddLapsAll, 2, 'range'));
    colormap(flipud(gray));
    title(sprintf('ROI %d - Odd Laps (All)', neuronIdx));
    ylabel('Odd Lap Index');
    hold on; xline([40 80 120 160], 'k--');
    if ~isempty(flaggedInOdd)
        scatter(repmat(nBins, size(flaggedInOdd)), flaggedInOdd, 40, 'm', 'filled');
    end
    
    subplot(2, 2, 2);
    imagesc(normalize(evenLapsAll, 2, 'range'));
    colormap(flipud(gray));
    title('Even Laps (All)');
    ylabel('Even Lap Index');
    hold on; xline([40 80 120 160], 'k--');
    if ~isempty(flaggedInEven)
        scatter(repmat(nBins, size(flaggedInEven)), flaggedInEven, 40, 'm', 'filled');
    end
    
    subplot(2, 2, 3); hold on;
    fill([xBins fliplr(xBins)], [mOddAll+sOddAll, fliplr(mOddAll-sOddAll)], [0.85 0.33 0.1], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([xBins fliplr(xBins)], [mEvenAll+sEvenAll, fliplr(mEvenAll-sEvenAll)], [0 0.45 0.74], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(xBins, mOddAll, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'DisplayName', 'Odd');
    plot(xBins, mEvenAll, 'Color', [0 0.45 0.74], 'LineWidth', 2, 'DisplayName', 'Even');

    title(sprintf('All laps included (Odd: %d, Even: %d)', size(oddIdxAll,1), size(evenIdxAll,1)));
    ylabel('\DeltaF/F');
    xline([40 80 120 160], 'k--');
    
    
    subplot(2, 2, 4); hold on;
    fill([xBins fliplr(xBins)], [mOddClean+sOddClean, fliplr(mOddClean-sOddClean)], [0.85 0.33 0.1], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([xBins fliplr(xBins)], [mEvenClean+sEvenClean, fliplr(mEvenClean-sEvenClean)], [0 0.45 0.74], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(xBins, mOddClean, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'DisplayName', 'Odd');
    plot(xBins, mEvenClean, 'Color', [0 0.45 0.74], 'LineWidth', 2, 'DisplayName', 'Even');

    title(sprintf('Mean Profile: Clean Laps Only (Odd: %d, Even: %d)', size(oddLapsClean,1), size(evenLapsClean,1)));
    xlabel('Position (cm)'); ylabel('\DeltaF/F'); 
    legend();
    xline([40 80 120 160], 'k--', 'HandleVisibility', 'off');
    
    
    exportgraphics(fig, pdfPath, 'Append', true);
    close(fig);
end

fprintf('PDF Generated Successfully: %s\n', pdfPath);
end