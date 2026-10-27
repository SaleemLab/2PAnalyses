% function plotIdentityDiscriminationExamples(RegionData, regionLabel)
%     % Selects and plots example ROIs illustrating landmark identity
%     % (non-)discrimination: a typical low-SI cell (representative of the
%     % bulk of the population) and a higher-SI cell (representative of the
%     % minority that does show identity preference), using the current
%     % peak-per-occurrence, averaged SI metric.
% 
%     landmarkCentres  = [40, 80, 120, 160];
%     gratingLandmarks = [40, 120];
%     plaidLandmarks   = [80, 160];
%     tolerance = 15;
%     windowRadius = 15;
%     gratingColor = [0.20 0.40 0.80];
%     plaidColor   = [0.90 0.50 0.10];
% 
%     w_space = gausswin(15); w_space = w_space / sum(w_space);
%     allowedLandmarkBins = [];
%     for c = landmarkCentres
%         allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)];
%     end
%     allowedLandmarkBins = unique(allowedLandmarkBins);
% 
%     allCandidates = struct('SI', {}, 'isGrating', {}, 'session', {}, 'roi', {}, ...
%         'mean', {}, 'sem', {}, 'numBins', {});
% 
%     for s = 1:length(RegionData)
%         sess = RegionData(s);
%         if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
%            ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
%             continue;
%         end
% 
%         lapActivity = sess.ConditionData.Baseline.LapActivity;
%         [numROIsTotal, numLaps, numBins] = size(lapActivity);
%         validSearchBins = allowedLandmarkBins(allowedLandmarkBins >= 1 & allowedLandmarkBins <= numBins);
% 
%         smoothedActivity = lapActivity;
%         for iCell = 1:numROIsTotal
%             for iLap = 1:numLaps
%                 trace = squeeze(lapActivity(iCell, iLap, :));
%                 if all(isnan(trace)), continue; end
%                 nanMask = isnan(trace); trace(nanMask) = 0;
%                 smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
%                 smoothedActivity(iCell, iLap, :) = smoothed;
%             end
%         end
% 
%         roisToAnalyze = sess.FilteredROIs;
%         for i = 1:length(roisToAnalyze)
%             roiIdx = roisToAnalyze(i);
% 
%             allLaps = squeeze(smoothedActivity(roiIdx, :, :));
%             if size(allLaps, 1) < 10, continue; end
%             meanOdd  = mean(allLaps(1:2:end, :), 1, 'omitnan');
%             meanEven = mean(allLaps(2:2:end, :), 1, 'omitnan');
%             if all(isnan(meanOdd)) || all(isnan(meanEven)), continue; end
% 
%             minOdd = min(meanOdd, [], 'omitnan');
%             maxOdd = max(meanOdd, [], 'omitnan');
%             rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
%             normOdd  = (meanOdd - minOdd) ./ rangeOdd;
%             normEven = (meanEven - minOdd) ./ rangeOdd;
% 
%             [~, maxIdxInSearch] = max(normOdd(validSearchBins));
%             prefBin = validSearchBins(maxIdxInSearch);
%             [~, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
%             primaryLandmark = landmarkCentres(closestLandmarkIdx);
% 
%             numBinsThisROI = length(normEven);
% 
%             gPeaks = nan(1, length(gratingLandmarks));
%             for gi = 1:length(gratingLandmarks)
%                 thisWindow = (gratingLandmarks(gi) - windowRadius):(gratingLandmarks(gi) + windowRadius);
%                 thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
%                 if isempty(thisWindow), continue; end
%                 gPeaks(gi) = max(max(0, normEven(thisWindow)), [], 'omitnan');
%             end
%             actG = mean(gPeaks, 'omitnan');
% 
%             pPeaks = nan(1, length(plaidLandmarks));
%             for pi = 1:length(plaidLandmarks)
%                 thisWindow = (plaidLandmarks(pi) - windowRadius):(plaidLandmarks(pi) + windowRadius);
%                 thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
%                 if isempty(thisWindow), continue; end
%                 pPeaks(pi) = max(max(0, normEven(thisWindow)), [], 'omitnan');
%             end
%             actP = mean(pPeaks, 'omitnan');
% 
%             actPeak = max(0, normEven(prefBin));
%             isGrating = ismember(primaryLandmark, gratingLandmarks);
%             actAlt = actP; if ~isGrating, actAlt = actG; end
%             if isnan(actAlt) || (actPeak + actAlt) <= 0, continue; end
%             SI = (actPeak - actAlt) / (actPeak + actAlt);
%             if isnan(SI), continue; end
%             SI = max(0, min(1, SI));
% 
%             fullMean = mean(allLaps, 1, 'omitnan');
%             fullSEM  = std(allLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(allLaps), 1));
% 
%             allCandidates(end+1) = struct('SI', SI, 'isGrating', isGrating, ...
%                 'session', s, 'roi', roiIdx, 'mean', fullMean, 'sem', fullSEM, ...
%                 'numBins', numBinsThisROI); %#ok<AGROW>
%         end
%     end
% 
%     if isempty(allCandidates)
%         error('No valid candidate ROIs found.');
%     end
% 
%     allSI = [allCandidates.SI];
%     medianSI = median(allSI);
% 
%     % "Typical" example: closest to the population median SI
%     [~, idxTypical] = min(abs(allSI - medianSI));
% 
%     % Higher-SI example: for contrast, a cell in the upper range but not
%     % the single most extreme outlier (75th percentile of SI)
%     targetHighSI = prctile(allSI, yyyy);
%     [~, idxHigh] = min(abs(allSI - targetHighSI));
% 
%     exTypical = allCandidates(idxTypical);
%     exHigh    = allCandidates(idxHigh);
% 
%     fig = figure('Color', 'w', 'Position', [200, 200, 700, 350]);
% 
%     subplot(1,2,1);
%     plotExample(exTypical, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, ...
%         sprintf('%s: typical (median SI)\nS%d R%d', regionLabel, exTypical.session, exTypical.roi));
% 
%     subplot(1,2,2);
%     plotExample(exHigh, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, ...
%         sprintf('%s: higher SI (90th pct.)\nS%d R%d', regionLabel, exHigh.session, exHigh.roi));
% 
%     fprintf('%s: median SI = %.3f, n = %d\n', regionLabel, medianSI, length(allSI));
%     fprintf('  Typical example: S%d R%d, SI = %.3f\n', exTypical.session, exTypical.roi, exTypical.SI);
%     fprintf('  Higher example:  S%d R%d, SI = %.3f (90th pct. = %.3f)\n', exHigh.session, exHigh.roi, exHigh.SI, targetHighSI);
% 
%     outputDir = sprintf('Z:\\ibn-vision\\USERS\\Sonali\\Figures\\ThesisFigs\\ResultsChapter2-RSP-PostExp\\Section2_Fig3.4.0\\identity_discrimination_examples');
%     if ~exist(outputDir, 'dir'), mkdir(outputDir); end
%     saveFigureFormats(fig, fullfile(outputDir, sprintf('identity_examples_%s', regionLabel)));
% end
% 
% function plotExample(ex, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, titleStr)
%     hold on;
%     xPos = 1:ex.numBins;
% 
%     for c = landmarkCentres
%         col = plaidColor; if ismember(c, gratingLandmarks), col = gratingColor; end
%         xline(c, '-', 'Color', [col, 0.5], 'LineWidth', 1.5);
%     end
% 
%     upper = ex.mean + ex.sem;
%     lower = ex.mean - ex.sem;
%     fill([xPos, fliplr(xPos)], [upper, fliplr(lower)], [0.7 0.7 0.7], ...
%         'EdgeColor', 'none', 'FaceAlpha', 0.5);
%     plot(xPos, ex.mean, 'k', 'LineWidth', 1.8);
% 
%     text(0.05, 0.92, sprintf('SI = %.2f', ex.SI), 'Units', 'normalized', ...
%         'FontSize', 11, 'FontWeight', 'bold');
% 
%     xlabel('Position (cm)');
%     ylabel('\Delta F/F (norm.)');
%     title(titleStr, 'FontWeight', 'normal');
%     box off;
%     if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
% end


function plotIdentityDiscriminationExamples(RegionData, regionLabel)
    % Selects and plots example ROIs at SI = 0 (no discrimination) and
    % SI = 0.5 (moderate discrimination), using the current peak-per-
    % occurrence, averaged SI metric.

    landmarkCentres  = [40, 80, 120, 160];
    gratingLandmarks = [40, 120];
    plaidLandmarks   = [80, 160];
    tolerance = 15;
    windowRadius = 15;
    gratingColor = [0.20 0.40 0.80];
    plaidColor   = [0.90 0.50 0.10];

    w_space = gausswin(15); w_space = w_space / sum(w_space);
    allowedLandmarkBins = [];
    for c = landmarkCentres
        allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)];
    end
    allowedLandmarkBins = unique(allowedLandmarkBins);

    allCandidates = struct('SI', {}, 'isGrating', {}, 'session', {}, 'roi', {}, ...
        'mean', {}, 'sem', {}, 'numBins', {});

    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end

        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [numROIsTotal, numLaps, numBins] = size(lapActivity);
        validSearchBins = allowedLandmarkBins(allowedLandmarkBins >= 1 & allowedLandmarkBins <= numBins);

        smoothedActivity = lapActivity;
        for iCell = 1:numROIsTotal
            for iLap = 1:numLaps
                trace = squeeze(lapActivity(iCell, iLap, :));
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
                smoothedActivity(iCell, iLap, :) = smoothed;
            end
        end

        roisToAnalyze = sess.FilteredROIs;
        for i = 1:length(roisToAnalyze)
            roiIdx = roisToAnalyze(i);

            allLaps = squeeze(smoothedActivity(roiIdx, :, :));
            if size(allLaps, 1) < 10, continue; end
            meanOdd  = mean(allLaps(1:2:end, :), 1, 'omitnan');
            meanEven = mean(allLaps(2:2:end, :), 1, 'omitnan');
            if all(isnan(meanOdd)) || all(isnan(meanEven)), continue; end

            minOdd = min(meanOdd, [], 'omitnan');
            maxOdd = max(meanOdd, [], 'omitnan');
            rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end
            normOdd  = (meanOdd - minOdd) ./ rangeOdd;
            normEven = (meanEven - minOdd) ./ rangeOdd;

            [~, maxIdxInSearch] = max(normOdd(validSearchBins));
            prefBin = validSearchBins(maxIdxInSearch);
            [~, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            primaryLandmark = landmarkCentres(closestLandmarkIdx);

            numBinsThisROI = length(normEven);

            gPeaks = nan(1, length(gratingLandmarks));
            for gi = 1:length(gratingLandmarks)
                thisWindow = (gratingLandmarks(gi) - windowRadius):(gratingLandmarks(gi) + windowRadius);
                thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
                if isempty(thisWindow), continue; end
                gPeaks(gi) = max(max(0, normEven(thisWindow)), [], 'omitnan');
            end
            actG = mean(gPeaks, 'omitnan');

            pPeaks = nan(1, length(plaidLandmarks));
            for pi = 1:length(plaidLandmarks)
                thisWindow = (plaidLandmarks(pi) - windowRadius):(plaidLandmarks(pi) + windowRadius);
                thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
                if isempty(thisWindow), continue; end
                pPeaks(pi) = max(max(0, normEven(thisWindow)), [], 'omitnan');
            end
            actP = mean(pPeaks, 'omitnan');

            actPeak = max(0, normEven(prefBin));
            isGrating = ismember(primaryLandmark, gratingLandmarks);
            actAlt = actP; if ~isGrating, actAlt = actG; end
            if isnan(actAlt) || (actPeak + actAlt) <= 0, continue; end
            SI = (actPeak - actAlt) / (actPeak + actAlt);
            if isnan(SI), continue; end
            SI = max(0, min(1, SI));

            fullMean = mean(allLaps, 1, 'omitnan');
            fullSEM  = std(allLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(allLaps), 1));

            allCandidates(end+1) = struct('SI', SI, 'isGrating', isGrating, ...
                'session', s, 'roi', roiIdx, 'mean', fullMean, 'sem', fullSEM, ...
                'numBins', numBinsThisROI); %#ok<AGROW>
        end
    end

    if isempty(allCandidates)
        error('No valid candidate ROIs found.');
    end

    allSI = [allCandidates.SI];

    % Example near SI = 0 (no discrimination)
    [~, idxZero] = min(abs(allSI - 0));
    % Example near SI = 0.5 (moderate discrimination)
    [~, idxHalf] = min(abs(allSI - 0.5));

    exZero = allCandidates(idxZero);
    exHalf = allCandidates(idxHalf);

    fig = figure('Color', 'w', 'Position', [200, 200, 700, 350]);

    subplot(1,2,1);
    plotExample(exZero, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, ...
        sprintf('%s: SI ~ 0\nS%d R%d', regionLabel, exZero.session, exZero.roi));

    subplot(1,2,2);
    plotExample(exHalf, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, ...
        sprintf('%s: SI ~ 0.5\nS%d R%d', regionLabel, exHalf.session, exHalf.roi));

    fprintf('%s: n = %d total ROIs\n', regionLabel, length(allSI));
    fprintf('  SI~0 example:   S%d R%d, SI = %.3f\n', exZero.session, exZero.roi, exZero.SI);
    fprintf('  SI~0.5 example: S%d R%d, SI = %.3f\n', exHalf.session, exHalf.roi, exHalf.SI);

    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\identity_discrimination_examples';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, sprintf('identity_examples_%s', regionLabel)));
end

function plotExample(ex, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, titleStr)
    hold on;
    xPos = 1:ex.numBins;

    for c = landmarkCentres
        col = plaidColor; if ismember(c, gratingLandmarks), col = gratingColor; end
        xline(c, '-', 'Color', [col, 0.5], 'LineWidth', 1.5);
    end

    upper = ex.mean + ex.sem;
    lower = ex.mean - ex.sem;
    fill([xPos, fliplr(xPos)], [upper, fliplr(lower)], [0.7 0.7 0.7], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xPos, ex.mean, 'k', 'LineWidth', 1.8);

    text(0.05, 0.92, sprintf('SI = %.2f', ex.SI), 'Units', 'normalized', ...
        'FontSize', 11, 'FontWeight', 'bold');

    xlabel('Position (cm)');
    ylabel('\Delta F/F (norm.)');
    title(titleStr, 'FontWeight', 'normal');
    box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
end