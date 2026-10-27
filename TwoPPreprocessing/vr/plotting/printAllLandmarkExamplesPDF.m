function printAllLandmarkExamplesPDF(RegionData, outputPath, regionLabel)
    % Prints every candidate ROI's tuning curve to a multi-page PDF,
    % using the CURRENT (peak-per-occurrence, averaged) SI calculation,
    % with session/ROI ID and SI annotated on each panel.

    if nargin < 3, regionLabel = 'Region'; end
    if nargin < 2 || isempty(outputPath)
        outputPath = sprintf('Z:\\ibn-vision\\USERS\\Sonali\\Figures\\ThesisFigs\\ResultsChapter2-RSP-PostExp\\Section2_Fig3.4.0\\landmark_selectivity_examples\\all_candidates_%s.pdf', regionLabel);
    end
    if exist(outputPath, 'file'), delete(outputPath); end

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

    [~, sortOrder] = sort([allCandidates.SI], 'descend');
    allCandidates = allCandidates(sortOrder);

    nPerPage = 20; nRows = 4; nCols = 5;
    nTotal = length(allCandidates);
    nPages = ceil(nTotal / nPerPage);

    fprintf('Printing %d candidates (%s) across %d pages to:\n%s\n', nTotal, regionLabel, nPages, outputPath);

    for pg = 1:nPages
        fig = figure('Color', 'w', 'Position', [50 50 1500 1100], 'Visible', 'off');
        idxStart = (pg-1)*nPerPage + 1;
        idxEnd = min(pg*nPerPage, nTotal);

        for k = idxStart:idxEnd
            ex = allCandidates(k);
            subplot(nRows, nCols, k - idxStart + 1);
            hold on;

            for c = landmarkCentres
                col = plaidColor; if ismember(c, gratingLandmarks), col = gratingColor; end
                xline(c, '-', 'Color', [col, 0.5], 'LineWidth', 1);
            end

            xPos = 1:ex.numBins;
            upper = ex.mean + ex.sem;
            lower = ex.mean - ex.sem;
            fill([xPos, fliplr(xPos)], [upper, fliplr(lower)], [0.7 0.7 0.7], ...
                'EdgeColor', 'none', 'FaceAlpha', 0.5);
            plot(xPos, ex.mean, 'k', 'LineWidth', 1);

            typeLabel = 'P'; if ex.isGrating, typeLabel = 'G'; end
            title(sprintf('S%d R%d | SI=%.2f (%s)', ex.session, ex.roi, ex.SI, typeLabel), ...
                'FontSize', 8, 'FontWeight', 'normal');
            set(gca, 'FontSize', 6, 'XTick', [], 'YTick', []);
            box off;
        end

        sgtitle(sprintf('%s — Page %d of %d (sorted by SI, descending, current peak-averaged metric)', regionLabel, pg, nPages));

        if pg == 1
            exportgraphics(fig, outputPath, 'ContentType', 'vector');
        else
            exportgraphics(fig, outputPath, 'ContentType', 'vector', 'Append', true);
        end
        close(fig);
    end

    fprintf('Done.\n');
end