function fig = plotSIExampleROI(RegionData, varargin)
    % Plots example ROI(s)' trial-averaged EVEN-LAP response for the
    % landmark identity Selectivity Index (SI), cross-normalized to the
    % ODD-lap range (matching how SI itself is computed in
    % extractSplitMetrics.m, and identical in structure to
    % computeSpatialModulationIdex.m's SMI normalization), smoothed with
    % a Gaussian window, with SI value(s) annotated.
    %
    % When using percentile-based selection, plots the N ROIs (default 5)
    % whose SI values are closest to the target percentile.
    %
    % Usage:
    %   plotSIExampleROI(RSPData);                          % 5 examples nearest median SI
    %   plotSIExampleROI(RSPData, 'Percentile', 75);        % 5 examples nearest 75th pctile
    %   plotSIExampleROI(RSPData, 'SessionROI', [4, 27]);   % single manual example
    %   plotSIExampleROI(RSPData, 'Percentile', 50, 'NumExamples', 3);

    p = inputParser;
    addParameter(p, 'Percentile', 50, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'SessionROI', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'Label', '', @ischar);
    addParameter(p, 'NumExamples', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, varargin{:});
    opts = p.Results;

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

    % --- Compute SI for every valid ROI, keeping session/ROI IDs ---
    allSI = []; allIsGrating = []; allSession = []; allROI = [];

    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
           ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end

        lapActivity = sess.ConditionData.Baseline.LapActivity;
        [numROIsTotal, numLaps, numBins] = size(lapActivity);
        validSearchBins = allowedLandmarkBins(allowedLandmarkBins >= 1 & allowedLandmarkBins <= numBins);

        hasSavedExclusions = isfield(sess, 'SMI') && isfield(sess.SMI, 'ExcludeEdgePeakCells') && ...
                             length(sess.SMI.ExcludeEdgePeakCells) == numROIsTotal;

        roisToAnalyze = sess.FilteredROIs;
        for i = 1:length(roisToAnalyze)
            roiIdx = roisToAnalyze(i);
            if hasSavedExclusions && sess.SMI.ExcludeEdgePeakCells(roiIdx), continue; end

            allLaps = squeeze(lapActivity(roiIdx, :, :));
            numBinsThisROI = size(allLaps, 2);

            smoothedAllLaps = allLaps;
            for iLap = 1:size(allLaps, 1)
                trace = allLaps(iLap, :);
                if all(isnan(trace)), continue; end
                nanMask = isnan(trace); trace(nanMask) = 0;
                smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
                smoothedAllLaps(iLap, :) = smoothed;
            end

            trainTrace = mean(smoothedAllLaps(1:2:end, :), 1, 'omitnan');
            testTrace  = mean(smoothedAllLaps(2:2:end, :), 1, 'omitnan');
            if all(isnan(trainTrace)) || all(isnan(testTrace)), continue; end

            minOdd = min(trainTrace, [], 'omitnan');
            maxOdd = max(trainTrace, [], 'omitnan');
            rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end

            normTrainOdd = (trainTrace - minOdd) ./ rangeOdd;
            normTestEven = (testTrace  - minOdd) ./ rangeOdd;

            [~, maxIdxInSearch] = max(normTrainOdd(validSearchBins));
            prefBin = validSearchBins(maxIdxInSearch);
            [~, closestLandmarkIdx] = min(abs(landmarkCentres - prefBin));
            primaryLandmark = landmarkCentres(closestLandmarkIdx);

            gPeaks = nan(1, length(gratingLandmarks));
            for gi = 1:length(gratingLandmarks)
                thisWindow = (gratingLandmarks(gi) - windowRadius):(gratingLandmarks(gi) + windowRadius);
                thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
                if isempty(thisWindow), continue; end
                gPeaks(gi) = max(max(0, normTestEven(thisWindow)), [], 'omitnan');
            end
            actG = mean(gPeaks, 'omitnan');

            pPeaks = nan(1, length(plaidLandmarks));
            for pi = 1:length(plaidLandmarks)
                thisWindow = (plaidLandmarks(pi) - windowRadius):(plaidLandmarks(pi) + windowRadius);
                thisWindow = thisWindow(thisWindow >= 1 & thisWindow <= numBinsThisROI);
                if isempty(thisWindow), continue; end
                pPeaks(pi) = max(max(0, normTestEven(thisWindow)), [], 'omitnan');
            end
            actP = mean(pPeaks, 'omitnan');

            actPeak = max(0, normTestEven(prefBin));
            isGrating = ismember(primaryLandmark, gratingLandmarks);
            actAlt = actP; if ~isGrating, actAlt = actG; end
            if isnan(actAlt) || (actPeak + actAlt) <= 0, continue; end
            SI = (actPeak - actAlt) / (actPeak + actAlt);
            if isnan(SI), continue; end
            SI = max(0, min(1, SI));

            allSI(end+1)        = SI; %#ok<AGROW>
            allIsGrating(end+1) = isGrating; %#ok<AGROW>
            allSession(end+1)   = s; %#ok<AGROW>
            allROI(end+1)       = roiIdx; %#ok<AGROW>
        end
    end

    if isempty(allSI)
        error('No valid SI data computed from RegionData.');
    end

    landmarkTicks = landmarkCentres;

    if ~isempty(opts.SessionROI)
        sessionNum = opts.SessionROI(1);
        roiIdx     = opts.SessionROI(2);
        matchIdx = find(allSession == sessionNum & allROI == roiIdx, 1);
        if isempty(matchIdx)
            error('Requested Session %d, ROI %d not found among valid candidates.', sessionNum, roiIdx);
        end
        siValue   = allSI(matchIdx);
        isGrating = allIsGrating(matchIdx);

        if isempty(opts.Label)
            titleLabel = 'Example ROI';
        else
            titleLabel = opts.Label;
        end

        fig = figure('Color', 'w', 'Position', [200, 200, 450, 350]);
        plotOneSIROI(RegionData, sessionNum, roiIdx, siValue, isGrating, titleLabel, w_space, landmarkTicks, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor);
        return;
    end

    % --- Percentile mode ---
    targetVal = prctile(allSI, opts.Percentile);
    distances = abs(allSI - targetVal);
    [~, sortedIdx] = sort(distances, 'ascend');

    nExamples = min(opts.NumExamples, length(sortedIdx));
    chosenIdx = sortedIdx(1:nExamples);

    if isempty(opts.Label)
        baseLabel = sprintf('%gth pctile', opts.Percentile);
    else
        baseLabel = opts.Label;
    end

    fig = figure('Color', 'w', 'Position', [100, 200, 220*nExamples, 300]);
    for k = 1:nExamples
        idx = chosenIdx(k);
        sessionNum = allSession(idx);
        roiIdx     = allROI(idx);
        siValue    = allSI(idx);
        isGrating  = allIsGrating(idx);

        subplot(1, nExamples, k);
        plotOneSIROI(RegionData, sessionNum, roiIdx, siValue, isGrating, sprintf('%s (#%d)', baseLabel, k), ...
            w_space, landmarkTicks, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor);

        fprintf('%s example %d: S%d R%d, SI = %.3f (target = %.3f)\n', ...
            baseLabel, k, sessionNum, roiIdx, siValue, targetVal);
    end
    sgtitle(sprintf('%s — %d examples nearest target SI = %.3f', baseLabel, nExamples, targetVal));
end

function plotOneSIROI(RegionData, sessionNum, roiIdx, siValue, isGrating, titleLabel, w_space, landmarkTicks, ...
    landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor)

    sess = RegionData(sessionNum);
    lapActivity = sess.ConditionData.Baseline.LapActivity;
    allLaps = squeeze(lapActivity(roiIdx, :, :));   % laps x bins

    numBins = size(allLaps, 2);
    smoothedAllLaps = allLaps;
    for iLap = 1:size(allLaps, 1)
        trace = allLaps(iLap, :);
        if all(isnan(trace)), continue; end
        nanMask = isnan(trace); trace(nanMask) = 0;
        smoothed = filtfilt(w_space, 1, trace); smoothed(nanMask) = NaN;
        smoothedAllLaps(iLap, :) = smoothed;
    end

    oddLaps  = smoothedAllLaps(1:2:end, :);
    evenLaps = smoothedAllLaps(2:2:end, :);

    % Cross-normalize using the ODD-lap range, same as SI/SMI computation
    meanOdd = mean(oddLaps, 1, 'omitnan');
    minOdd = min(meanOdd, [], 'omitnan');
    maxOdd = max(meanOdd, [], 'omitnan');
    rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end

    normEvenLaps = (evenLaps - minOdd) ./ rangeOdd;

    fullMean = mean(normEvenLaps, 1, 'omitnan');
    fullSEM  = std(normEvenLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(normEvenLaps), 1));
    xPos = 1:numBins;

    hold on;

    for c = landmarkCentres
        col = plaidColor; if ismember(c, gratingLandmarks), col = gratingColor; end
        xline(c, '-', 'Color', [col, 0.5], 'LineWidth', 1);
    end

    upper = fullMean + fullSEM;
    lower = fullMean - fullSEM;
    fill([xPos, fliplr(xPos)], [upper, fliplr(lower)], [0.7 0.7 0.7], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xPos, fullMean, 'k', 'LineWidth', 1.8);

    typeLabel = 'Plaid-preferring'; if isGrating, typeLabel = 'Grating-preferring'; end
    text(0.05, 0.92, sprintf('SI = %.2f (%s)', siValue, typeLabel), 'Units', 'normalized', ...
        'FontSize', 9, 'FontWeight', 'bold');

    validTicks = landmarkTicks(landmarkTicks >= 1 & landmarkTicks <= numBins);
    xticks(validTicks);
    xlim([1, numBins]);

    xlabel('Position (cm)');
    ylabel('\Delta F/F (norm.)');
    title(sprintf('%s\nS%d R%d', titleLabel, sessionNum, roiIdx), 'FontWeight', 'normal', 'FontSize', 9);
    box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
end