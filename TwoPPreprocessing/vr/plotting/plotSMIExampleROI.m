function fig = plotSMIExampleROI(RegionData, varargin)
    % Plots example ROI(s)' trial-averaged EVEN-LAP response, cross-
    % normalized to the ODD-lap range (matching how SMI itself was
    % computed in computeSpatialModulationIdex.m), smoothed with a
    % Gaussian window, with SMI value(s) annotated.
    %
    % When using percentile-based selection, plots the 5 ROIs whose SMI
    % values are closest to the target percentile (not just the single
    % nearest one), shown as a row of subplots, so a single noisy cell
    % doesn't stand in for the population at that percentile.
    %
    % Usage:
    %   % Auto-select 5 examples nearest the median (50th percentile):
    %   plotSMIExampleROI(RSPData);
    %
    %   % Auto-select 5 examples nearest the 75th percentile:
    %   plotSMIExampleROI(RSPData, 'Percentile', 75);
    %
    %   % Manually specify one exact session/ROI (single plot, no subplot row):
    %   plotSMIExampleROI(RSPData, 'SessionROI', [4, 27]);
    %
    %   % Add a custom title label:
    %   plotSMIExampleROI(RSPData, 'Percentile', 25, 'Label', 'RSP: 25th pctile');
    %
    %   % Change how many examples are shown per percentile call (default 5):
    %   plotSMIExampleROI(RSPData, 'Percentile', 50, 'NumExamples', 3);

    p = inputParser;
    addParameter(p, 'Percentile', 50, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'SessionROI', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'Label', '', @ischar);
    addParameter(p, 'NumExamples', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, varargin{:});
    opts = p.Results;

    w_space = gausswin(15); w_space = w_space / sum(w_space);
    landmarkTicks = [40, 80, 120, 160];

    if ~isempty(opts.SessionROI)
        % Manual single-ROI mode
        sessionNum = opts.SessionROI(1);
        roiIdx     = opts.SessionROI(2);
        smiValue   = RegionData(sessionNum).SMI.SMI(roiIdx);

        if isempty(opts.Label)
            titleLabel = 'Example ROI';
        else
            titleLabel = opts.Label;
        end

        fig = figure('Color', 'w', 'Position', [200, 200, 450, 350]);
        plotOneROI(RegionData, sessionNum, roiIdx, smiValue, titleLabel, w_space, landmarkTicks);
        return;
    end

    % --- Percentile mode: find the N ROIs closest to the target percentile ---
    pooledSMI = [];
    sessionIDArr = [];
    roiIDArr = [];
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMI = sess.SMI.SMI;
            theseROIs = sess.FilteredROIs(:);
            theseSMI = rawSMI(theseROIs);
            validMask = ~isnan(theseSMI);
            pooledSMI = [pooledSMI; theseSMI(validMask)];
            sessionIDArr = [sessionIDArr; repmat(s, sum(validMask), 1)];
            roiIDArr = [roiIDArr; theseROIs(validMask)];
        end
    end

    if isempty(pooledSMI)
        error('No valid SMI data found in RegionData.');
    end

    targetVal = prctile(pooledSMI, opts.Percentile);
    distances = abs(pooledSMI - targetVal);
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
        sessionNum = sessionIDArr(chosenIdx(k));
        roiIdx     = roiIDArr(chosenIdx(k));
        smiValue   = pooledSMI(chosenIdx(k));

        subplot(1, nExamples, k);
        plotOneROI(RegionData, sessionNum, roiIdx, smiValue, sprintf('%s (#%d)', baseLabel, k), w_space, landmarkTicks);

        fprintf('%s example %d: S%d R%d, SMI = %.3f (target = %.3f)\n', ...
            baseLabel, k, sessionNum, roiIdx, smiValue, targetVal);
    end
    sgtitle(sprintf('%s — %d examples nearest target SMI = %.3f', baseLabel, nExamples, targetVal));
end

function plotOneROI(RegionData, sessionNum, roiIdx, smiValue, titleLabel, w_space, landmarkTicks)
    sess = RegionData(sessionNum);
    lapActivity = sess.ConditionData.Baseline.LapActivity;
    allLaps = squeeze(lapActivity(roiIdx, :, :));   % laps x bins

    % Smooth each lap individually (odd AND even both needed: odd sets
    % the normalization range, even is what gets plotted), matching
    % computeSpatialModulationIdex.m
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

    % CHANGED: cross-normalize using the ODD-lap range, matching how SMI
    % itself was computed (odd sets the scale, even is tested on that scale)
    meanOdd = mean(oddLaps, 1, 'omitnan');
    minOdd = min(meanOdd, [], 'omitnan');
    maxOdd = max(meanOdd, [], 'omitnan');
    rangeOdd = maxOdd - minOdd; if rangeOdd == 0, rangeOdd = 1; end

    % Normalize each individual even lap using the odd-lap range, so SEM
    % reflects lap-to-lap variability on the same normalized scale
    normEvenLaps = (evenLaps - minOdd) ./ rangeOdd;

    fullMean = mean(normEvenLaps, 1, 'omitnan');
    fullSEM  = std(normEvenLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(normEvenLaps), 1));
    xPos = 1:numBins;

    hold on;
    upper = fullMean + fullSEM;
    lower = fullMean - fullSEM;
    fill([xPos, fliplr(xPos)], [upper, fliplr(lower)], [0.7 0.7 0.7], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(xPos, fullMean, 'k', 'LineWidth', 1.8);

    text(0.05, 0.92, sprintf('SMI = %.3f', smiValue), 'Units', 'normalized', ...
        'FontSize', 10, 'FontWeight', 'bold');

    validTicks = landmarkTicks(landmarkTicks >= 1 & landmarkTicks <= numBins);
    xticks(validTicks);
    xlim([1, numBins]);

    xlabel('Position (cm)');
    ylabel('\Delta F/F (norm.)');   % CHANGED: reflects normalized scale
    title(sprintf('%s\nS%d R%d', titleLabel, sessionNum, roiIdx), 'FontWeight', 'normal', 'FontSize', 9);
    box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
end