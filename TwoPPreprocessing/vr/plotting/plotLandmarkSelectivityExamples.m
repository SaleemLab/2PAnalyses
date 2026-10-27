function plotLandmarkSelectivityExamples(RSPData, varargin)
    % Selects and plots example boutons illustrating strong grating-,
    % strong plaid-, and low-selectivity responses, in the style of
    % Fischer et al. 2020 Figure 2G.
    %
    % Optional name-value pairs to manually specify examples instead of
    % auto-selecting by SI:
    %   'GratingSessionROI', [sessionNum, roiIdx]
    %   'PlaidSessionROI',   [sessionNum, roiIdx]
    %   'LowSessionROI',     [sessionNum, roiIdx]
    %
    % Example usage:
    %   plotLandmarkSelectivityExamples(RSPData, 'GratingSessionROI', [4, 27]);

    p = inputParser;
    addParameter(p, 'GratingSessionROI', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'PlaidSessionROI',   [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'LowSessionROI',     [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    parse(p, varargin{:});
    manualPicks = p.Results;

    landmarkCentres  = [40, 80, 120, 160];
    gratingLandmarks = [40, 120];
    plaidLandmarks   = [80, 160];
    tolerance = 15;
    windowRadius = 15;

    gratingColor = [0.20 0.40 0.80];  % blue
    plaidColor   = [0.90 0.50 0.10];  % orange

    w_space = gausswin(15); w_space = w_space / sum(w_space);

    gratingBins = []; for g = gratingLandmarks, gratingBins = [gratingBins, (g - windowRadius):(g + windowRadius)]; end
    plaidBins   = []; for pl = plaidLandmarks,  plaidBins   = [plaidBins,   (pl - windowRadius):(pl + windowRadius)]; end
    allowedLandmarkBins = [];
    for c = landmarkCentres
        allowedLandmarkBins = [allowedLandmarkBins, (c - tolerance):(c + tolerance)];
    end
    allowedLandmarkBins = unique(allowedLandmarkBins);

    % --- Pass 1: collect every candidate (needed even if manually picking,
    % since we still need each ROI's SI/mean/sem computed the same way) ---
    allCandidates = struct('SI', {}, 'isGrating', {}, 'session', {}, 'roi', {}, ...
        'mean', {}, 'sem', {}, 'numBins', {});

    for s = 1:length(RSPData)
        sess = RSPData(s);
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

        % NOTE: when manually specifying ROIs, we still need to loop over
        % all ROIs in this session (not just FilteredROIs) in case the
        % manually requested ROI didn't pass the standard filtering.
        roisToAnalyze = 1:numROIsTotal;
        for i = 1:length(roisToAnalyze)
            roiIdx = roisToAnalyze(i);

            allLaps = squeeze(smoothedActivity(roiIdx, :, :));   % laps x bins
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

            actG = mean(max(0, normEven(gratingBins)), 'omitnan');
            actP = mean(max(0, normEven(plaidBins)), 'omitnan');
            actPeak = max(0, normEven(prefBin));

            isGrating = ismember(primaryLandmark, gratingLandmarks);
            actAlt = actP; if ~isGrating, actAlt = actG; end
            if (actPeak + actAlt) <= 0, continue; end
            SI = (actPeak - actAlt) / (actPeak + actAlt);
            if isnan(SI), continue; end
            SI = max(0, min(1, SI));

            fullMean = mean(allLaps, 1, 'omitnan');
            fullSEM  = std(allLaps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(allLaps), 1));

            allCandidates(end+1) = struct('SI', SI, 'isGrating', isGrating, ...
                'session', s, 'roi', roiIdx, 'mean', fullMean, 'sem', fullSEM, ...
                'numBins', numBins); %#ok<AGROW>
        end
    end

    if isempty(allCandidates)
        error('No valid candidate boutons found.');
    end

    sessionArr = [allCandidates.session];
    roiArr     = [allCandidates.roi];
    allSI      = [allCandidates.SI];
    isGratingArr = [allCandidates.isGrating];

    % --- Selection: manual pick if given, else auto-select by SI ---
    exGratingHigh = selectCandidate(allCandidates, sessionArr, roiArr, allSI, isGratingArr, ...
        manualPicks.GratingSessionROI, 'max', true);
    exPlaidHigh = selectCandidate(allCandidates, sessionArr, roiArr, allSI, isGratingArr, ...
        manualPicks.PlaidSessionROI, 'max', false);
    exLow = selectCandidate(allCandidates, sessionArr, roiArr, allSI, isGratingArr, ...
        manualPicks.LowSessionROI, 'min', []);

    % --- Plotting ---
    fig = figure('Color', 'w', 'Position', [200, 200, 1250, 350]);

    subplot(1,3,1);
    plotExampleTrace(exGratingHigh, landmarkCentres, gratingLandmarks, plaidLandmarks, ...
        gratingColor, plaidColor, sprintf('Grating-selective (Session %d, ROI %d)', exGratingHigh.session, exGratingHigh.roi));

    subplot(1,3,2);
    plotExampleTrace(exPlaidHigh, landmarkCentres, gratingLandmarks, plaidLandmarks, ...
        gratingColor, plaidColor, sprintf('Plaid-selective (Session %d, ROI %d)', exPlaidHigh.session, exPlaidHigh.roi));

    subplot(1,3,3);
    plotExampleTrace(exLow, landmarkCentres, gratingLandmarks, plaidLandmarks, ...
        gratingColor, plaidColor, sprintf('Less selective (Session %d, ROI %d)', exLow.session, exLow.roi));

    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\landmark_selectivity_examples';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, 'example_grating_plaid_lowSI_boutons'));
end

function ex = selectCandidate(allCandidates, sessionArr, roiArr, allSI, isGratingArr, manualPick, extremeType, gratingFilter)
    if ~isempty(manualPick)
        matchIdx = find(sessionArr == manualPick(1) & roiArr == manualPick(2), 1);
        if isempty(matchIdx)
            error('Requested Session %d, ROI %d not found among valid candidates (check it passed the >=10 lap filter).', ...
                manualPick(1), manualPick(2));
        end
        ex = allCandidates(matchIdx);
        return;
    end

    if isempty(gratingFilter)
        mask = true(size(allSI));   % pool both types (used for 'low' selection)
    else
        mask = isGratingArr == gratingFilter;
    end

    scores = allSI;
    scores(~mask) = NaN;
    if strcmp(extremeType, 'max')
        [~, idx] = max(scores);
    else
        [~, idx] = min(scores);
    end
    ex = allCandidates(idx);
end

function plotExampleTrace(ex, landmarkCentres, gratingLandmarks, plaidLandmarks, gratingColor, plaidColor, titleStr)
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