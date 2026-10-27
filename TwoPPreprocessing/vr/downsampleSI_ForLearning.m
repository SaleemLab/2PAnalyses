function downsampledResults = downsampleSI_ForLearning(RSPDataAcrossDays, varargin)
% downsampleSI_ForLearning: Matches lap counts across days/sessions before
% computing Skaggs spatial information (SI), corrected against a shuffled
% null, so SI can be fairly compared across days despite differing lap
% counts. Mirrors downsampleCVR2_ForLearning_Fast in structure: repeated
% random lap subsampling, pooled null across repeats, one p-value/z-score
% per cell from the averaged true SI vs. the pooled null.
%
% Shuffle is generated fresh via circular shift of each lap's tuning
% curve (same approach as the EV null), since no precomputed shuffle
% matrix is available in RSPDataAcrossDays.

    p = inputParser;
    addRequired(p, 'RSPDataAcrossDays', @isstruct);
    addParameter(p, 'nShufflesPerRepeat', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1)==0);
    addParameter(p, 'nRepeats', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1)==0);
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1)==0);
    parse(p, RSPDataAcrossDays, varargin{:});

    nShuffles  = p.Results.nShufflesPerRepeat;
    nRepeats   = p.Results.nRepeats;
    applySmooth = p.Results.ApplySmoothing;

    rng(p.Results.RandomSeed);

    allLaps = arrayfun(@(s) s.ConditionData.Baseline.NumLaps, RSPDataAcrossDays);
    targetNumLaps = min(allLaps);
    fprintf('Downsampling every session to %d laps (min across all sessions/days) for SI.\n', targetNumLaps);

    if applySmooth
        w = gausswin(15); w = w / sum(w);
    end

    % Skaggs SI formula, same as getSpatialInfo_WithShuffle
    calcSI = @(tc) sum((1/length(tc)) .* (tc ./ (mean(tc) + eps)) .* log2((tc ./ (mean(tc) + eps)) + eps));

    nSessions = length(RSPDataAcrossDays);
    downsampledResults = struct();
    totalTicStart = tic;

    for s = 1:nSessions
        sessTic = tic;
        thisSession = RSPDataAcrossDays(s);
        lapActivity = thisSession.ConditionData.Baseline.LapActivity; % [nROIs x nLaps x nPositions]
        [nROIs, numLaps, numPositions] = size(lapActivity);

        if numLaps < targetNumLaps
            warning('Mouse %s Day %d has fewer laps (%d) than target (%d) -- skipping.', ...
                thisSession.MouseID, thisSession.Day, numLaps, targetNumLaps);
            downsampledResults(s).meanSI = [];
            downsampledResults(s).zSI = [];
            downsampledResults(s).isSignificant = [];
            continue;
        end

        % Smoothing, matched to the EV pipeline
        if applySmooth
            for iCell = 1:nROIs
                for iLap = 1:numLaps
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

        rawSI_reps         = nan(nROIs, nRepeats);
        pooledNullSI_perCell = nan(nROIs, nRepeats * nShuffles);

        for r = 1:nRepeats
            subLapIdx  = randperm(numLaps, targetNumLaps);
            subActivity = lapActivity(:, subLapIdx, :); % [nROIs x targetNumLaps x numPositions]

            % Mean tuning curve per cell for this repeat's subsample
            meanTuning = squeeze(mean(subActivity, 2, 'omitnan')); % [nROIs x numPositions]
            if nROIs == 1, meanTuning = meanTuning'; end
            meanTuning(meanTuning < 0) = 0; % negative values (z-scored data) clipped, as in original function

            nullVals_thisRep = nan(nROIs, nShuffles);

            for iCell = 1:nROIs
                tc = meanTuning(iCell, :);
                if all(isnan(tc)) || all(tc == 0), continue; end

                rawSI_reps(iCell, r) = calcSI(tc);

                % CHANGED: vectorized shuffle across all shuffles at once
                % for this cell's mean tuning curve (same circular-shift
                % logic as the EV null, applied here to the already-
                % averaged tuning curve rather than per-lap, since SI is
                % computed on the mean curve, not per-lap).
                cellActivity = squeeze(subActivity(iCell, :, :)); % [targetNumLaps x numPositions]
                validLapMask = ~all(isnan(cellActivity), 2);
                randomShifts = randi(numPositions, targetNumLaps, nShuffles);

                shuffledStack = nan(targetNumLaps, numPositions, nShuffles);
                for iLap = 1:targetNumLaps
                    if ~validLapMask(iLap)
                        shuffledStack(iLap, :, :) = repmat(cellActivity(iLap, :)', 1, nShuffles);
                        continue;
                    end
                    lapData = cellActivity(iLap, :);
                    shiftAmounts = randomShifts(iLap, :);
                    idxBase = (0:numPositions-1)';
                    shiftedIdx = mod(idxBase + shiftAmounts - 1, numPositions) + 1;
                    shuffledStack(iLap, :, :) = reshape(lapData(shiftedIdx), 1, numPositions, nShuffles);
                end

                for thisShuff = 1:nShuffles
                    shuffledMeanTuning = mean(shuffledStack(:, :, thisShuff), 1, 'omitnan');
                    shuffledMeanTuning(shuffledMeanTuning < 0) = 0;
                    nullVals_thisRep(iCell, thisShuff) = calcSI(shuffledMeanTuning);
                end
            end

            nullCols = (r-1)*nShuffles + (1:nShuffles);
            pooledNullSI_perCell(:, nullCols) = nullVals_thisRep;
        end

        meanSI_avg = mean(rawSI_reps, 2, 'omitnan');
        mu_s    = mean(pooledNullSI_perCell, 2, 'omitnan');
        sigma_s = std(pooledNullSI_perCell, 0, 2, 'omitnan');
        zSI = (meanSI_avg - mu_s) ./ (sigma_s + eps);

        siThreshold95 = prctile(pooledNullSI_perCell, 95, 2);
        isSISignificant = meanSI_avg > siThreshold95;

        downsampledResults(s).meanSI        = meanSI_avg;
        downsampledResults(s).zSI           = zSI;
        downsampledResults(s).isSignificant = isSISignificant;
        downsampledResults(s).targetNumLaps = targetNumLaps;
        downsampledResults(s).nShufflesTotal = nRepeats * nShuffles;
        downsampledResults(s).nRepeats      = nRepeats;

        sessElapsed = toc(sessTic);
        fprintf('Mouse %s Day %d (%d ROIs): SI done in %.1f sec. [%d/%d sessions]\n', ...
            thisSession.MouseID, thisSession.Day, nROIs, sessElapsed, s, nSessions);
    end

    fprintf('Total elapsed: %.1f min\n', toc(totalTicStart)/60);
end