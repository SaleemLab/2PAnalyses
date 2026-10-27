function downsampledResults = downsampleCVR2_ForLearning(RSPDataAcrossDays, varargin)
% downsampleCVR2_ForLearning_Fast: Matches lap counts across days before
% computing cross-validated explained variance (cvEV).
%
% CHANGED (this version):
%   - targetNumLaps = (true min lap count) - LapSlack, so every session,
%     including the one at the original floor, has genuine repeat-to-
%     repeat resampling variability (a session at exactly the target lap
%     count would otherwise draw the identical laps every repeat).
%   - nFolds default changed from 5 to 3, since at ~20 laps a 5-fold
%     split leaves only ~4 laps per test fold, a noisy estimate; 3-fold
%     leaves ~7 laps per test fold, more stable, at some cost to
%     training-set size.
%   - nRepeats default increased from 50 to 75, since more repeats
%     directly tightens the averaged EV estimate (reduces noise from
%     which laps got subsampled); this is independent of the fold-size
%     fix and complements it.

    p = inputParser;
    addRequired(p, 'RSPDataAcrossDays', @isstruct);
    addParameter(p, 'nFolds', 3, @(x) isnumeric(x) && isscalar(x) && x >= 2 && mod(x,1)==0);
    addParameter(p, 'nShufflesPerRepeat', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1)==0);
    addParameter(p, 'nRepeats', 75, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1)==0);
    addParameter(p, 'ApplySmoothing', true, @islogical);
    addParameter(p, 'RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1)==0);
    addParameter(p, 'LapSlack', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1)==0);
    parse(p, RSPDataAcrossDays, varargin{:});

    nFolds     = p.Results.nFolds;
    nShuffles  = p.Results.nShufflesPerRepeat;
    nRepeats   = p.Results.nRepeats;
    applySmooth = p.Results.ApplySmoothing;
    lapSlack   = p.Results.LapSlack;

    rng(p.Results.RandomSeed);

    allLaps = arrayfun(@(s) s.ConditionData.Baseline.NumLaps, RSPDataAcrossDays);
    trueMinLaps = min(allLaps);
    targetNumLaps = trueMinLaps - lapSlack;
    if targetNumLaps < 2
        error('True minimum lap count (%d) minus LapSlack (%d) is too small (<2) to run cross-validation. Reduce LapSlack.', trueMinLaps, lapSlack);
    end
    fprintf('True minimum laps: %d. Downsampling to %d laps (%d lap(s) slack), %d-fold CV, %d repeats.\n', ...
        trueMinLaps, targetNumLaps, lapSlack, nFolds, nRepeats);

    if applySmooth
        w = gausswin(15); w = w / sum(w);
    end

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
            downsampledResults(s).meanExpVar = [];
            downsampledResults(s).pValues = [];
            continue;
        end

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

        meanExpVar_reps   = nan(nROIs, nRepeats);
        pooledNull_perCell = nan(nROIs, nRepeats * nShuffles);

        for r = 1:nRepeats
            subLapIdx  = randperm(numLaps, targetNumLaps);
            subActivity = lapActivity(:, subLapIdx, :);

            cvExpVar_thisRep = nan(nROIs, nFolds);
            nullVals_thisRep = nan(nROIs, nShuffles);

            for iCell = 1:nROIs
                cellActivity = squeeze(subActivity(iCell, :, :));
                if all(isnan(cellActivity), 'all'); continue; end

                cvExpVar_thisRep(iCell, :) = getAllExpVar(cellActivity, nFolds);

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
                    shuffledActivity = shuffledStack(:, :, thisShuff);
                    shuffledFolds = getAllExpVar(shuffledActivity, nFolds);
                    nullVals_thisRep(iCell, thisShuff) = mean(shuffledFolds, 'omitnan');
                end
            end

            meanExpVar_reps(:, r) = mean(cvExpVar_thisRep, 2, 'omitnan');
            nullCols = (r-1)*nShuffles + (1:nShuffles);
            pooledNull_perCell(:, nullCols) = nullVals_thisRep;
        end

        meanExpVar_avg = mean(meanExpVar_reps, 2, 'omitnan');
        countBetterNulls = sum(pooledNull_perCell >= meanExpVar_avg, 2, 'omitnan');
        totalNullCount   = sum(~isnan(pooledNull_perCell), 2);
        totalNullCount(totalNullCount == 0) = 1;
        pValues = (countBetterNulls + 1) ./ (totalNullCount + 1);

        downsampledResults(s).meanExpVar = meanExpVar_avg;
        downsampledResults(s).pValues    = pValues;
        downsampledResults(s).targetNumLaps = targetNumLaps;
        downsampledResults(s).nFolds        = nFolds;
        downsampledResults(s).nShufflesTotal = nRepeats * nShuffles;
        downsampledResults(s).nRepeats      = nRepeats;

        sessElapsed = toc(sessTic);
        fprintf('Mouse %s Day %d (%d ROIs): done in %.1f sec. [%d/%d sessions]\n', ...
            thisSession.MouseID, thisSession.Day, nROIs, sessElapsed, s, nSessions);
    end

    fprintf('Total elapsed: %.1f min\n', toc(totalTicStart)/60);
end