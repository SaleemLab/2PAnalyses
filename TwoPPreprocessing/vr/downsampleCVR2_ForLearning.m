%% Downsampling control: match lap counts across days before computing cvEV
% Mirrors getCrossValidatedExplainedVariance / getAllExpVar, but subsamples
% laps to a common target count first, repeating many times to average out
% subsampling noise.

nFolds     = 5;      % match whatever you used originally
nShuffles  = 200;    % reduced from 1000 for tractability across repeats -- raise if you can afford the compute
nRepeats   = 50;     % number of random lap subsamples per session

% 1) Find the minimum lap count across ALL sessions you're comparing
allLaps = arrayfun(@(s) s.ConditionData.Baseline.NumLaps, RSPDataAcrossDays);
targetNumLaps = min(allLaps);
fprintf('Downsampling every session to %d laps (min across all sessions/days).\n', targetNumLaps);

nSessions = length(RSPDataAcrossDays);
downsampledResults = struct();

for s = 1:nSessions
    thisSession = RSPDataAcrossDays(s);
    lapActivity = thisSession.ConditionData.Baseline.LapActivity; % [nROIs x nLaps x nPositions] -- confirm field name/shape
    [nROIs, numLaps, numPositions] = size(lapActivity);

    if numLaps < targetNumLaps
        warning('Mouse %s Day %d has fewer laps (%d) than target (%d) -- skipping.', ...
            thisSession.MouseID, thisSession.Day, numLaps, targetNumLaps);
        downsampledResults(s).meanExpVar = [];
        downsampledResults(s).pValues = [];
        continue;
    end

    meanExpVar_reps = nan(nROIs, nRepeats);
    pValues_reps    = nan(nROIs, nRepeats);

    for r = 1:nRepeats
        subLapIdx  = randperm(numLaps, targetNumLaps);
        subActivity = lapActivity(:, subLapIdx, :);

        cvExpVar_thisRep     = nan(nROIs, nFolds);
        cvExpVarNull_thisRep = nan(nROIs, nShuffles);

        for iCell = 1:nROIs
            cellActivity = squeeze(subActivity(iCell, :, :)); % [Laps x Position]
            if all(isnan(cellActivity), 'all'); continue; end

            cvExpVar_thisRep(iCell, :) = getAllExpVar(cellActivity, nFolds);

            for thisShuff = 1:nShuffles
                shuffledActivity = nan(targetNumLaps, numPositions);
                rs = RandStream('mt19937ar', 'Seed', thisShuff + r*10000);
                for iLap = 1:targetNumLaps
                    lapData = cellActivity(iLap, :);
                    if ~all(isnan(lapData))
                        randomShift = randi(rs, numPositions);
                        shuffledActivity(iLap, :) = circshift(lapData, randomShift);
                    else
                        shuffledActivity(iLap, :) = lapData;
                    end
                end
                shuffledFolds = getAllExpVar(shuffledActivity, nFolds);
                cvExpVarNull_thisRep(iCell, thisShuff) = mean(shuffledFolds, 'omitnan');
            end
        end

        meanExpVar_thisRep = mean(cvExpVar_thisRep, 2, 'omitnan');
        countBetterNulls   = sum(cvExpVarNull_thisRep >= meanExpVar_thisRep, 2, 'omitnan');
        pValues_thisRep    = (countBetterNulls + 1) / (nShuffles + 1);

        meanExpVar_reps(:, r) = meanExpVar_thisRep;
        pValues_reps(:, r)    = pValues_thisRep;
    end

    % average across repeats -> stable per-ROI estimate at matched N
    downsampledResults(s).meanExpVar = mean(meanExpVar_reps, 2, 'omitnan');
    downsampledResults(s).pValues    = mean(pValues_reps, 2, 'omitnan');

    fprintf('Mouse %s Day %d: downsampled to %d laps (%d repeats) done.\n', ...
        thisSession.MouseID, thisSession.Day, targetNumLaps, nRepeats);
end