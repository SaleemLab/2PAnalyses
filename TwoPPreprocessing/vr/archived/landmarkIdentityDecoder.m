%% Landmark Identity Decoding - RSP vs VISp
% Train on odd laps, test on even laps
% Decode: is the animal at landmark type A (40,120cm) or type B (80,160cm)?

function decodingResults = landmarkIdentityDecoder(RegionData, regionName)

if nargin < 2, regionName = 'Region'; end

% landmark windows ±10 bins
lmA_bins = {35:45, 115:125};  % type A: 40cm and 120cm
lmB_bins = {75:85, 155:165};  % type B: 80cm and 160cm

allAcc = [];

for s = 1:length(RegionData)
    sess = RegionData(s);
    
    if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
       ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
        continue;
    end
    
    lapActivity = sess.ConditionData.Baseline.LapActivity;
    [numROIsTotal, numLaps, numBins] = size(lapActivity);
    if numBins < 170, continue; end
    
    % smooth
    w = gausswin(15); w = w/sum(w);
    smoothed = lapActivity;
    for iCell = 1:numROIsTotal
        for iLap = 1:numLaps
            tr = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(tr)), continue; end
            nm = isnan(tr); tr(nm) = 0;
            sm = filtfilt(w, 1, tr); sm(nm) = NaN;
            smoothed(iCell, iLap, :) = sm;
        end
    end
    
    roiActivity = smoothed(sess.FilteredROIs, :, :);  % nCells x nLaps x nBins
    nCells = length(sess.FilteredROIs);
    
    oddLaps  = 1:2:numLaps;
    evenLaps = 2:2:numLaps;
    
    % --- TRAINING: mean response at each landmark on odd laps ---
    % type A tuning: mean across A landmark windows
    trainA = [];
    trainB = [];
    for lm = 1:2
        % mean across bins in this landmark window, per lap
        actA = squeeze(mean(roiActivity(:, oddLaps, lmA_bins{lm}), 3, 'omitnan'));  % nCells x nOddLaps
        actB = squeeze(mean(roiActivity(:, oddLaps, lmB_bins{lm}), 3, 'omitnan'));  % nCells x nOddLaps
        trainA = [trainA, actA];  % nCells x (2*nOddLaps)
        trainB = [trainB, actB];
    end
    
    % mean tuning at each landmark type on odd laps
    meanA_train = mean(trainA, 2, 'omitnan');  % nCells x 1
    meanB_train = mean(trainB, 2, 'omitnan');  % nCells x 1
    
    % --- TESTING: decode landmark identity on even laps ---
    correct = 0; total = 0;
    
    for lm = 1:2
        % test activity at A landmarks on even laps
        testA = squeeze(mean(roiActivity(:, evenLaps, lmA_bins{lm}), 3, 'omitnan'));  % nCells x nEvenLaps
        testB = squeeze(mean(roiActivity(:, evenLaps, lmB_bins{lm}), 3, 'omitnan'));
        
        for iLap = 1:length(evenLaps)
            popVec = testA(:, iLap);
            if sum(~isnan(popVec)) < 3, continue; end
            
            % nearest centroid classifier
            distA = norm(popVec(~isnan(popVec)) - meanA_train(~isnan(popVec)));
            distB = norm(popVec(~isnan(popVec)) - meanB_train(~isnan(popVec)));
            
            if distA < distB
                correct = correct + 1;  % correctly decoded as A
            end
            total = total + 1;
        end
        
        for iLap = 1:length(evenLaps)
            popVec = testB(:, iLap);
            if sum(~isnan(popVec)) < 3, continue; end
            
            distA = norm(popVec(~isnan(popVec)) - meanA_train(~isnan(popVec)));
            distB = norm(popVec(~isnan(popVec)) - meanB_train(~isnan(popVec)));
            
            if distB < distA
                correct = correct + 1;  % correctly decoded as B
            end
            total = total + 1;
        end
    end
    
    sessionAcc = correct / total;
    allAcc = [allAcc; sessionAcc];
    fprintf('Session %d: accuracy = %.1f%%\n', s, sessionAcc*100);
end

decodingResults.accuracy = allAcc;
decodingResults.meanAcc  = mean(allAcc, 'omitnan');
decodingResults.semAcc   = std(allAcc, 'omitnan') / sqrt(numel(allAcc));

fprintf('\n=== %s Landmark Identity Decoding ===\n', regionName);
fprintf('Mean accuracy: %.1f%% ± %.1f%% SEM\n', decodingResults.meanAcc*100, decodingResults.semAcc*100);
fprintf('Chance level: 50%%\n');

end

