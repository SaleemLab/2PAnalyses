function reliabilityIdx = computeVisualReliabilityIndex(iROI, allRFMapping)
% 75th percentile of pairwise trial correlations, computed across ALL
% conditions concatenated (grid + blank) -- stimulus-independent version.

trialMatrix      = allRFMapping(iROI).baselineSubtracted;
bTrialsCorrected  = allRFMapping(iROI).baselineSubtractedBlank;

if isempty(trialMatrix) || isempty(bTrialsCorrected)
    reliabilityIdx = NaN;
    return;
end

allConditions = trialMatrix(:);
allConditions = allConditions(~cellfun(@isempty, allConditions));
allConditions{end+1} = double(bTrialsCorrected);

% Find the minimum trial count across conditions so we can align trial
% INDEX (not identity) across conditions -- i.e. trial "repeat 1" from
% every condition gets concatenated together, forming one long trace
% per repeat number. This mirrors "de-randomized" in the quoted method.
nTrialsPerCond = cellfun(@(x) size(x,1), allConditions);
minTrials = min(nTrialsPerCond);
if minTrials < 4
    reliabilityIdx = NaN;
    return;
end

% Build [minTrials x totalTimepoints] matrix: each row = one repeat,
% concatenated across all conditions in a fixed order
concatTraces = [];
for c = 1:numel(allConditions)
    trials = allConditions{c};
    concatTraces = [concatTraces, trials(1:minTrials, :)]; %#ok<AGROW>
end

% Pairwise Pearson correlation across all trial-repeats
R = corr(concatTraces', 'rows', 'pairwise');  % [minTrials x minTrials]
offDiag = R(~eye(size(R)));
reliabilityIdx = prctile(offDiag, 75);
end