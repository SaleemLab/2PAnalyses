function reliabilityIdx = computeVisualReliabilityIndex_DotFields_fromScalars(thisBouton, allDotUnits, si)
% Stimulus-independent-style reliability metric using ONLY the scalar
% per-trial values already stored in allDotUnits (no script changes needed).
% si = 1 (stationary) or 2 (locomotion)

alltraces   = allDotUnits(thisBouton).alltraces(:, si);   % cell, nSpeeds x 1
blankTrials = allDotUnits(thisBouton).blankTrials{si};    % vector

allConditions = alltraces(:);
allConditions = allConditions(~cellfun(@isempty, allConditions));
if ~isempty(blankTrials)
    allConditions{end+1} = blankTrials(:)';
end

% Clean NaNs within each condition
allConditions = cellfun(@(x) x(~isnan(x)), allConditions, 'UniformOutput', false);

nTrialsPerCond = cellfun(@numel, allConditions);
minTrials = min(nTrialsPerCond);
if isempty(minTrials) || minTrials < 4
    reliabilityIdx = NaN;
    return;
end

% Build [minTrials x nConditions] matrix: row = one repeat's profile
% across all conditions
profileMatrix = nan(minTrials, numel(allConditions));
for c = 1:numel(allConditions)
    vals = allConditions{c};
    profileMatrix(:, c) = vals(1:minTrials);
end

R = corr(profileMatrix', 'rows', 'pairwise');   % correlate REPEATS (rows), not conditions
offDiag = R(~eye(size(R)));
reliabilityIdx = prctile(offDiag, 75);
end