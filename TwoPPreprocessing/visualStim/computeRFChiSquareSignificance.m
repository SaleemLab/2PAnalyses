function [chi2stat, pValue] = computeRFChiSquareSignificance(stimMatrix, trialResponse, nShuffles)
% computeRFChiSquareSignificance
%   Chi-square test for a significant receptive field, following the
%   shuffle-based approach: chi2 = sum_i sum_j (Ei-Oi,j)^2 / Ei
%
% Inputs:
%   stimMatrix    - [gridY x gridX x nTrials], binary/contrast stimulus,
%                    e.g. 1 where a grid location was "on" that trial
%   trialResponse - [nTrials x 1], one response value per trial
%                    (e.g. average dF/F or deconvolved activity in the
%                    response window for that trial)
%   nShuffles     - number of shuffles for the null distribution (default 1000)
%
% Outputs:
%   chi2stat - observed chi-square statistic
%   pValue   - proportion of shuffled chi2 stats >= observed

if nargin < 3, nShuffles = 1000; end

[gridY, gridX, nTrials] = size(stimMatrix);
gridLinIdx = nan(nTrials, 1);

% Identify which single grid location was active on each trial
for t = 1:nTrials
    activeIdx = find(stimMatrix(:,:,t) ~= 0);   % non-zero = stimulated location
    if isscalar(activeIdx)
        gridLinIdx(t) = activeIdx;
    elseif numel(activeIdx) > 1
        % more than one location active this trial - take the first,
        % or adapt this depending on your stimulus design
        gridLinIdx(t) = activeIdx(1);
    end
end

validTrials = ~isnan(gridLinIdx) & ~isnan(trialResponse);
gridLinIdx = gridLinIdx(validTrials);
resp = trialResponse(validTrials);

chi2stat = localChiSquare(gridLinIdx, resp, gridY*gridX);

% --- Build null distribution by shuffling grid location labels ---
nullStats = nan(nShuffles, 1);
for s = 1:nShuffles
    shuffledIdx = gridLinIdx(randperm(numel(gridLinIdx)));
    nullStats(s) = localChiSquare(shuffledIdx, resp, gridY*gridX);
end

pValue = mean(nullStats >= chi2stat);

end

% ------------------------------------------------------------------
function chi2 = localChiSquare(locIdx, resp, nLocations)
    grandMean = mean(resp);
    chi2 = 0;
    for loc = 1:nLocations
        trialsAtLoc = (locIdx == loc);
        m = sum(trialsAtLoc);
        if m == 0
            continue
        end
        Oi = mean(resp(trialsAtLoc));
        Ei = grandMean;
        chi2 = chi2 + m * ((Ei - Oi)^2) / Ei;
    end
end