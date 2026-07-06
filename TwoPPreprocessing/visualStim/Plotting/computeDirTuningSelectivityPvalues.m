function allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning, responsivenessField, nShuffle)
% COMPUTEDIRTUNINGSELECTIVITYPVALUES  Permutation test for OSI/DSI.
%
%   allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning)
%   allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning, responsivenessField, nShuffle)
%
% WHY THIS EXISTS: OSI/DSI computed from meanDirResponse alone are
% unprotected point estimates. Two failure modes were found empirically
% today: (1) isTunedCV only gates on a p-value with no minimum effect
% size, so weakly/non-tuned boutons can still pass; (2) half-wave
% rectification (needed since a negative "vector length" is undefined)
% means that for a mostly-noise bouton, ~half its 8 direction means go
% negative by chance and get zeroed out, leaving only 1-2 small positive
% survivors -- the vector sum of just those can look artificially
% concentrated (near-1 OSI/DSI) despite there being no real tuning.
%
% This function directly tests whether a bouton's REAL OSI/DSI is bigger
% than what you'd get from PURE NOISE with the same trial counts, by
% shuffling which direction each trial's response is assigned to
% (destroying any true direction/orientation structure) many times, and
% comparing the real value against that null distribution.
%
% INPUTS:
%   allDirTuning        - struct array; must have trialMeanResp,
%                         stimValues per bouton.
%   responsivenessField - (optional, default 'isResponsive_ttest') only
%                         boutons passing this are tested (saves time).
%   nShuffle             - (optional, default 500) number of shuffles.
%
% OUTPUT: adds these fields per bouton:
%   .OSI_pval, .DSI_pval  - permutation p-values
%   .isOSI_reliable       - OSI_pval < 0.05
%   .isDSI_reliable       - DSI_pval < 0.05

if nargin < 2 || isempty(responsivenessField)
    responsivenessField = 'isResponsive_ttest';
end
if nargin < 3 || isempty(nShuffle)
    nShuffle = 500;
end

nBoutonsTotal = numel(allDirTuning);
OSI_pval = nan(nBoutonsTotal, 1);
DSI_pval = nan(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    s = allDirTuning(b);
    if ~isfield(s, responsivenessField) || ~s.(responsivenessField)
        continue;
    end
    if ~isfield(s, 'trialMeanResp') || ~isfield(s, 'stimValues')
        continue;
    end

    thetaDeg = s.stimValues(:)';
    nDir = numel(thetaDeg);

    % pool all trials with a direction label, same pattern as the ANOVA setup
    y = []; grp = [];
    for d = 1:nDir
        vals = s.trialMeanResp{d}(:);
        vals = vals(~isnan(vals));
        y   = [y; vals];
        grp = [grp; repmat(d, numel(vals), 1)];
    end
    if numel(y) < nDir * 2
        continue;
    end

    [realOSI, realDSI] = localComputeIndices(y, grp, thetaDeg, nDir);
    if isnan(realOSI) && isnan(realDSI)
        continue;
    end

    nullOSI = nan(nShuffle, 1);
    nullDSI = nan(nShuffle, 1);
    for sh = 1:nShuffle
        shuffledGrp = grp(randperm(numel(grp)));
        [nullOSI(sh), nullDSI(sh)] = localComputeIndices(y, shuffledGrp, thetaDeg, nDir);
    end

    OSI_pval(b) = (sum(nullOSI >= realOSI, 'omitnan') + 1) / (nShuffle + 1);
    DSI_pval(b) = (sum(nullDSI >= realDSI, 'omitnan') + 1) / (nShuffle + 1);

    if mod(b, 50) == 0
        fprintf('Permutation testing OSI/DSI: %d / %d boutons...\n', b, nBoutonsTotal);
    end
end

for b = 1:nBoutonsTotal
    allDirTuning(b).OSI_pval       = OSI_pval(b);
    allDirTuning(b).DSI_pval       = DSI_pval(b);
    allDirTuning(b).isOSI_reliable = OSI_pval(b) < 0.05;
    allDirTuning(b).isDSI_reliable = DSI_pval(b) < 0.05;
end

fprintf('\n%d boutons have a reliable (permutation-tested) OSI.\n', sum(OSI_pval < 0.05, 'omitnan'));
fprintf('%d boutons have a reliable (permutation-tested) DSI.\n', sum(DSI_pval < 0.05, 'omitnan'));

end

function [osi, dsi] = localComputeIndices(y, grp, thetaDeg, nDir)
    R = nan(1, nDir);
    for d = 1:nDir
        R(d) = mean(y(grp == d), 'omitnan');
    end
    R_rect = max(R, 0);

    if sum(R_rect) == 0 || all(isnan(R_rect))
        osi = NaN; dsi = NaN;
        return;
    end

    thetaRad = deg2rad(thetaDeg);
    vecOSI = sum(R_rect .* exp(1i * 2 * thetaRad));
    osi = abs(vecOSI) / sum(R_rect);

    vecDSI = sum(R_rect .* exp(1i * thetaRad));
    dsi = abs(vecDSI) / sum(R_rect);
end
