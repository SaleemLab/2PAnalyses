function allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning, responsivenessField, nShuffle)
% COMPUTEDIRTUNINGSELECTIVITYPVALUES  Permutation test for OSI/DSI.
%
%   allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning)
%   allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning, responsivenessField, nShuffle)
%
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
OSIsimple_pval = nan(nBoutonsTotal, 1);
DSIsimple_pval = nan(nBoutonsTotal, 1);

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

    [realOSI, realDSI, realOSIsimple, realDSIsimple] = localComputeIndices(y, grp, thetaDeg, nDir);
    if isnan(realOSI) && isnan(realDSI) && isnan(realOSIsimple) && isnan(realDSIsimple)
        continue;
    end

    nullOSI = nan(nShuffle, 1);
    nullDSI = nan(nShuffle, 1);
    nullOSIsimple = nan(nShuffle, 1);
    nullDSIsimple = nan(nShuffle, 1);
    for sh = 1:nShuffle
        shuffledGrp = grp(randperm(numel(grp)));
        [nullOSI(sh), nullDSI(sh), nullOSIsimple(sh), nullDSIsimple(sh)] = localComputeIndices(y, shuffledGrp, thetaDeg, nDir);
    end

    OSI_pval(b) = (sum(nullOSI >= realOSI, 'omitnan') + 1) / (nShuffle + 1);
    DSI_pval(b) = (sum(nullDSI >= realDSI, 'omitnan') + 1) / (nShuffle + 1);
    OSIsimple_pval(b) = (sum(nullOSIsimple >= realOSIsimple, 'omitnan') + 1) / (nShuffle + 1);
    DSIsimple_pval(b) = (sum(nullDSIsimple >= realDSIsimple, 'omitnan') + 1) / (nShuffle + 1);

    if mod(b, 50) == 0
        fprintf('Permutation testing OSI/DSI: %d / %d boutons...\n', b, nBoutonsTotal);
    end
end

for b = 1:nBoutonsTotal
    allDirTuning(b).OSI_pval             = OSI_pval(b);
    allDirTuning(b).DSI_pval             = DSI_pval(b);
    allDirTuning(b).OSIsimple_pval       = OSIsimple_pval(b);
    allDirTuning(b).DSIsimple_pval       = DSIsimple_pval(b);
    allDirTuning(b).isOSI_reliable       = OSI_pval(b) < 0.05;
    allDirTuning(b).isDSI_reliable       = DSI_pval(b) < 0.05;
    allDirTuning(b).isOSIsimple_reliable = OSIsimple_pval(b) < 0.05;
    allDirTuning(b).isDSIsimple_reliable = DSIsimple_pval(b) < 0.05;
end

fprintf('\n%d boutons have a reliable (permutation-tested) vector OSI.\n', sum(OSI_pval < 0.05, 'omitnan'));
fprintf('%d boutons have a reliable (permutation-tested) vector DSI.\n', sum(DSI_pval < 0.05, 'omitnan'));
fprintf('%d boutons have a reliable (permutation-tested) simple-ratio OSI.\n', sum(OSIsimple_pval < 0.05, 'omitnan'));
fprintf('%d boutons have a reliable (permutation-tested) simple-ratio DSI.\n', sum(DSIsimple_pval < 0.05, 'omitnan'));

end

function [osi, dsi, osiSimple, dsiSimple] = localComputeIndices(y, grp, thetaDeg, nDir)
    R = nan(1, nDir);
    for d = 1:nDir
        R(d) = mean(y(grp == d), 'omitnan');
    end
    R_rect = max(R, 0);

    if sum(R_rect) == 0 || all(isnan(R_rect))
        osi = NaN; dsi = NaN; osiSimple = NaN; dsiSimple = NaN;
        return;
    end

    thetaRad = deg2rad(thetaDeg);
    vecOSI = sum(R_rect .* exp(1i * 2 * thetaRad));
    osi = abs(vecOSI) / sum(R_rect);

    vecDSI = sum(R_rect .* exp(1i * thetaRad));
    dsi = abs(vecDSI) / sum(R_rect);

    % --- simple ratio OSI: fold opposite-direction pairs into 4 orientations ---
    osiSimple = NaN;
    thetaDegMod = mod(round(thetaDeg), 360);
    orientBins = mod(thetaDegMod, 180);
    uniqueOrients = unique(orientBins);
    if numel(uniqueOrients) == 4
        orientResponse = nan(1, 4);
        for oi = 1:4
            orientResponse(oi) = mean(R_rect(orientBins == uniqueOrients(oi)), 'omitnan');
        end
        [Rpref, prefIdx] = max(orientResponse);
        orthoIdx = mod(prefIdx - 1 + 2, 4) + 1;
        Rortho = orientResponse(orthoIdx);
        if (Rpref + Rortho) > 0
            osiSimple = (Rpref - Rortho) / (Rpref + Rortho);
        end
    end

    % simple ratio DSI: preferred direction vs its exact opposite (180 deg away) ---
    dsiSimple = NaN;
    [thetaSorted, sortIdx] = sort(thetaDegMod);
    Rsorted = R_rect(sortIdx);
    spacing = mode(diff(thetaSorted));
    stepsFor180 = round(180 / spacing);
    if abs(stepsFor180 * spacing - 180) <= 1e-6 && mod(nDir, 2) == 0
        [Rpref, prefIdx] = max(Rsorted);
        oppIdx = mod(prefIdx - 1 + stepsFor180, nDir) + 1;
        Ropp = Rsorted(oppIdx);
        if (Rpref + Ropp) > 0
            dsiSimple = (Rpref - Ropp) / (Rpref + Ropp);
        end
    end
end
