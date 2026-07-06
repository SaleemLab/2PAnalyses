% DirTuning_SummaryFigures.m
%
% Produces final summary figures for OSI_simple and DSI_simple (the
% SIMPLE RATIO method: (Rpref-Rortho)/(Rpref+Rortho) for OSI,
% (Rpref-Ropp)/(Rpref+Ropp) for DSI -- NOT the vector method):
%   1. OSI_simple distribution, with permutation-test-reliable boutons
%      highlighted separately from unreliable ones.
%   2. DSI_simple distribution, same treatment.
%   3. Polar scatter: each RELIABLE bouton plotted at its preferred
%      orientation angle, with radius = OSI_simple (angle doubled for
%      display, since orientation is periodic over 180 degrees).
%   4. Polar scatter: same idea for DSI_simple / preferred direction
%      (NOT doubled -- direction spans the full 360 degrees).
%
% "Reliable" = passed the permutation test in
% computeDirTuningSelectivityPvalues.m (OSIsimple_pval/DSIsimple_pval <
% 0.05) -- a direct test of whether the observed OSI_simple/DSI_simple
% exceeds what pure noise would produce with the same trial counts.
%
% Requires: allDirTuning must already have OSI_simple, DSI_simple
% (computeDirTuningOSI / computeDirTuningDSI), preferred angle fields,
% AND OSIsimple_pval/DSIsimple_pval/isOSIsimple_reliable/
% isDSIsimple_reliable (computeDirTuningSelectivityPvalues, updated
% version that also tests the simple ratio metrics).
%
% NOTE ON PREFERRED ANGLE: the simple-ratio method doesn't save its own
% "preferred angle" field directly (only vector method does, via
% prefOrientationDeg/prefDirectionDeg). This script recomputes the
% simple-ratio preferred orientation/direction directly from
% meanDirResponse for the polar scatter plots, so the angle used matches
% the OSI_simple/DSI_simple value being plotted.

if ~isfield(allDirTuning, 'isOSIsimple_reliable') || ~isfield(allDirTuning, 'isDSIsimple_reliable')
    error(['allDirTuning is missing isOSIsimple_reliable/isDSIsimple_reliable -- run the updated ' ...
        'computeDirTuningSelectivityPvalues(allDirTuning) (with simple-ratio testing) first.']);
end

nBoutonsTotal = numel(allDirTuning);
allOSIsimple = [allDirTuning.OSI_simple];
allDSIsimple = [allDirTuning.DSI_simple];
isOSIrel     = logical([allDirTuning.isOSIsimple_reliable]);
isDSIrel     = logical([allDirTuning.isDSIsimple_reliable]);

validOSI = ~isnan(allOSIsimple);
validDSI = ~isnan(allDSIsimple);

fprintf('OSI_simple: %d valid, %d reliable (permutation p<0.05)\n', sum(validOSI), sum(isOSIrel & validOSI));
fprintf('DSI_simple: %d valid, %d reliable (permutation p<0.05)\n', sum(validDSI), sum(isDSIrel & validDSI));

%% ===================== recompute simple-ratio preferred angles (for polar plots only) =====================
prefOrientSimple = nan(1, nBoutonsTotal);
prefDirSimple    = nan(1, nBoutonsTotal);

for b = 1:nBoutonsTotal
    s = allDirTuning(b);
    if ~isfield(s, 'stimValues') || ~isfield(s, 'meanDirResponse')
        continue;
    end
    thetaDeg = mod(round(s.stimValues(:)'), 360);
    R = s.meanDirResponse(:)';
    R_rect = max(R, 0);
    nDir = numel(thetaDeg);

    % preferred ORIENTATION (folded, for OSI_simple polar plot)
    orientBins = mod(thetaDeg, 180);
    uniqueOrients = unique(orientBins);
    if numel(uniqueOrients) == 4
        orientResponse = nan(1, 4);
        for oi = 1:4
            orientResponse(oi) = mean(R_rect(orientBins == uniqueOrients(oi)), 'omitnan');
        end
        [~, prefIdx] = max(orientResponse);
        prefOrientSimple(b) = uniqueOrients(prefIdx);
    end

    % preferred DIRECTION (raw, for DSI_simple polar plot)
    [~, prefIdxDir] = max(R_rect);
    prefDirSimple(b) = thetaDeg(prefIdxDir);
end

%%  OSI_simple distribution 
figure('Position', [100 100 500 400]);
edges = -1:0.1:1; % simple ratio can go slightly negative if pref happens to be lower than ortho after rectification quirks
histogram(allOSIsimple(validOSI & ~isOSIrel), edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6, 'DisplayName', 'Not reliable (shuffle p\geq0.05)');
hold on;
histogram(allOSIsimple(validOSI & isOSIrel), edges, 'FaceColor', [0.85 0.33 0.10], 'FaceAlpha', 0.8, 'DisplayName', 'Reliable (shuffle p<0.05)');
xlabel('OSI (simple ratio)'); ylabel('Number of boutons');
title(sprintf('OSI\\_simple distribution (n=%d valid, %d reliable)', sum(validOSI), sum(isOSIrel & validOSI)));
legend('Location', 'northwest');

%%  DSI_simple distribution 
figure('Position', [100 100 500 400]);
histogram(allDSIsimple(validDSI & ~isDSIrel), edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6, 'DisplayName', 'Not reliable (shuffle p\geq0.05)');
hold on;
histogram(allDSIsimple(validDSI & isDSIrel), edges, 'FaceColor', [0 0.45 0.74], 'FaceAlpha', 0.8, 'DisplayName', 'Reliable (shuffle p<0.05)');
xlabel('DSI (simple ratio)'); ylabel('Number of boutons');
title(sprintf('DSI\\_simple distribution (n=%d valid, %d reliable)', sum(validDSI), sum(isDSIrel & validDSI)));
legend('Location', 'northwest');

%% preferred orientation polar histogram 
% Polar, not linear, since orientation is circular data -- a linear
% histogram can visually split a real cluster that straddles the
% wrap-around point (e.g. boutons preferring ~350-10 degrees would show
% as two separate small bumps at opposite edges instead of one
% continuous cluster). Angle doubled for display since orientation is
% periodic over 180 degrees, not 360.
relOSI = validOSI & isOSIrel & ~isnan(prefOrientSimple);

figure('Position', [100 100 500 500]);
polarhistogram(deg2rad(prefOrientSimple(relOSI)) * 2, 16);
title(sprintf('Preferred orientation distribution (reliable boutons, n=%d)', sum(relOSI)));

%% preferred direction polar histogram 
relDSI = validDSI & isDSIrel & ~isnan(prefDirSimple);

figure('Position', [100 100 500 500]);
polarhistogram(deg2rad(prefDirSimple(relDSI)), 16); % full 360 range, no doubling needed for direction
title(sprintf('Preferred direction distribution (reliable boutons, n=%d)', sum(relDSI)));
