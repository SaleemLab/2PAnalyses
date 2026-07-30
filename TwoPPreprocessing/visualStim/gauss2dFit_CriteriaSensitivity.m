% RFMapping_CriteriaSensitivity.m
%
% Breaks down, gate by gate, how many currently-excluded fits would
% become "trusted" if that ONE gate were dropped (all others held fixed)
% -- i.e. the marginal contribution of each individual criterion to the
% final trusted-set size. Also reports how many fits each gate flags in
% total (regardless of overlap with other failures), since a gate that
% rarely fires ALONE might still flag a lot of fits that also fail
% elsewhere -- those wouldn't be recovered by dropping that gate alone.
%
% Six gates total: {BeyondRange, DegenerateSigma, CoverageFail} x {Az, El}.
% "isTrusted" in the main script = TRUE only if ALL SIX pass.
%
% REQUIRES: gaussFitResults (from RFMapping_Gaussian2DFit.m -- run that
% first and keep it in the workspace).

if ~exist('gaussFitResults', 'var')
    error('gaussFitResults not found -- run RFMapping_Gaussian2DFit.m first and keep it in the workspace.');
end
if ~exist('coverageThreshold', 'var')
    coverageThreshold = 0.95; % matches default in RFMapping_Gaussian2DFit.m -- override if you used a different value
end

nFits = numel(gaussFitResults);

% Pull per-bouton values
beyondAz = [gaussFitResults.isBeyondRangeAz];
beyondEl = [gaussFitResults.isBeyondRangeEl];
degenAz  = [gaussFitResults.isDegenerateSigmaAz];
degenEl  = [gaussFitResults.isDegenerateSigmaEl];
pctAz    = [gaussFitResults.pctWithinStepAz];
pctEl    = [gaussFitResults.pctWithinStepEl];

% NaN coverage (bootstrap couldn't gather >=20 valid resamples) counts as a FAIL, matching how
% isCoverageOK behaves in the main script (comparisons against NaN are false in MATLAB, so this
% is automatic, but made explicit here for clarity).
covFailAz = ~(pctAz >= coverageThreshold);
covFailEl = ~(pctEl >= coverageThreshold);

% 6 individual gates
gateNames = {'BeyondRange-Az', 'BeyondRange-El', 'DegenerateSigma-Az', 'DegenerateSigma-El', ...
             'CoverageFail-Az', 'CoverageFail-El'};
F = [beyondAz(:), beyondEl(:), degenAz(:), degenEl(:), covFailAz(:), covFailEl(:)];
nGates = size(F, 2);

currentlyTrusted = ~any(F, 2);
storedTrusted = [gaussFitResults.isTrusted];
if ~isequal(currentlyTrusted(:), storedTrusted(:))
    warning(['Recomputed trusted set does not match gaussFit_isTrusted stored in gaussFitResults -- ' ...
        'check coverageThreshold above matches the value used in the main script run.']);
end

fprintf('\n=== Criteria sensitivity (n = %d fitted ROIs) ===\n', nFits);
fprintf('Currently trusted (all 6 gates pass): %d / %d (%.1f%%)\n\n', ...
    sum(currentlyTrusted), nFits, 100*sum(currentlyTrusted)/nFits);

fprintf('%-22s %10s %14s\n', 'Gate', 'FlagsAny', 'GainIfDropped');
totalFlagged = nan(1, nGates);
soleFail     = nan(1, nGates);
for g = 1:nGates
    otherGates  = setdiff(1:nGates, g);
    failsOthers = any(F(:, otherGates), 2);
    totalFlagged(g) = sum(F(:, g));
    soleFail(g)     = sum(F(:, g) & ~failsOthers); % fails ONLY this gate -> would flip to trusted if dropped
    fprintf('%-22s %10d %14d\n', gateNames{g}, totalFlagged(g), soleFail(g));
end

fprintf('\n"FlagsAny"       = how many fits this gate excludes, counting overlap with other failed gates.\n');
fprintf('"GainIfDropped"  = how many fits fail ONLY this gate -- i.e. how many would move from\n');
fprintf('                   excluded to trusted if this ONE gate were dropped, holding all others fixed.\n');
fprintf('                   This is the true marginal value of each gate: a fit failing multiple gates\n');
fprintf('                   stays excluded regardless of any single one being dropped.\n');

%% pairwise overlap of failures (which gates tend to co-occur)
fprintf('\n=== Pairwise co-occurrence of gate failures (n fits failing BOTH gates) ===\n');
fprintf('%-24s', '');
for g = 1:nGates
    fprintf('%16s', gateNames{g});
end
fprintf('\n');
for g1 = 1:nGates
    fprintf('%-24s', gateNames{g1});
    for g2 = 1:nGates
        fprintf('%16d', sum(F(:,g1) & F(:,g2)));
    end
    fprintf('\n');
end
fprintf('(diagonal = total flagged by that gate alone, regardless of overlap; matches "FlagsAny" above)\n');

%% bar chart of marginal gain per gate
figure('Color', 'w', 'Name', 'Criteria sensitivity');
subplot(1,2,1);
bar(totalFlagged);
set(gca, 'XTick', 1:nGates, 'XTickLabel', gateNames, 'XTickLabelRotation', 30);
ylabel('# fits flagged (any overlap)');
title('Total flagged per gate');

subplot(1,2,2);
bar(soleFail);
set(gca, 'XTick', 1:nGates, 'XTickLabel', gateNames, 'XTickLabelRotation', 30);
ylabel('# fits gained if dropped alone');
title('Marginal gain if gate dropped');
sgtitle(sprintf('Currently trusted: %d / %d', sum(currentlyTrusted), nFits));

%% scenario: drop an entire criterion TYPE (both axes) at once -- more realistic than
%% dropping just one axis of one criterion, since in practice you'd likely loosen/tighten
%% a criterion symmetrically across axes rather than asymmetrically
fprintf('\n=== Scenario: drop an entire criterion TYPE (both axes) at once ===\n');
scenarioNames      = {'Drop BeyondRange (both axes)', 'Drop DegenerateSigma (both axes)', 'Drop Coverage (both axes)'};
scenarioGateGroups = {[1 2], [3 4], [5 6]};
for s = 1:numel(scenarioGateGroups)
    keepGates      = setdiff(1:nGates, scenarioGateGroups{s});
    wouldBeTrusted = ~any(F(:, keepGates), 2);
    gained         = sum(wouldBeTrusted) - sum(currentlyTrusted);
    fprintf('%-32s : %d / %d trusted (currently %d) -- gain of %d\n', ...
        scenarioNames{s}, sum(wouldBeTrusted), nFits, sum(currentlyTrusted), gained);
end

%% cross-check against isResponsive: of the fits that WOULD be gained under each scenario,
%% how many are actually visually responsive vs. just noise that cleared the R^2 floor?
%% (relevant since the main script now fits ALL ROIs, not just responsive ones)
isResp = [gaussFitResults.isResponsive];
fprintf('\n=== Of currently trusted, how many are isResponsive? ===\n');
fprintf('%d / %d currently trusted fits are isResponsive (%.1f%%).\n', ...
    sum(currentlyTrusted & isResp(:)), sum(currentlyTrusted), 100*sum(currentlyTrusted & isResp(:))/max(sum(currentlyTrusted),1));
fprintf('If this is well below 100%%, your trusted set includes non-responsive boutons whose noise\n');
fprintf('happened to clear every QC gate -- worth ANDing with isResponsive downstream regardless of\n');
fprintf('which gates you keep.\n');
