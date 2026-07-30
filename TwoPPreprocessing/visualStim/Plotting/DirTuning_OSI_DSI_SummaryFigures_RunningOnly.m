% DirTuning_OSI_DSI_SummaryFigures_RunningOnly.m
%
% Population-level OSI/DSI summary for RUNNING-ONLY boutons, gated on
% isTunedCVR2_run (cross-validated R^2, running trials only).
%
% Requires allDirTuning already has, per bouton:
%   isTunedCVR2_run, cvR2_run, minTrialUsed_run   (from the running-only
%       pooling + cvR2 script)
%   OSI_simple, prefOrientationDeg_simple          (from
%       computeDirTuningOSI(allDirTuning, 'isTunedCVR2_run'))
%   DSI_simple, prefDirectionDeg_simple             (from
%       computeDirTuningDSI(allDirTuning, 'isTunedCVR2_run'))
%
% If you haven't already run the OSI/DSI compute step on this struct,
% uncomment the two lines below.

% allDirTuning = computeDirTuningOSI(allDirTuning, 'isTunedCVR2_run', false);
% allDirTuning = computeDirTuningDSI(allDirTuning, 'isTunedCVR2_run', false);

selectivityThresh = 0.2; % common OSI/DSI cutoff for calling a bouton "selective" -- adjust to match your convention

selectedIdx = find([allDirTuning.isTunedCVR2_run]);
fprintf('%d / %d boutons pass isTunedCVR2_run.\n', numel(selectedIdx), numel(allDirTuning));

if isempty(selectedIdx)
    error('No boutons pass isTunedCVR2_run -- nothing to plot.');
end

OSI_sel  = [allDirTuning(selectedIdx).OSI_simple]';
DSI_sel  = [allDirTuning(selectedIdx).DSI_simple]';
prefOri_sel = [allDirTuning(selectedIdx).prefOrientationDeg_simple]';
prefDir_sel = [allDirTuning(selectedIdx).prefDirectionDeg_simple]';
cvR2_sel = [allDirTuning(selectedIdx).cvR2_run]';
minTrial_sel = [allDirTuning(selectedIdx).minTrialUsed_run]';

%% ===================== Figure 1: OSI/DSI population summary =====================
figure('Color', 'w', 'Position', [80 80 1300 750]);

subplot(2,3,1);
histogram(OSI_sel(~isnan(OSI_sel)), 20);
xline(selectivityThresh, 'r--', 'LineWidth', 1.2);
xlabel('OSI (simple ratio)'); ylabel('Number of boutons');
title(sprintf('OSI distribution (n=%d)', sum(~isnan(OSI_sel))));

subplot(2,3,2);
histogram(DSI_sel(~isnan(DSI_sel)), 20);
xline(selectivityThresh, 'r--', 'LineWidth', 1.2);
xlabel('DSI (simple ratio)'); ylabel('Number of boutons');
title(sprintf('DSI distribution (n=%d)', sum(~isnan(DSI_sel))));

subplot(2,3,3);
bothValid = ~isnan(OSI_sel) & ~isnan(DSI_sel);
scatter(OSI_sel(bothValid), DSI_sel(bothValid), 15, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
xline(selectivityThresh, 'k--'); yline(selectivityThresh, 'k--');
xlabel('OSI'); ylabel('DSI');
title(sprintf('OSI vs DSI relationship (n=%d)', sum(bothValid)));
axis square; xlim([0 1]); ylim([0 1]);

subplot(2,3,4);
polarhistogram(deg2rad(prefOri_sel(~isnan(prefOri_sel))) * 2, 16);
title('Preferred orientation (doubled for display)');

subplot(2,3,5);
polarhistogram(deg2rad(prefDir_sel(~isnan(prefDir_sel))), 16);
title('Preferred direction');

subplot(2,3,6);
isOS = OSI_sel > selectivityThresh;
isDS = DSI_sel > selectivityThresh;
nOSonly    = sum(isOS & ~isDS);
nDSonly    = sum(~isOS & isDS);
nBoth      = sum(isOS & isDS);
nNeither   = sum(~isOS & ~isDS);
bar([nOSonly, nDSonly, nBoth, nNeither]);
set(gca, 'XTickLabel', {'OS only', 'DS only', 'Both', 'Neither'});
ylabel('Number of boutons');
title(sprintf('OS/DS breakdown (thresh=%.2f)', selectivityThresh));

sgtitle(sprintf('OSI/DSI summary, RUNNING ONLY (isTunedCVR2\\_run, n=%d boutons)', numel(selectedIdx)), ...
    'Interpreter', 'tex');

%% ===================== Figure 2: diagnostic checks =====================
% Worth checking these BEFORE trusting the population numbers above:
% does OSI/DSI correlate with tuning reliability (cvR2) or with how many
% running trials a bouton had (minTrialUsed_run)? If selectivity tracks
% trial count rather than cvR2, that's a sign of a small-sample bias
% (noisier tuning curves can look artificially more "selective") rather
% than a real effect.

figure('Color', 'w', 'Position', [80 80 1000 700]);

subplot(2,2,1);
scatter(cvR2_sel, OSI_sel, 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Cross-validated R^2 (running only)'); ylabel('OSI');
title('OSI vs tuning reliability');

subplot(2,2,2);
scatter(cvR2_sel, DSI_sel, 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Cross-validated R^2 (running only)'); ylabel('DSI');
title('DSI vs tuning reliability');

subplot(2,2,3);
scatter(minTrial_sel, OSI_sel, 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Trials/direction used (after downsampling)'); ylabel('OSI');
title('OSI vs running-trial count');

subplot(2,2,4);
scatter(minTrial_sel, DSI_sel, 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Trials/direction used (after downsampling)'); ylabel('DSI');
title('DSI vs running-trial count');
% set(gca, 'XTick', 1:4, 'XTickLabel', {'OS only', 'DS only', 'Both', 'Neither'});
sgtitle('Diagnostics: does selectivity track reliability or trial count?', 'FontWeight', 'bold');

fprintf('\nSelectivity summary (thresh=%.2f): %d OS-only, %d DS-only, %d both, %d neither (n=%d)\n', ...
    selectivityThresh, nOSonly, nDSonly, nBoth, nNeither, numel(selectedIdx));
