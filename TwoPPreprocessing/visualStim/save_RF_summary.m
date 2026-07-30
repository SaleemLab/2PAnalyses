%% ================== SAVE RF SUMMARY ==================
% Paste this at the end of your existing pipeline, AFTER isResp / gaussFitResults /
% trustedMask / trustedIROI / sessionLabels have all been computed for the current
% dataset (V1 somas, V1 boutons, or RSP boutons). Just change datasetName each run.
%
% Key fix vs. raw trustedMask: gaussFitResults / trustedMask only cover the RESPONSIVE
% subset of ROIs (that's what got fit), so their length/order does NOT match
% isResponsive / sessionLabels (which cover ALL ROIs). We re-expand trusted status back
% out to full-ROI-length using iROI, so everything in the saved summary is indexed
% consistently 1:numBoutons and can be compared/plotted together later.

datasetName = 'RSP_boutons';   % <-- CHANGE THIS EACH RUN: 'V1_somas' / 'V1_boutons' / 'RSP_boutons'

isTrustedFull = false(1, numBoutons);
isTrustedFull(trustedIROI) = true;   % trustedIROI = actual ROI indices (into allRFMapping) that passed

summary = struct();
summary.datasetName    = datasetName;
summary.nTotal         = length(allMeanCrossValR2);                 % all ROIs recorded (responsive + not)
summary.isResponsive   = isResp;                      % logical, 1 x numBoutons
summary.isTrustedFull  = isTrustedFull;                % logical, 1 x numBoutons -- responsive AND good/trusted fit
summary.sessionLabels  = sessionLabels;                % cell, numBoutons x 1

summaryDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Supp_Section1_Fig4_1_VISp\summaries';
if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end
save(fullfile(summaryDir, sprintf('%s_summary.mat', datasetName)), 'summary');

fprintf('[%s] saved: %d total | %d responsive (%.1f%%) | %d trusted fit (%.1f%%)\n', ...
    datasetName, summary.nTotal, ...
    sum(summary.isResponsive), 100*sum(summary.isResponsive)/summary.nTotal, ...
    sum(summary.isTrustedFull), 100*sum(summary.isTrustedFull)/summary.nTotal);

%% ---- If you're NOT rerunning V1 and just have the numbers written down ----
% Skip everything above and build the summary struct manually instead, e.g.:
%
summary = struct();
summary.datasetName   = 'V1_somas';
summary.nTotal        = 214;                 % <-- your total n
summary.nResponsive   = 130;                 % <-- your responsive n
summary.nTrusted      = 98;                  % <-- your trusted/good-fit n
summary.sessionLabels = {};                  % leave empty -- no per-session scatter for this one
save(fullfile(summaryDir, 'V1_somas_summary.mat'), 'summary');
%
% The plotting script below handles either format (isResponsive vector OR nResponsive scalar).
