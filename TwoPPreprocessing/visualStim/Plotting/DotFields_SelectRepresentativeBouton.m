% DotFields_SelectRepresentativeBouton.m
%
% Picks TWO representative boutons to illustrate the running-vs-
% stationary effect from both directions:
%   - a "running-preferring" example: dynamic range much bigger during
%     running than stationary
%   - a "stationary-preferring" example: dynamic range much bigger
%     during stationary than running (the less common, opposite pattern)
% Both are drawn from boutons with solid cross-validated tuning in BOTH
% states (dual R^2 filter), so neither pick is just a noisy fluke.

r2_thresh  = 0.1;
r2p_thresh = 0.05;

validIdx_stat = find(cat(1,allDotUnits.statR2) > r2_thresh & cat(1,allDotUnits.statR2_pval) < r2p_thresh);
validIdx_run  = find(cat(1,allDotUnits.runR2)  > r2_thresh & cat(1,allDotUnits.runR2_pval)  < r2p_thresh);
validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));

if isempty(validIdx_both)
    error('No boutons pass the dual R^2 filter -- cannot pick a clean representative example.');
end

dynRangeGap = [allDotUnits(validIdx_both).dynamicRange_run] - [allDotUnits(validIdx_both).dynamicRange_stat];
[~, sortIdx] = sort(dynRangeGap, 'descend');
rankedCandidates = validIdx_both(sortIdx); % most running-preferring first, most stationary-preferring last

fprintf('=== Top 5 RUNNING-preferring candidates (run dynamic range >> stat) ===\n');
for i = 1:min(5, numel(rankedCandidates))
    b = rankedCandidates(i);
    fprintf('  %d) Bouton index %d | %s | ROI %d | statR2=%.2f runR2=%.2f | dynRange stat=%.3f run=%.3f\n', ...
        i, b, allDotUnits(b).sessionLabel, allDotUnits(b).roiIdx, ...
        allDotUnits(b).statR2, allDotUnits(b).runR2, ...
        allDotUnits(b).dynamicRange_stat, allDotUnits(b).dynamicRange_run);
end

fprintf('\n=== Top 5 STATIONARY-preferring candidates (stat dynamic range >> run) ===\n');
for i = 1:min(5, numel(rankedCandidates))
    b = rankedCandidates(end - i + 1);
    fprintf('  %d) Bouton index %d | %s | ROI %d | statR2=%.2f runR2=%.2f | dynRange stat=%.3f run=%.3f\n', ...
        i, b, allDotUnits(b).sessionLabel, allDotUnits(b).roiIdx, ...
        allDotUnits(b).statR2, allDotUnits(b).runR2, ...
        allDotUnits(b).dynamicRange_stat, allDotUnits(b).dynamicRange_run);
end

runPreferringIdx  = rankedCandidates(1);
statPreferringIdx = rankedCandidates(end);

fprintf('\nSelected RUNNING-preferring:    bouton %d -- mouseID=%s, sessionName=%s, roiIdx=%d\n', ...
    runPreferringIdx, allDotUnits(runPreferringIdx).mouseID, allDotUnits(runPreferringIdx).sessionName, ...
    allDotUnits(runPreferringIdx).roiIdx);
fprintf('Selected STATIONARY-preferring: bouton %d -- mouseID=%s, sessionName=%s, roiIdx=%d\n', ...
    statPreferringIdx, allDotUnits(statPreferringIdx).mouseID, allDotUnits(statPreferringIdx).sessionName, ...
    allDotUnits(statPreferringIdx).roiIdx);

fprintf('\nCalls:\n');
fprintf('plotDotFieldsExampleBoutonStateSplit(''%s'', ''%s'', %d)  %% running-preferring\n', ...
    allDotUnits(runPreferringIdx).mouseID, allDotUnits(runPreferringIdx).sessionName, allDotUnits(runPreferringIdx).roiIdx);
fprintf('plotDotFieldsExampleBoutonStateSplit(''%s'', ''%s'', %d)  %% stationary-preferring\n', ...
    allDotUnits(statPreferringIdx).mouseID, allDotUnits(statPreferringIdx).sessionName, allDotUnits(statPreferringIdx).roiIdx);

%%
plotDotFieldsExampleBoutonStateSplit(allDotUnits(runPreferringIdx).mouseID, ...
    allDotUnits(runPreferringIdx).sessionName, allDotUnits(runPreferringIdx).roiIdx);

plotDotFieldsExampleBoutonStateSplit(allDotUnits(statPreferringIdx).mouseID, ...
    allDotUnits(statPreferringIdx).sessionName, allDotUnits(statPreferringIdx).roiIdx);