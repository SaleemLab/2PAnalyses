% RFMapping_SelectSplitOnlyExamples.m
%
% Picks a few example boutons from the "split-only" category (responsive
% in stationary OR running, but NOT when all trials are combined) for
% visual spot-checking -- to see whether these look like genuine,
% spatially coherent RFs or noise-driven artifacts.

isComb   = [allRFBehaviorSplit.isResponsive_combined];
isStat   = [allRFBehaviorSplit.isResponsive_stat];
isRun    = [allRFBehaviorSplit.isResponsive_run];
isEither = isStat | isRun;

splitOnlyIdx = find(~isComb & isEither);
fprintf('%d boutons in the split-only category.\n', numel(splitOnlyIdx));

nExamples = min(6, numel(splitOnlyIdx));
exampleIdx = splitOnlyIdx(1:nExamples);

fprintf('\nExample split-only boutons to plot:\n');
for i = 1:nExamples
    b = exampleIdx(i);
    whichState = '';
    if isStat(b), whichState = [whichState 'STAT ']; end
    if isRun(b),  whichState = [whichState 'RUN'];  end
    fprintf('  %d) %s | ROI %d | responsive during: %s\n', ...
        i, allRFBehaviorSplit(b).sessionLabel, allRFBehaviorSplit(b).roiIdx, whichState);
    fprintf('     plotRFExampleBoutonStateSplit(''%s'', ''%s'', %d)\n', ...
        allRFBehaviorSplit(b).mouseID, allRFBehaviorSplit(b).sessionName, allRFBehaviorSplit(b).roiIdx);
end

%% plot them all
for i = 1:nExamples
    b = exampleIdx(i);
    plotRFExampleBoutonStateSplit(allRFBehaviorSplit(b).mouseID, allRFBehaviorSplit(b).sessionName, allRFBehaviorSplit(b).roiIdx);
end
