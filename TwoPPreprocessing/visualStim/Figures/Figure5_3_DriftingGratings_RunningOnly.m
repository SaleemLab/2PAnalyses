
% run this to load all variables 
%DirTuning_IncTuningCurveAnalysis_RunningOnly_2PData
DirTuning_IncTuningCurveAnalysis_PoolStates_2PData
allDirTuning = computeDirTuningOSI(allDirTuning, 'isTunedCVR2_comb');
allDirTuning = computeDirTuningDSI(allDirTuning, 'isTunedCVR2_comb');


%%
% allDirTuning = computeDirTuningOSI(allDirTuning);
% allDirTuning = computeDirTuningDSI(allDirTuning);

% plotExampleDirTuningBouton(allDirTuning, 1871)
% plotExampleDirTuningBouton(allDirTuning, 789)
% plotExampleDirTuningBouton(allDirTuning, 1838)

%DirTuning_PlotAllSelectedBoutons

% is it necessary to shuffle?
% allDirTuning = computeDirTuningSelectivityPvalues(allDirTuning); % now also tests OSI_simple/DSI_simple significance
% DirTuning_SummaryFigures


DirTuning_PlotAllSelectedBoutons_PooledStates

plotExampleDirTuningBouton(allDirTuning, 372, [], 30)


plotExampleDirTuningBouton(allDirTuning, 615, [], 30)


plotExampleDirTuningBouton(allDirTuning, 501, [], 10)