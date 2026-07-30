% Manipulation corridor 
pairs=struct;
pairs.M25132 = ['20260313', '20260228']; 
pairs.M26003 = ['20260324', '20260325'];

ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

RSPDataSpks = getTuningDataByCondition(ExpRSPSessions, 'signalToUse', 'spks');
% Filter rois based on critera defined; here i include all rois or only the
% ones that do not have peak responses in the edges?

RSPDataSpks = appendFilteredROIs(RSPDataSpks,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true, 'FilterDuplicateBoutons', true);


plotThreeConditions_DifferenceIncluded(RSPDataSpks, 'RSP', ...
    'TypeToPlot', 'Boutons', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.8\poplation_andDifference\conditions_pop_with_diff_dff')



 [incIdx_o2, decIdx_o2, deltas_o2, incTable_o2, decTable_o2, ...
          incIdx_o3, decIdx_o3, deltas_o3, incTable_o3, decTable_o3, ...
          figHandle] = findLandmarkModulatedROIs(RSPDataSpks, 'Baseline', 'Omit_2', 'Omit_3');


 [suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPDataSpks, 'RSP', ...
    'zThresh', -0.5, 'windowCm', 15);

% Return
sessionIDs = allSessionIDs;
roiIDs     = allROIIDs;

%%
