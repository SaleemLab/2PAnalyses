% pool 'experinenced' rsp sessions across sessions (where fovs are not
% matched); day 5 onwards; pool baseline trials accross manipulation runs 
%% Grab session to load with unique FOVs [5; 200] 
pairs=struct;
pairs.M25132 = ['20260226', '20260228', '20260313'];
pairs.M25133 = '20260224';
pairs.M26003 = ['20260322', '20260324', '20260325'];

ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
% bin data by condition and load inclusion critera
RSPData = getTuningDataByCondition(ExpRSPSessions);
% Filter rois based on critera defined 
RSPData = appendFilteredROIs(RSPData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1);
%Total Raw ROIs Processed: 2777; Total ROIs Retained: 669; Overall Exclusion Rate:   75.91%



%% plotting functions 
% plot odd even and save in pdf png and svg 
plotPooledPopulation_OddEven(RSPData, 'RSP', ...
    'DaysToPlot',[5 200], ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section1\Fig2.2_ExpSetup');


% plot even and conditions in pdf png and svg 
plotPooledPopulation_AcrossConditions(RSPData, 'RSP', ...
    'SavePath', '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1\RSPBoutons_ExpMice_5_200\RSPpooled_Conditions_ExpMice_OddEvenStable');

plotAndSaveFilteredTuningCurves(RSPData, true)

% load smi across sessions 
plotFilteredSMIDistributions(RSPData)


% this includes the difference plots
plotThreeConditions_DifferenceIncluded(RSPData, 'RSP', ...
    'TypeToPlot', 'Boutons', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\conditionsAndDiff_pooledAcrossMice_3Conditions_withoutHalves')
%%