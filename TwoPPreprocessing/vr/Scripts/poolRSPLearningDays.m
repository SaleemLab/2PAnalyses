%poolRSPLearningDays
% days 1 to 5 plus day 200 included for m25132 (3 sessions) and m26003 (2 sessions)
% Grab sessions to load 
pairs=struct;
pairs.M25132 = ['20260219','20260220','20260221','20260223', '20260226', '20260228', '20260313'];
pairs.M25133 = ['20260219','20260220','20260221','20260223','20260224'];
pairs.M26003 = ['20260316','20260317', '20260320','20260321','20260322', '20260324', '20260325'];



RSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

%% Get spatial tuning data 
RSPDataAcrossDays = getTuningDataByCondition(RSPSessions);
% Filter rois based on critera defined 
RSPDataAcrossDays = appendFilteredROIs(RSPDataAcrossDays, 'RhoHalvesThreshold', 0.6, 'CvExpVarThreshold', 0.1, 'FilterDuplicateBoutons', true); %Total Raw ROIs Processed: 2777; Total ROIs Retained: 669; Overall Exclusion Rate:   75.91%

%% plots 
% finds sort index from odd and plots even 
plotPooledPopulation_DayWise_Even(RSPDataAcrossDays, ...
    'RSP', ...
    'DaysToPlot',[1 2 3 4 5 200], ...
    'SavePath', '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter3\RSPBoutons_Learning\RSPpooled_Learning');


RSPDataAcrossDays = getTuningDataByCondition(RSPSessions);
RSPDataAcrossDays = computeSMIAcrossSessions(RSPDataAcrossDays);
pooledRSP = poolSMIAcrossSessions(RSPDataAcrossDays);
plotPooledSMI_CDF(pooledRSP, 'RSP', 'Days');

%%
