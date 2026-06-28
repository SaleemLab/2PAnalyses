% Questions: should i swap between spikes and dff for rsp and visp 
pairs=struct;
pairs.M25132 = ['20260226', '20260228', '20260313'];
pairs.M25133 = '20260224';
pairs.M26003 = ['20260322', '20260324', '20260325'];

ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

RSPData = getTuningDataByCondition(ExpRSPSessions);
% Filter rois based on critera defined; this alos include rois that peak;
% this exlcudes approx 300 rois that peak in the beginning or end of the
% average 
RSPData = appendFilteredROIs(RSPData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);
% this will be a slightly smaller pool; to extract rois that peak in the
% start and end of the lap; 

%% apply the same filter to visp data; include all sessions with baseline trials; 
% Load visp data as well for the cgfs
pairs=struct; 
pairs.M26005 = ['20260305', '20260306', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M25126 = ['20260311', '20260312', '20260313']; % unique fovs 
 

VISpSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

% 
%VISpDataSpks = getTuningDataByCondition(VISpSessions, 'signalToUse', 'spks');
VISpData = getTuningDataByCondition(VISpSessions, 'signalToUse', 'dFFNeuropilCorrected');


VISpData = appendFilteredROIs(VISpData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);
%VISpDataSpks = appendFilteredROIs(VISpDataSpks,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);


%% plots 

% spatial averaging of activity - to give an overview of where they peak,
% how similar the peaks are relative to the idential landmark and 
alignToPeakTuningCurvePosition_Pooled(RSPData, 'RSP');
alignToPeakTuningCurvePosition_Pooled(VISpData, 'VISp');

% Plot smi cumulative probablity and histograms 
compareAndPlot_SMI_RSP_vs_VISp(RSPData, VISpData)

% plot peak smis for boutons and somas 
% this analysis splits boutons/somas into two smi categories 
% and plots the position of the peak responses as histograms 
plotSMI_PeakDistributions_VISP_RSP(RSPData, VISpData)


% Splits corridor into two segments (30-100; 101-170) and finds the
% preferred peak location. plots as histograms and runs stats 
plotPeaksCorridorSplit(RSPData, VISpData)

% do they have a preference for one or the other landmark identity? 
% this function computs a new index: 
% It finds the preferred peak lap position bin on odd laps; find peak
% activity on even laps. If an roi peaks at grating, it computes the mean
% across the alternative landmark (plaid 1 and plaid 2) and then computes
% fmi = (peak - meanAlternative)/(peak + meanAlternative) 
% bounds the index strictly between 0 and 1.

%1.0 = Perfectly selective and stable.
%0.0 = Not selective or unstable.
%Anything in between = Moderately selective.

plotLandmarkIdentityPreference(RSPData, VISpData)


% autocorrelation across position bins
plotSpatialAutocorrelation_Pooled(RSPData, VISpData);

%

%% test- where are they suppressed? 



[suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPData, 'RSP', ...
    'zThresh', -1, 'windowCm', 10);

% Check a specific bouton
row = suppTable(suppTable.suppLM4 == 1, :);
disp(row)


% VISp - spikes; need to load spike data from another script (temp)
[suppROIs_v1, suppTable_v1, fig2] = findLandmarkSuppressedROIs(VISpDataSpks, 'V1', ...
'zThresh', -1, 'windowCm', 10, 'TypeToPlot', 'Somas', 'daysToPlot',  200 );

% VISp - dff 
[suppROIs_v1_df, suppTable_v1_df, fig3] = findLandmarkSuppressedROIs(VISpData, 'V1', ...
'zThresh', -1, 'windowCm', 10, 'TypeToPlot', 'Somas', 'daysToPlot',  200 );
