pairs=struct; 
pairs.M26005 = ['20260305', '20260306', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M25126 = ['20260311', '20260312', '20260313']; % unique fovs 
 

VISpSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

% edit this function to include spikes
VISpDataSpks = getTuningDataByCondition(VISpSessions, 'signalToUse', 'spks');

% these critera were included to plot all rois 
VISpDataSpks = appendFilteredROIs(VISpDataSpks, 'UseExpVar_SigNullDist', true, 'cvExpvarThreshold', 0.1, 'ExpVarSigThreshold', 0.01, 'FilterSomasByRF', false);

% filter used to plot the 3 conditions summary for aman
% VISpDataSpks = appendFilteredROIs(VISpDataSpks, 'UseExpVar_SigNullDist', true, 'cvExpvarThreshold', 0.1, 'ExpVarSigThreshold', 0.01, 'UseHalves', true, 'RhoHalvesThreshold', 0.8, 'FilterSomasByRF', true);

% 
% plotPooledPopulation_AcrossConditions(VISpDataSpks, 'V1', ...
%     'TypeToPlot', 'Somas', ...
%     'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_5_Section3\pooled_5conditions\pooled_5conditions');

plotConditions_DifferenceIncluded(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_5_Section3\pooled_5conditions_withdiff\pooled_5conditions_withdiff');


% plotBackgroundROIs_DifferenceOnly(VISpDataSpks, 'V1', ...
%     'TypeToPlot', 'Somas', ...
%     'RequiredConditions', {'Swap_2_3', 'Swap_3_4', 'Omit_2', 'Omit_3', 'Omit_4'},  ...
%     'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_5_Section3\pooled_5conditions_withDiff\pooled_5conditions_withDiff');

