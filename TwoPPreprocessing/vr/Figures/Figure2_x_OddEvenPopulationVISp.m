pairs=struct; 
pairs.M26005 = ['20260305', '20260306', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M26004 = ['20260305', '20260307', '20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs
pairs.M25131 = ['20260312', '20260313', '20260314', '20260318', '20260321', '20260322']; % unique fovs 
pairs.M25126 = ['20260311', '20260312', '20260313']; % unique fovs 
 


VISpSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
VISpData = getTuningDataByCondition(VISpSessions);


% Filter rois based on critera defined 

VISpData = appendFilteredROIs(VISpData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', false);

plotPooledPopulation_OddEven(VISpData, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Section3_Fig2.3\OddEvenPop-Baseline\odd_even_stableSomas_dffNeu');


%% spikes 
VISpDataSpks = getTuningDataByCondition(VISpSessions, 'signalToUse', 'spks');

VISpDataSpks = appendFilteredROIs(VISpDataSpks, 'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);

plotPooledPopulation_OddEven(VISpDataSpks, 'V1', ...
    'TypeToPlot', 'Somas', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Section3_Fig2.3\OddEvenPop-Baseline\odd_even_stableSomas_spks');


