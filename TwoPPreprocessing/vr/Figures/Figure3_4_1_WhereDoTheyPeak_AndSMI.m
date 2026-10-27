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
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\smi_RSP_VISP_dff';
compareAndPlot_SMI_RSP_vs_VISp(RSPData, VISpData, outputDir)



%% EXAMPLES FOR smi rsp and v1 

figLow = plotSMIExampleROI(RSPData,'Percentile', 25);
% figLow = plotRoiSpatialTuning(VISpDataSpks, 22, 11, 'Low_SMI');
saveFigureFormats(figLow, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SMI\ExampleTuning_25th_RSP');

figMid = plotSMIExampleROI(RSPData,'Percentile', 50 ); % good but slightly lower than the 50th
saveFigureFormats(figMid, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SMI\ExampleTuning_50th_RSP');

figHigh = plotSMIExampleROI(RSPData,'Percentile', 75 );
saveFigureFormats(figHigh, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SMI\ExampleTuning_75th_RSP');


figLowV1 = plotSMIExampleROI(VISpData,'Percentile', 25);
% figLow = plotRoiSpatialTuning(VISpDataSpks, 22, 11, 'Low_SMI');
saveFigureFormats(figLowV1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SMI\ExampleTuning_25th_V1');

figMidV1 = plotSMIExampleROI(VISpData,'Percentile', 50 ); % good but slightly lower than the 50th
saveFigureFormats(figMidV1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SMI\ExampleTuning_50th_V1');

figHighV1 = plotSMIExampleROI(VISpData,'Percentile', 75 );
saveFigureFormats(figHighV1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SMI\ExampleTuning_75th_V1');

%%

% plot smr cumulative and histogram 
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\smr_RSP_VISP_dffNeu';
compareAndPlot_SMR_RSP_vs_VISp(RSPData, VISpData, outputDir)

% plot peak smis for boutons and somas 
% this analysis splits boutons/somas into two smi categories 
% and plots the position of the peak responses as histograms 
plotSMI_PeakDistributions_VISP_RSP(RSPData, 'RSP')


% Splits corridor into two segments (30-100; 101-170) and finds the
% preferred peak location. plots as histograms and runs stats 
plotPeaksCorridorSplit(RSPData, VISpData)

% Original behavior, all cells with a valid peak (no SMI filtering):
plotPeaksCorridorSplit(RSPData, VISpData);

% Restrict to modulated cells only, e.g. SMI > 0.25:
plotPeaksCorridorSplit(RSPData, VISpData, 'SMIThreshold', 0.25);

% Or the stricter cutoff, matching Diamanti's paper:
plotPeaksCorridorSplit(RSPData, VISpData, 'SMIThreshold', 0.5);

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

lRSP = plotSIExampleROI(RSPData,'Percentile', 25);
saveFigureFormats(lRSP, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SI\ExampleTuning_25th_RSP');

MRSP = plotSIExampleROI(RSPData,'Percentile', 50 ); % good but slightly lower than the 50th
saveFigureFormats(MRSP, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SI\ExampleTuning_50th_RSP');

HRSP = plotSIExampleROI(RSPData,'Percentile', 75 );
saveFigureFormats(HRSP, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SI\ExampleTuning_75th_RSP');


lV1 = plotSIExampleROI(VISpData,'Percentile', 25);
% figLow = plotRoiSpatialTuning(VISpDataSpks, 22, 11, 'Low_SMI');
saveFigureFormats(lV1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SI\ExampleTuning_25th_V1');

MV1 = plotSIExampleROI(VISpData,'Percentile', 50 ); % good but slightly lower than the 50th
saveFigureFormats(MV1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SI\ExampleTuning_50th_V1');

HV1 = plotSIExampleROI(VISpData,'Percentile', 75 );
saveFigureFormats(HV1, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\exampleROIs_SI\ExampleTuning_75th_V1');



landmarkCentres  = [40, 80, 120, 160];
gratingLandmarks = [40, 120];
plaidLandmarks   = [80, 160];
tolerance = 15;        % Max distance from center to count as a landmark peak (+/- 15cm)
windowRadius = 15;     % Radius around centers to measure activity (+/- 15cm)
plotExampleLandmarkTraces(RSPData, landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius, 'RSP');
plotExampleLandmarkTraces(VISpData, landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius, 'VISp');
plotExampleLandmarkTraces(RSPData, landmarkCentres, gratingLandmarks, plaidLandmarks, tolerance, windowRadius, 'RSP', [0 0 0 0.5 1], 1);
% autocorrelation across position bins
plotSpatialAutocorrelation_Pooled(RSPData, VISpData);

%% test- where are they suppressed? 



[suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPData, 'RSP', ...
    'zThresh', -0.5, 'windowCm', 15);

% Check a specific bouton
row = suppTable(suppTable.suppLM4 == 1, :);
disp(row)


% VISp - spikes; need to load spike data from another script (temp)
[suppROIs_v1, suppTable_v1, fig2] = findLandmarkSuppressedROIs(VISpDataSpks, 'V1', ...
'zThresh', -1, 'windowCm', 10, 'TypeToPlot', 'Somas', 'daysToPlot',  200 );

% VISp - dff 
[suppROIs_v1_df, suppTable_v1_df, fig3] = findLandmarkSuppressedROIs(VISpData, 'V1', ...
'zThresh', -1, 'windowCm', 10, 'TypeToPlot', 'Somas', 'daysToPlot',  200 );


%%
pooledSMI = [];
pooledRp  = [];
pooledRn  = [];

for s = 1:length(VISpData)
    sess = VISpData(s);
    if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
        smiVals = sess.SMI.SMI(sess.FilteredROIs);
        rpVals  = sess.SMI.Rp(sess.FilteredROIs);   % <-- check this field name
        rnVals  = sess.SMI.Rn(sess.FilteredROIs);   % <-- check this field name

        validMask = ~isnan(smiVals);
        pooledSMI = [pooledSMI; smiVals(validMask)];
        pooledRp  = [pooledRp;  rpVals(validMask)];
        pooledRn  = [pooledRn;  rnVals(validMask)];
    end
end

nTotal = length(pooledSMI);
nNegRp = sum(pooledRp < 0);
nNegRn = sum(pooledRn < 0);

fprintf('RSP boutons pooled: n = %d\n', nTotal);
fprintf('Negative Rp: %d/%d (%.1f%%)\n', nNegRp, nTotal, nNegRp/nTotal*100);
fprintf('Negative Rn: %d/%d (%.1f%%)\n', nNegRn, nTotal, nNegRn/nTotal*100);

outOfRange = abs(pooledSMI) > 1;
fprintf('SMI outside [-1,1]: %d/%d (%.1f%%)\n', sum(outOfRange), nTotal, sum(outOfRange)/nTotal*100);
fprintf('SMI actual range: %.3f to %.3f\n', min(pooledSMI), max(pooledSMI));

