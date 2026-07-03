% Fig 3.3: Plot ROI examples chosen from the 25th, 50th and 75th quantile
% of the SMI distribution; And overall odd and even population heatmaps 

pairs=struct;
pairs.M25132 = ['20260226', '20260228', '20260313'];
pairs.M25133 = '20260224';
pairs.M26003 = ['20260322', '20260324', '20260325'];

ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

RSPData = getTuningDataByCondition(ExpRSPSessions, 'signalToUse');
% Filter rois based on critera defined; here i include all rois or only the
% ones that do not have peak responses in the edges?

RSPData = appendFilteredROIs(RSPData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);

plotPooledPopulation_OddEven(RSPData, 'RSP', ...
    'DaysToPlot',[5 200], ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.3\OddEvenPop\odd_even_stableBoutons');




%% 
pooledSMIs = [];
pooledROIs = [];
pooledMice = {};
pooledSessions = {};
pooledSessIndices = []; % 

for iSess = 1:length(RSPData)
    thisSession = RSPData(iSess);
    if isempty(thisSession.FilteredROIs), continue; end
    
    keptROIs = thisSession.FilteredROIs(:);
    smiValues = thisSession.SMI.SMI(keptROIs);
    
    validMask = ~isnan(smiValues);
    if ~any(validMask), continue; end
    
    pooledSMIs        = [pooledSMIs; smiValues(validMask)];
    pooledROIs        = [pooledROIs; keptROIs(validMask)];
    pooledMice        = [pooledMice; repmat({thisSession.MouseID}, sum(validMask), 1)];
    pooledSessions    = [pooledSessions; repmat({thisSession.Session}, sum(validMask), 1)];
    pooledSessIndices = [pooledSessIndices; repmat(iSess, sum(validMask), 1)]; % NEW
end

%%  find global quantiles and list alternative candidates
if ~isempty(pooledSMIs)
    sortedSMIs = sort(pooledSMIs);
    N = length(sortedSMIs);
    
    idx25 = max(1, round(0.25 * N));
    idx50 = max(1, round(0.50 * N));
    idx75 = max(1, round(0.80 * N));
    
    qVals = [sortedSMIs(idx25), sortedSMIs(idx50), sortedSMIs(idx75)];
    quantLabels = {'25th Pct (Low SMI)', '50th Pct (Mid SMI)', '75th Pct (High SMI)'};
    
    % Define how close an alternative SMI must be to the exact target value
    smiTolerance = 0.01; 
    
    fprintf('=== QUANTILE TARGET CANDIDATES ===\n\n');
    for q = 1:3
        targetSMI = qVals(q);
        
        % Find all indices within the tolerance window
        candidateIndices = find(abs(pooledSMIs - targetSMI) <= smiTolerance);
        
        % Sort candidates so the ones closest to the true quantile value appear first
        [~, sortOrder] = sort(abs(pooledSMIs(candidateIndices) - targetSMI));
        candidateIndices = candidateIndices(sortOrder);
        
        fprintf('=========================================\n');
        fprintf('  %s (Target SMI: %.3f)\n', quantLabels{q}, targetSMI);

        
        % Display up to 10 alternative options for you to inspect
        numOptionsToDisplay = min(10, length(candidateIndices));
        for iOpt = 1:numOptionsToDisplay
            idx = candidateIndices(iOpt);
            
            fprintf('  Option %d [SMI: %.3f]:\n', iOpt, pooledSMIs(idx));
            fprintf('     Session Index: %d\n', pooledSessIndices(idx));
            fprintf('     ROI Number:    %d\n', pooledROIs(idx));
            fprintf('     Mouse | Sess:  %s | %s\n\n', pooledMice{idx}, pooledSessions{idx});
        end
    end
else
    warning('No valid filtered SMI values found across any loaded sessions.');
end
%% plot them 
% figLow = plotRoiSpatialTuning(RSPData, 2, 271, 'Low_SMI');
figLow = plotRoiSpatialTuning(RSPData, 2,271, 'Low_SMI');
saveFigureFormats(figLow, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.3\exampleROIs\ExampleTuning_Low');

figMid = plotRoiSpatialTuning(RSPData, 5, 29, 'Middle_Significant');
saveFigureFormats(figMid, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.3\exampleROIs\ExampleTuning_Middle');

figHigh = plotRoiSpatialTuning(RSPData, 1, 81, 'High_Significant');
saveFigureFormats(figHigh, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.3\exampleROIs\ExampleTuning_High');

