
% Filter table for EVERYTHING (Boutons AND Somas)
% masterTable = filterMasterTable('MouseID', {'M25131', 'M26005', 'M26004'}, 'Session', {'20260318', '20260321', '20260322'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', 200); %'DayOfExperience', 100
% masterTable = filterMasterTable('MouseID', {'M25132', 'M26003'}, 'Session', {'20260312', '20260313', '20260325', '20260324'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', 200); %'DayOfExperience', 100
masterTable = filterMasterTable('MouseID', {'M25132', 'M25133', 'M26003'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', 5);

% It will use 
allData = getTuningData(masterTable, 'applySmoothing', true); 

% allData = getTuningDataByCondition(masterTable);
% Plot distributions per mouse for all critera (excluding matched bouton
% correlation

% plotConditionTuning(allData, targetArea, days, savePath)
% plotOddEvenBaseline(allData, 'RSP', '200', '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\OddEvenSomaBaselineCondition.png')
% plotDistributionsROISelectionCriteria(allData)

%Plot the main grid
plotConditionTuning(allData, 'RSP', 200, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\BoutonsAcrossCondition.png')

plotTuningHeatmaps(allData, 'DaysToPlot', [1 2 3 4 5], 'TypeToPlot', 'bouton');
plotEvenTuningHeatmaps(allData)
plotOddEvenTuning(allData, 'RSP', 5)

%  Plot Bouton vs Soma for one mouse
plotBoutonSomaComparison(allData, 'M25037', 'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\M25041_Comparison.png');


%  Plot RSP Population Pooled
plotPooledPopulation(allData, 'RSP', 'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\RSP_PooledNewCorridor2026.png', 'DaysToPlot', [5], 'TypeToPlot', 'Boutons');

plotModulation_Histogram(allData, 'DaysToPlot', 1:5, ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\Modulation_Hist_Days1to5_newmice.png');

plotModulation_CDF(allData);