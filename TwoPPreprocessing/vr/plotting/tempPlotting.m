
% Filter table for EVERYTHING (Boutons AND Somas)
masterTable = filterMasterTable('Exclude', 0, 'Suite2PPreprocessing', 1, 'TargetArea', 'RSP',  'DayOfExperience', [1,2,3,4,5]); %'DayOfExperience', 100


% It will use 
allData = getTuningData(masterTable, 'applySmoothing', true); 


% Plot distributions per mouse for all critera (excluding matched bouton
% correlation
plotDistributionsROISelectionCriteria(allData)

%Plot the main grid
plotTuningHeatmaps(allData);

%  Plot Bouton vs Soma for one mouse
plotBoutonSomaComparison(allData, 'M25037', 'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\M25041_Comparison.png');


%  Plot RSP Population Pooled
plotPooledPopulation(allData, 'RSP', 'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\RSP_Pooled100somasDataClub2026.png', 'DaysToPlot', 100, 'TypeToPlot', 'Somas');

plotModulation_Histogram(allData, 'DaysToPlot', 1:5, ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\Modulation_Hist_Days1to5.png');

plotModulation_CDF(allData, 'DaysToPlot', [1,3,5], 'TypeToPlot', 'Boutons', ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\Modulation_CDF_Boutons_D1vsD5.png');