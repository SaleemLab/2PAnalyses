
% Filter table for EVERYTHING (Boutons AND Somas)
masterTable = filterMasterTable('MouseID', {'M25131', 'M26005', 'M26004'}, 'Session', {'20260305','20260306', '20260311','20260318', '20260321', '20260322'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', 200); %'DayOfExperience', 100
% masterTable = filterMasterTable('MouseID', {'M25132', 'M26003'}, 'Session', {'20260312', '20260313', '20260325', '20260324'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', 200); %'DayOfExperience', 100
%masterTable = filterMasterTable('MouseID', {'M25057', 'M25058'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', [4 5]);

% It will use 
%allData = getTuningData(masterTable, 'applySmoothing', true); 

allData = getTuningDataByCondition(masterTable);
% Plot distributions per mouse for all critera (excluding matched bouton
% correlation

% plotConditionTuning(allData, targetArea, days, savePath)
% plotOddEvenBaseline(allData, 'RSP', '200', '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\OddEvenSomaBaselineCondition.png')
% plotDistributionsROISelectionCriteria(allData)

%Plot the main grid
% plotConditionTuning(allData, 'RSP', 200, '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\BoutonsAcrossCondition.png')

plotTuningHeatmaps(allData, 'DaysToPlot', [4 5]);
plotEvenTuningHeatmaps(allData)
plotOddEvenTuning(allData, 'ENTm', [4 5])

%  Plot Bouton vs Soma for one mouse
plotBoutonSomaComparison(allData, 'M25037', 'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\M25041_Comparison.png');


%  Plot RSP Population Pooled
plotPooledPopulation(allData, 'ENTm', 'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\NewRSPBoutonMice\RSP_PooledNewCorridor2026_entm.png', 'DaysToPlot', [4 5], 'TypeToPlot', 'Boutons');

plotModulation_Histogram(allData, 'DaysToPlot', 1:5, ...
    'SavePath', 'Z:\ibn-vision\USERS\Sonali\Figures\Modulation_Hist_Days1to5_newmice.png');

plotModulation_CDF(allData);

%
pooledRSP_Filtered = pooledRSP;

%
daysToKeep = {'Day5', 'Day200'};
allDaysFound = fieldnames(pooledRSP.RSP.Days);

% 3. Remove all other days (Day 2, Day 4, etc.)
for i = 1:length(allDaysFound)
    if ~ismember(allDaysFound{i}, daysToKeep)
        pooledRSP_Filtered.RSP.Days = rmfield(pooledRSP_Filtered.RSP.Days, allDaysFound{i});
    end
end

% 4. Call the plotting function with the filtered data
plotPooledSMI_CDF(pooledRSP_Filtered, 'RSP', 'Days');

%%
% 1. Setup the figure
figure('Color', 'w', 'Position', [200 200 600 500]);
hold on;

% 2. Define colors for RSP Days
rspColors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.494, 0.184, 0.556]; 
daysToShow = {'Day1', 'Day3', 'Day5'};
legendEntries = [];

% 3. Plot the "Expert" V1 (Day 200) as a Reference
if isfield(pooledVISp.V1.Days, 'Day200')
    v1SMI = pooledVISp.V1.Days.Day200;
    [f, x] = ecdf(v1SMI);
    hV1 = plot(x, f, '--k', 'LineWidth', 2); % Dashed black line
    legendEntries(end+1) = hV1;
    v1Label = sprintf('V1 Expert (n=%d)', length(v1SMI));
else
    v1Label = 'V1 Day 200 not found';
end

% 4. Plot RSP Days 1, 3, and 5
rspLabels = {};
for d = 1:length(daysToShow)
    dayField = daysToShow{d};
    if isfield(pooledRSP.RSP.Days, dayField)
        rspSMI = pooledRSP.RSP.Days.(dayField);
        [f, x] = ecdf(rspSMI);
        h = plot(x, f, '-', 'Color', rspColors(d,:), 'LineWidth', 2.5);
        legendEntries(end+1) = h;
        rspLabels{end+1} = sprintf('RSP %s (n=%d)', dayField, length(rspSMI));
    end
end


line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
line([-1 1], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
grid on; xlim([-1 1]); ylim([0 1]);
xlabel('Spatial Modulation Index (SMI)');
ylabel('Cumulative Probability');
title('RSP vs VISp');

% 6. Legend
legend(legendEntries, [{v1Label}, rspLabels], 'Location', 'southeast', 'Box', 'off');

%%
% 1. Setup the figure
figure('Color', 'w', 'Position', [200 200 600 500]);
hold on;

% 2. Define colors for RSP Days
rspColors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.494, 0.184, 0.556]; 
daysToShow = {'Day1', 'Day3', 'Day5'};
legendEntries = [];
displayLabels = {};

% Helper to clean data (Removes NaNs and Infs)
cleanData = @(data) data(isfinite(data) & data < 0.999);

% 3. Plot the pooled V1 Total as a Reference
if isfield(pooledVISp.V1, 'AllSMI')
    v1SMI = cleanData(pooledVISp.V1.AllSMI); % Clean data
    [f, x] = ecdf(v1SMI);
    hV1 = plot(x, f, '--k', 'LineWidth', 2); 
    legendEntries(end+1) = hV1;
    displayLabels{end+1} = sprintf('VISp > Day5 (n=%d)', length(v1SMI));
end

% 4. Plot the pooled RSP Total
if isfield(pooledRSP.RSP, 'AllSMI')
    rspTotalSMI = cleanData(pooledRSP.RSP.AllSMI); % Clean data
    [f, x] = ecdf(rspTotalSMI);
    hRSPTotal = plot(x, f, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2); 
    legendEntries(end+1) = hRSPTotal;
    displayLabels{end+1} = sprintf('RSP > Day5 (n=%d)', length(rspTotalSMI));
end



% 6. Formatting
line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
line([-1 1], [0.5 0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
grid on; 
xlim([-1 1]); 
ylim([0 1]);
xlabel('Spatial Modulation Index (SMI)');
ylabel('Cumulative Probability');
title('RSP vs VISp');

% 7. Legend
if ~isempty(legendEntries)
    legend(legendEntries, displayLabels, 'Location', 'southeast', 'Box', 'off');
end