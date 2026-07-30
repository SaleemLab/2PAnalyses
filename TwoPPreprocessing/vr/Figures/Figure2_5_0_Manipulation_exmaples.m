% Manipulation corridor 

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_5_0_Section3\example_rois_spks';
%%
resp018_004 = load("Z:\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260318\M26004_20260318_Response_M26004_LandManipCorridor_20260318_CombinedRuns.mat");
sfi018_004 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_sessionFileInfo.mat");


%
resp024 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns.mat");
sfi024 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_sessionFileInfo.mat");

figA1 = plotSingleROI_ConditionStack(sfi018_004.sessionFileInfo, resp018_004, 2, 'spks', 1, 1);
saveFigureFormats(figA1, fullfile(outputDir, 'visual-background'));


figA2 = plotSingleROI_ConditionStack(sfi018_004.sessionFileInfo, resp018_004, 22, 'spks', 1, 1);
saveFigureFormats(figA2, fullfile(outputDir, 'swap23'));

figA3 = plotSingleROI_ConditionStack(sfi018_004.sessionFileInfo, resp018_004, 24, 'spks', 1, 1);

figA4 = plotSingleROI_ConditionStack(sfi018_004.sessionFileInfo, resp018_004, 12, 'spks', 1, 1);

figA5 = plotSingleROI_ConditionStack(sfi018_004.sessionFileInfo, resp018_004, 10, 'spks', 1, 1);
saveFigureFormats(figA5, fullfile(outputDir, 'no_spatial_structure'));

figA6 = plotSingleROI_ConditionStack(sfi024.sessionFileInfo, resp024, 29, 'spks', 1, 1);
saveFigureFormats(figA6, fullfile(outputDir, 'visual-omission-present'));

figA7 = plotSingleROI_ConditionStack(sfi024.sessionFileInfo, resp024, 58, 'spks', 1, 1);
saveFigureFormats(figA7, fullfile(outputDir, 'swap23_magnitudehigh'));


