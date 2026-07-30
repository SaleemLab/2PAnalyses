% Manipulation corridor 

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.8.0_examples\example_rois_spks';
%%
resp025 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat");
sfi025 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_sessionFileInfo.mat");


%
resp024 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns.mat");
sfi024 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_sessionFileInfo.mat");

figA1 = plotSingleROI_ConditionStack(sfi025.sessionFileInfo, resp025, 1, 'spks', 1, 1);
saveFigureFormats(figA1, fullfile(outputDir, 'visual'));


figA2 = plotSingleROI_ConditionStack(sfi025.sessionFileInfo, resp025, 22, 'spks', 1, 1);
saveFigureFormats(figA2, fullfile(outputDir, 'swap23'));

figA3 = plotSingleROI_ConditionStack(sfi025.sessionFileInfo, resp025, 24, 'spks', 1, 1);

figA4 = plotSingleROI_ConditionStack(sfi025.sessionFileInfo, resp025, 12, 'spks', 1, 1);

figA5 = plotSingleROI_ConditionStack(sfi025.sessionFileInfo, resp025, 10, 'spks', 1, 1);
saveFigureFormats(figA5, fullfile(outputDir, 'no_spatial_structure'));

figA6 = plotSingleROI_ConditionStack(sfi024.sessionFileInfo, resp024, 29, 'spks', 1, 1);
saveFigureFormats(figA6, fullfile(outputDir, 'visual-omission-present'));

figA7 = plotSingleROI_ConditionStack(sfi024.sessionFileInfo, resp024, 58, 'spks', 1, 1);
saveFigureFormats(figA7, fullfile(outputDir, 'swap23_magnitudehigh'));


