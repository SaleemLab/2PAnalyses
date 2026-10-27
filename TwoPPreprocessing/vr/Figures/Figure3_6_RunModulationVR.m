pairs=struct;
pairs.M25132 = ['20260226', '20260228', '20260313'];
pairs.M25133 = '20260224';
pairs.M26003 = ['20260322', '20260324', '20260325'];

ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

RSPData = getTuningDataByCondition(ExpRSPSessions);
% Filter rois based on critera defined; this alos include rois that peak;
% this exlcudes approx 300 rois that peak in the beginning or end of the

% selecting the same rois as the where do they peak figure
RSPData = appendFilteredROIs(RSPData,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);


%% Load running data from an example sessions and plot 

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.6\runBeh_medianAcrossPos_withlapspeedpos';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    

% resp_1 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260220\M25132_20260220_Response_M25132_BaselineCorridor_20260220_00001.mat");
% figA1 = plotSpeedTrajectories(resp_1);
% set(figA1, 'Visible', 'on');
% saveFigureFormats(figA1, fullfile(outputDir, 'runBeh_example1'));

% exclude this
% resp_1 = load("Z:\ibn-vision\DATA\SUBJECTS\M25133\Analysis\20260224\M25133_20260224_Response_M25133_BaselineCorridor_20260224_CombinedRuns.mat");
% figA1 = plotSpeedTrajectories(resp_1);
% set(figA1, 'Visible', 'on');
% saveFigureFormats(figA1, fullfile(outputDir, 'runBeh_example1'));

% this includes low medium high using the other cersion (cant remember
% what) 


% resp_2 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_CombinedRuns.mat");
% figA2 = plotSpeedTrajectories(resp_2);
% set(figA2, 'Visible', 'on');
% saveFigureFormats(figA2, fullfile(outputDir, 'runBeh_example2_m25132_20260226'));
% 
% 
% resp_3 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_BaselineCorridor_20260322_00002.mat");
% figA3 = plotSpeedTrajectories(resp_3);
% set(figA2, 'Visible', 'on');
% saveFigureFormats(figA3, fullfile(outputDir, 'runBeh_example3_m26003_20260322'));
% 
% resp_4 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260313\M25132_20260313_Response_M25132_LandManipCorridor_20260313_CombinedRuns.mat");
% figA4 = plotSpeedTrajectories(resp_4);
% set(figA4, 'Visible', 'on');
% saveFigureFormats(figA4, fullfile(outputDir, 'runBeh_example4_m25132_20260313'));
% 
% 
% 
% figA5 = plotSpeedTrajectories(resp_6);
% set(figA5, 'Visible', 'on');
% saveFigureFormats(figA5, fullfile(outputDir, 'runBeh_example5_m26003_20260325'));


%% Used for thesis to get the speed plots 


resp_6 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat"); % used 
resp_4 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260313\M25132_20260313_Response_M25132_LandManipCorridor_20260313_CombinedRuns.mat");
resp_3 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_BaselineCorridor_20260322_00002.mat"); 
resp_1 = load("Z:\ibn-vision\DATA\SUBJECTS\M25133\Analysis\20260224\M25133_20260224_Response_M25133_BaselineCorridor_20260224_CombinedRuns.mat"); % dont use 
resp_7 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_CombinedRuns.mat"); % used 
resp_8 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns.mat")


% figA2 = plotSpeedTrajectories_positionMedian(resp_7);
% set(figA2, 'Visible', 'on');
% saveFigureFormats(figA2, fullfile(outputDir, 'runBeh_example2_m25132_20260226'));
% 
% 
% 
% figA3 = plotSpeedTrajectories_positionMedian(resp_3);
% set(figA2, 'Visible', 'on');
% saveFigureFormats(figA3, fullfile(outputDir, 'runBeh_example3_m26003_20260322'));
% 
% 
% figA4 = plotSpeedTrajectories_positionMedian(resp_4);
% set(figA4, 'Visible', 'on');
% saveFigureFormats(figA4, fullfile(outputDir, 'runBeh_example4_m25132_20260313'));
% 
% 
% figA5 = plotSpeedTrajectories_positionMedian(resp_6);
% set(figA5, 'Visible', 'on');
% saveFigureFormats(figA5, fullfile(outputDir, 'runBeh_example5_m26003_20260325'));



%% next draw the heatmaps and activity above 
% load example sessions
% resp_4 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_CombinedRuns.mat");
% load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_sessionFileInfo.mat");
% saveSpeedPositionActivity_AllROIs(sessionFileInfo, resp_4)
% clear sessionFileInfo
% 
% resp_5 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_BaselineCorridor_20260322_00002.mat");
% load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_sessionFileInfo.mat")
% saveSpeedPositionActivity_AllROIs(sessionFileInfo, resp_5)

% resp_6 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat");
% think about what examples to show later on 
% load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_sessionFileInfo.mat")
% saveSpeedPositionActivity_AllROIs(sessionFileInfo, resp_6)


savePath = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section3_Fig3.6\speedPosActivity_roisegs';


% No Difference 
figB1 = plotSpeedPositionActivity_ForROI(resp_4,27, 1, 0);
saveFigureFormats(figB1, fullfile(savePath, 'speedPosAct_meanSE_For_speedranges'));

plotSpeedPositionActivity_ForROI(resp_5, 8); % High-plaid; low-grating
plotSpeedPositionActivity_ForROI(resp_5, 10);
plotSpeedPositionActivity_ForROI(resp_4, 335); % High speed preference for
plotSpeedPositionActivity_ForROI(resp_4, 1);  % Low-medium Speed Preferred
plotSpeedPositionActivity_ForROI(resp_4, 417); % Medium Speed Preferred

% 
%saveSpeedPositionAndTuningSummaries(sessionFileInfo, resp_4, respGray.response, respDark.response)

%% low 
% figA2 = plotSpeedTrajectories_positionMedian(resp_7)
% set(figA2, 'Visible', 'on');
% saveFigureFormats(figA2, fullfile(outputDir, 'low_trej'));
% 
% figA22 = plotSpeedPositionActivity_ForROI(resp_7, 1)
% set(figA22, 'Visible', 'on');
% saveFigureFormats(figA22, fullfile(outputDir, 'low_trej_dff'));
%% high
figA3 = plotSpeedTrajectories_positionMedian(resp_7, 1);
set(figA3, 'Visible', 'on');
saveFigureFormats(figA3, fullfile(outputDir, 'high_trej'));

figA33 = plotSpeedPositionActivity_ForROI(resp_7, 335)
set(figA33, 'Visible', 'on');
saveFigureFormats(figA33, fullfile(outputDir, 'high_trej_dff'));
%% keep 


figA4 = plotSpeedTrajectories_positionMedian(resp_6);
set(figA4, 'Visible', 'on');
saveFigureFormats(figA4, fullfile(outputDir, 'similar_trej'));

figA44 = plotSpeedPositionActivity_ForROI(resp_6, 83);
set(figA44, 'Visible', 'on');
saveFigureFormats(figA44, fullfile(outputDir, 'similar_trej_dff'));


%% low

figA5 = plotSpeedTrajectories_positionMedian(resp_8);
set(figA5, 'Visible', 'on');
saveFigureFormats(figA5, fullfile(outputDir, 'low_trej'));

figA55 = plotSpeedPositionActivity_ForROI(resp_8, 29);
set(figA55, 'Visible', 'on');
saveFigureFormats(figA55, fullfile(outputDir, 'low_trej_dff'));

%%




%%
% %% viewAllROIs_SpeedPositionActivity.m
% % Shows plotSpeedPositionActivity_ForROI one ROI at a time. Press any key
% % (with the figure or command window focused) to close the current ROI
% % and move to the next. Nothing is saved.
% 
% numROIs = size(resp_8.speedPositionActivity.lowHigh.matrixLow, 1);
% fprintf('Viewing %d ROIs — press any key to advance, Ctrl+C to stop.\n', numROIs);
% 
% for targetROI = 1:numROIs
%     try
%         figHandle = plotSpeedPositionActivity_ForROI(resp_8, targetROI);
%         fprintf('ROI %d / %d — press any key for next\n', targetROI, numROIs);
%         pause(2)
%         close(figHandle);
%     catch ME
%         fprintf('  Skipped ROI %d — %s\n', targetROI, ME.message);
%     end
% end