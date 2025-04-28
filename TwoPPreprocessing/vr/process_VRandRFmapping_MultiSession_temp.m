% process_all_mice_sessions.m
% Script to process all mice and sessions

clear; clc;

% Define mouse and session information 
mouseInfo = {
    'M25037', {'20250428'}, {'M25037_VRCorr_20250428_00001', 'M25037_RFMapping_20250428_00001', 'M25037_GrayScreen_20250428_00001'};
    'M25038', {'20250428'}, {'M25038_VRCorr_20250428_00001', 'M25038_RFMapping_20250428_00001', 'M25038_GrayScreen_20250428_00001'};
};

% Processing parameters 
interpRate = 60;        % Hz
rfpreStimTime = 0.5;     % seconds
planeNums = 8;      % original planes in recording
channelsSaved = 2;       % Channel to align frame times to
method = 2;            % Method for PSTH extraction
frameRate = 7.28;       % Degrees/bin or other relevant stimulus unit
applyNeuropilCorrection = true;
calculateDFF = true;
pdthreshold = 5;
rfpostStimTime = 3; 

% Loop through each mouse
for thisMouse = 1:size(mouseInfo,1)
    
    mousenumber = mouseInfo{thisMouse,1};
    sessionNames = mouseInfo{thisMouse,2};
    stimlists = mouseInfo{thisMouse,3};
    fprintf('Processing mouse: %s\n', mousenumber);
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
        vrStimName = stimlists{1};
        rfStimName = stimlists{2};
        fprintf('  Session: %s | Stimulus: %s\n', sessionName, stimName);
        try
            % Get session file paths
            sessionFileInfo = get2PsessionFilePaths(mousenumber, sessionName, {stimName});
            % Load metadata; all general
            sessionFileInfo = get2PMetadata(sessionFileInfo);
            % Get frame times for two channels; all general
            sessionFileInfo = get2PFrameTimes_TwoChannels(sessionFileInfo, planeNums, channelsSaved);
            % Process peripheral files (wheel, position etc); all general
            sessionFileInfo = processPeripheralFiles(sessionFileInfo);
            % Find VRLag
            findBonsaiPeripheralLag(sessionFileInfo, 1, 60, vrStimName);
            % Merge Bonsai and Suite2p files; all general
            sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo);
            %-------- VR ------------
            % Align bonsai to peripheral data 
            [vrbonsaiData, sessionFileInfo]  = alignVRBonsaiToPeripheralData(sessionFileInfo,vrStimName);
            % Interpolate VR data streams 
            [vrprocessedTwoPData, vrbonsaiData, vrperipheralData, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo,interpRate,'TwoPFrameTime', vrStimName, true);
            % Extract VR peripheral data
            [vrresponse, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo,   vrStimName, true);
            % Get twopframe position bins 
            [vrresponse, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, vrStimName);
            % Get lap position activity 
            [vrresponse, sessionFileInfo] = getLapPositionActivityV2(sessionFileInfo, 'F', false, vrStimName , false, applyNeuropilCorrection, calculateDFF);
            % Plot population summary 
            plotSortedPopulationResponse_OddEven(sessionFileInfo, vrresponse, true)

            %--------- RF -------------
            % Get stimulus events from Bonsai file
            [rfbonsaiData, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, rfStimName, 'StimulusParams');
            % Get stimulus times
            [rfbonsaiData, sessionFileInfo] = getStimTimes(sessionFileInfo, rfStimName, pdthreshold);
            % Get 2P frames by trial
            [rfresponse, sessionFileInfo] = get2PFramesByTrialV3(sessionFileInfo, rfStimName, rfpostStimTime, rfpreStimTime);
            % Group trials
            [rfbonsaiData, sessionFileInfo] = getTrialGroupsV2(sessionFileInfo, rfStimName);
            % PSTHs
            [rfresponse, sessionFileInfo] = getTrialResponsePSTHsV4(sessionFileInfo, rfStimName, method, interpRate, frameRate, applyNeuropilCorrection, calculateDFF);
            % Plot individual rois PSTH
            %plotAllROIPSTHsByPosition(sessionFileInfo, rfresponse);
            % Plot per roi 
            plotRFGrid_byPosition_ROIs(sessionFileInfo, rfStimName);
            
            % ---- VR+RF Plotting -----
            % Combined RF and VR roi-wise plot 
            
            fprintf('    Done!\n'); 
        catch ME
            warning('    Error processing %s %s: %s', mousenumber, sessionName, ME.message);
        end
        
    end
end

disp('All mice and sessions processed and figures saved in session-wise analysis folders');
