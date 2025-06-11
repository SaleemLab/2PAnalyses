% process_all_mice_sessions.m
% Script to process all mice and sessions

clear; clc;

% Define mouse and session information 
mouseInfo = {
%      
%      'M25057', {'20250604'}, {'M25057_VRCorr_20250604_00001',...
%                                'M25057_RFMapping_20250604_00001',...
%                               'M25057_GrayScreen_20250604_00001',...
%                               'M25057_DirTuning_20250604_00001',...
%                               'M25057_CheckerBoard_20250604_00001',...
%                               'M25057_DotMotion_SpeedTuning_20250604_00001',...
%                               'M25057_DriftingBar_20250604_00001',...
%                               'M25057_FullFieldFlash_20250604_00001'...
%      };

     'M25058', {'20250604'}, {'M25058_VRCorr_20250604_00001',...
                            'M25058_VRCorr_20250604_00002',...
                          'M25058_RFMapping_20250604_00001',...
                          'M25058_GrayScreen_20250604_00001',...
                          'M25058_DirTuning_20250604_00001',...
                          'M25058_CheckerBoard_20250604_00001',...
                          'M25058_DotMotion_SpeedTuning_20250604_00001',...
                          'M25058_DriftingBar_20250604_00002',...
                          'M25058_FullFieldFlash_20250604_00001'...
     };

%       'M25040', {'20250523'}, {'M25040_VRCorr_20250523_00001',...
%                             'M25040_RFMapping_20250523_00001',...
%                             'M25040_GrayScreen_20250523_00001',...
%                             'M25040_DirTuning_20250523_00001',...
%                             'M25040_CheckerBoard_20250523_00001',...
%                             'M25040_DotMotion_SpeedTuning_20250523_00001',...
%                             'M25040_DriftingBar_20250523_00001',...
%                             'M25040_FullFieldFlash_20250523_00001',...
%              };

    
%        'M25012', {'20250522'}, {'M25012_VRCorr_20250522_00001',...
%                             'M25012_VRCorr_20250522_00002',...
%                             'M25012_RFMapping_20250522_00001',...
%                             'M25012_GrayScreen_20250522_00001',...
%                             'M25012_DirTuning_20250522_00001',...
%                             'M25012_CheckerBoard_20250522_00001',...
%                             'M25012_DotMotion_SpeedTuning_20250522_00001',...
%                             'M25012_DriftingBar_20250522_00001',...
%                             'M25012_FullFieldFlash_20250522_00001',...
%        
%        };

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
pdthreshold = 10;
rfpostStimTime = 3; 

% Loop through each mouse
for thisMouse = 1:size(mouseInfo,1)
    
    mousenumber = mouseInfo{thisMouse,1};
    sessionNames = mouseInfo{thisMouse,2};
    stimlists = mouseInfo{thisMouse,3};
    fprintf('Processing mouse: %s\n', mousenumber);
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};
%         if strcmp(mousenumber, 'M25040') || strcmp(mousenumber, 'M25037')
        if mousenumber == 'M25040'
            vrStimNames = stimlists(1);
            rfStimName = stimlists{2};
        elseif mousenumber == 'M25057'
            vrStimNames = stimlists(1);
            rfStimName = stimlists{2};
        elseif mousenumber ==  'M25058'
            vrStimNames = stimlists(1:2);
            rfStimName = stimlists{3};
        end 
        
        disp(['  Session: %s | Stimulus: %s\n', sessionName]);
        try

            % Get session file paths
            sessionFileInfo = get2PsessionFilePaths(mousenumber, sessionName, stimlists,1);
            % Load metadata; all general
            sessionFileInfo = get2PMetadata(sessionFileInfo);
            % Get frame times for two channels; all general
            sessionFileInfo = get2PFrameTimes_TwoChannels(sessionFileInfo, planeNums, channelsSaved);
            % Process peripheral files (wheel, position etc); all general
            sessionFileInfo = processPeripheralFiles(sessionFileInfo);
            % Merge Bonsai and Suite2p files; all general
            sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo);

            for thisVRStim = 1:length(vrStimNames)
                
                vrStimName = vrStimNames(thisVRStim);
                if iscell(vrStimName)
                    vrStimName = vrStimName{1};  % extract the string from the cell
                else
                    vrStimName = vrStimName;    
                end
            
                %-------- VR ------------
                % Get bonsai file 
                [bonsaiData, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, vrStimName);
                % Find VRLag
                [bonsaiData, sessionFileInfo]=findBonsaiPeripheralLag(sessionFileInfo, 1, 60, vrStimName);
                % Align bonsai to peripheral data 
                [bonsaiData, sessionFileInfo]  = alignVRBonsaiToPeripheralData(sessionFileInfo,vrStimName);
                % Interpolate VR data streams 
                [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo,60,'TwoPFrameTime', vrStimName, true);
                [response, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo,   vrStimName, true);
                % Get twopframe position bins 
                [response, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, vrStimName);
                % Get lap position activity 
                [response, sessionFileInfo] = getLapPositionActivityV2(sessionFileInfo, 'F', false, vrStimName , false, true, true);
                % Plot population summary 
                plotSortedPopulationResponse_OddEven(sessionFileInfo, response, true)
                % Single neuron summary plot
                plotAllNeuronSummariesToPDF(sessionFileInfo, response, true)
            end
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
            plotAllROIPSTHsByPosition(sessionFileInfo, rfresponse);
            % Plot per roi 
            plotRFGrid_byPosition_ROIs(sessionFileInfo, rfStimName);
            
            % ---- VR+RF Plotting -----

            % Combined RF and VR roi-wise summary plot
            %plotROISummary_VRAndRFHeatmaps(sessionFileInfo, vrStimName, rfStimName)
            
            fprintf('    Done!\n'); 
        catch ME
            warning('    Error processing %s %s: %s', mousenumber, sessionName, ME.message);
        end
        
    end
end

disp('All mice and sessions processed and figures saved in session-wise analysis folders');
