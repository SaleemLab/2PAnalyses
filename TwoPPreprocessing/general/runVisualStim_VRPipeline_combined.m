% process_all_mice_sessions_combined_final.m
% Combined script with dynamic parameter loading based on TypeImaged from filteredTable.
clear; clc;
%% DEFINE MICE, SESSIONS, AND KEYWORDS
% Define all mice and their sessions to be processed
% The mouseInfo is generated from the filteredTable
vrKeywords = {'VRCorr', 'BaselineCorridor', 'LandManipCorridor'};
rfKeywords = {'RFMapping'};
doNotCombine = {'M25040_VRCorr_20250507_00001', 'M25040_VRCorr_20250507_00002', 'M25057_VRCorr_20250526_00001', 'M25057_VRCorr_20250526_00002', 'M25126_VRCorr_20260123_00001', 'M25126_VRCorrBaseline_20260123_00002', 'M25126_VRCorrWithManipulations_20260123_00003', 'M25132_BaselineCorridor_20260219_00001', 'M25132_BaselineCorridor_20260219_00002', ...
    'M26003_BaselineCorridor_20260322_00001',  'M26003_BaselineCorridor_20260322_00002','M26003_BaselineCorridor_20260324_00001', 'M25131_BaselineCorridor_20260421_00001', 'M26005_BaselineCorridor_20260421_00001'};
%changed to dff 2026 jan
% filteredTable now holds the key metadata (TypeImaged) needed for parameter setting.

filteredTable = filterMasterTable('Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'MouseID', {'M25131'}, ... 
    'Session', {'20260423'}); %rerun '20260318' m25131 (checking eye tracking bits here); m26003 20260321 (too)
mouseInfo = sessionsToProcess(filteredTable);
% newOrder = [3, 1, 2, 4, 5];
responseName = {'LandManipCorridor','LandManipCorridor','LandManipCorridor','LandManipCorridor'};
% 
% mouseInfo = mouseInfo(newOrder, :);



totalSessionsToProcess = 0;
sessionsProcessedCount = 0;

for i = 1:size(mouseInfo, 1)
    totalSessionsToProcess = totalSessionsToProcess + length(mouseInfo{i, 2});
end
%% DEFINE PARAMETER SETS BASED ON IMAGING TYPE
% Parameters for Somas Imaging;
paramsSomas.interpRate = 60;
%paramsSomas.frameRate = 30; %changed
paramsSomas.pdthreshold = 10;
paramsSomas.isZcorrected = false;
paramsSomas.zScoreProcessedSignals = true;
paramsSomas.useZScoredProcessedSignals = false; 
paramsSomas.applyTemporalSmoothing = true;
paramsSomas.prctlF = 8; % The percentile from which to take F0 (baseline F).
paramsSomas.windowSize = 60; % The rolling window over which to calculate F0.
paramsSomas.signalName = 'dFF';


% Parameters for Bouton Imaging
paramsBoutons.interpRate = 60;
%paramsBoutons.frameRate = 7.28;
paramsBoutons.pdthreshold = 10;
paramsBoutons.isZcorrected = true;
paramsBoutons.zScoreProcessedSignals = true;
paramsBoutons.useZScoredProcessedSignals = true;
paramsBoutons.applyTemporalSmoothing = true;
paramsBoutons.prctlF = 8; % The percentile from which to take F0 (baseline F).
paramsBoutons.windowSize = 60; % The rolling window over which to calculate F0.
paramsBoutons.signalName = 'dFFNeuropilCorrected';

% Common processing parameters
rfpreStimTime = 2;    % seconds
rfpostStimTime = 3;     % seconds
method = 2;             % Method for PSTH extraction (e.g., 2 for mean)

%% INITIALISE ERROR LOG
logFilePath = fullfile('Z:\ibn-vision\USERS\Sonali\errorLogs', 'runNewVRCorr__20260316.csv');
logHeaders = {'Timestamp', 'Mouse', 'Session', 'ErrorMessage', 'Function', 'LineNumber'};
% If the log file doesn't exist, create it with headers
if ~exist(logFilePath, 'file')
    logTable = table('Size', [0, numel(logHeaders)], 'VariableTypes', repmat("string", 1, numel(logHeaders)), 'VariableNames', logHeaders);
    writetable(logTable, logFilePath);
end
fprintf('Error logging enabled. Log will be saved to: %s\n', logFilePath);

%% START PROCESSING LOOP
% Loop through each mouse defined in mouseInfo
for thisMouse = 1:size(mouseInfo, 1)

    mousenumber = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};

    fprintf('Processing Mouse: %s\n', mousenumber);
    % Loop through each session for the current mouse
    for thisSession = 1:length(sessionNames)
        sessionName = sessionNames{thisSession};

        fprintf('\n-- Processing Session: %s --\n', sessionName);

        try

            sessionsProcessedCount = sessionsProcessedCount + 1;
            fprintf(' --------------- Session: %d of %d ------------ \n', ...
                sessionsProcessedCount, totalSessionsToProcess);

            % Find the row in filteredTable corresponding to the current mouse and session
            sessionIdx = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);

            if isempty(sessionIdx)
                error('Session info not found in filteredTable for %s, %s.', mousenumber, sessionName);
            end

            typeImaged = filteredTable.TypeImaged{sessionIdx}; % Get the imaging type



            fprintf('  Type Imaged: %s. Loading specific processing parameters.\n', typeImaged);

            % Assign parameters based on the imaging type
            if strcmpi(typeImaged, 'Somas')
                currentParams = paramsSomas;
            elseif strcmpi(typeImaged, 'Boutons')
                currentParams = paramsBoutons;
            else
                error('Unknown TypeImaged: %s. Cannot set processing parameters.', typeImaged);
            end


            interpRate = currentParams.interpRate;
            zScoreProcessedSignals = currentParams.zScoreProcessedSignals;
            useZScoredProcessedSignals = currentParams.useZScoredProcessedSignals; 
            isZcorrected = currentParams.isZcorrected;
            applyTemporalSmoothing = currentParams.applyTemporalSmoothing;
            prctlF = currentParams.prctlF;
            windowSize = currentParams.windowSize;
            pdThresholdForStimEvents = currentParams.pdthreshold;
            signalName = currentParams.signalName;

            % Get the list of stimuli for this specific session from
            % Suite2p; These stim list will be extracted
            % the order of stimuli that stimuli were concatenated before
            % feeding into Suite2p.
            stimList = getStimList(mousenumber, sessionName); % If the names of tif files have been changed after running through suite2p this will not work!

            fprintf('  Found stimuli: %s\n', strjoin(stimList, ', '));

            % Identify which stimuli are VR and which are others (including RF)
            vrStimNames = {};
            otherVisualStim = {};
            rfStimNames = {};
            % Will currently pull out all Corr stim names; baseline and
            % manipulations grouped together
            for thisStim = 1:length(stimList)
                if any(contains(stimList{thisStim}, vrKeywords, 'IgnoreCase', true))
                    vrStimNames{end+1} = stimList{thisStim};
                else
                    otherVisualStim{end+1} = stimList{thisStim};
                    if any(contains(stimList{thisStim}, rfKeywords, 'IgnoreCase', true))
                        rfStimNames{end+1} = stimList{thisStim};
                    end
                end
            end

            % A. General file processing for the entire session; includes
            % all stimuli TODO: Integrate sparsenoise
            sessionFileInfo = get2PsessionFilePaths(mousenumber, sessionName, stimList, 1); % Overwrite is a must @Sonali cant rememeber why; take a look..
            sessionFileInfo = get2PMetadata(sessionFileInfo); % This function will not run as the tif files have been moved to a different repo.
            [sessionFileInfo] = get2PFrameTimes_SpeedyVersion(sessionFileInfo, isZcorrected); % Uses dynamic planeNums, isZcorrected
            sessionFileInfo = processPeripheralFiles(sessionFileInfo);
            sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo);
            % [sessionFileInfo] = computeNeuropilCorrectionAndDFF_OnRawTraces(sessionFileInfo);

            
            [sessionFileInfo] = createSessionROIData(sessionFileInfo);
            % add eye tracking to bonsai.mat 


            % B. Process other visual stimuli that do not have bonsai data
            % saved like grayscreen or darkness..
            if ~isempty(otherVisualStim)

                ovIdx = find(~contains(otherVisualStim, {'RFMapping', 'DotMotion_SpeedTuning', 'archived'}, 'IgnoreCase', true) );
                if any(ovIdx)
                    disp('Found other visual stimulus file(s).\n');
                    for thisOtherVisualStim = 1:length(ovIdx)
                        otherVisualStimName = otherVisualStim{ovIdx(thisOtherVisualStim)};
                        fprintf('Processing Stimulus file: %s\n', otherVisualStimName);

                        try
                            % Loads Bonsai data files where appropriate
                            [~, sessionFileInfo]       = getTuningStimEventsBonsaiFile(sessionFileInfo, otherVisualStimName, true);
                            
                            [~, ~, ~, sessionFileInfo] = ...
                                resamplAndAlignVisualStim_BonsaiPeripheralSuite2P(sessionFileInfo,interpRate,'TwoPFrameTime', otherVisualStimName);  % Plotflag false, trimNaNs true and overwrite true;
                            [~, sessionFileInfo]       = computeNeuropilCorrectionAndDFF(sessionFileInfo, otherVisualStimName, zScoreProcessedSignals); % Overwrite is true
                            %[response] = getRunningSpeedTuningCurves(sessionFileInfo, stimName, useZScoredSignals, shuffle);
                            [~, sessionFileInfo] = getStimTimes(sessionFileInfo, otherVisualStimName, pdThresholdForStimEvents);
                        catch
                            fprintf('Missing bonsai data structure for stimulus file: %s\n', otherVisualStimName)
                        end
                    end
                end
            end

            % B.2 Process SparseNoiseTexture
            if ~isempty(otherVisualStim)
                snIdx = find(contains(otherVisualStim, 'SparseNoiseTexture', 'IgnoreCase', true));
                if ~isempty(snIdx)
                    thisSNName = otherVisualStim{snIdx(1)};
                    fprintf('Found Sparse Noise at index %d: %s\n', snIdx(1), thisSNName);
                    [~, sessionFileInfo] = getSNTFramesByTrial(sessionFileInfo, thisSNName, false); %Defaults for sparseNoise have been set tp 0 prestim and 3s post stim
                    [~, sessionFileInfo] = analyseSparseNoise(sessionFileInfo, signalName, 0); %plot flag is 0; use the same signal for SNT as for the vr
                end
            end


            %% B.3 Process RFMapping, Dot Fields, and DirTuning: Tuning Stimuli 
            if ~isempty(otherVisualStim)
        
                tuningStimIdx = find(contains(otherVisualStim, {'RFMapping', 'DotMotion_SpeedTuning', 'dirTuning'}, 'IgnoreCase', true));

                if ~isempty(tuningStimIdx)
                    %loop loop
                    for i = 1:length(tuningStimIdx)
                        thisName = otherVisualStim{tuningStimIdx(i)};
                        fprintf('Processing stimulus: %s\n', thisName);

                        % common processing 
                        [~, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, thisName, true);

                        [~, ~, ~, sessionFileInfo] = resamplAndAlignVisualStim_BonsaiPeripheralSuite2P(sessionFileInfo, interpRate, 'TwoPFrameTime', thisName);
                        % 20s window; 5th percentile to compute f0 
                        [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, thisName, zScoreProcessedSignals,false, 5,5,20);
                        [~, sessionFileInfo] = getStimTimes(sessionFileInfo, thisName, pdThresholdForStimEvents);
                       
                     
                        % specififc analysis
                        if contains(thisName, 'RFMapping')
                            fprintf('  Running RFMapping Analysis...\n');
                            [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                            [~, sessionFileInfo] = getRFMappingTrialGroups(sessionFileInfo, thisName);
                            [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                            % plotRFGrid_byPosition_ROIs(sessionFileInfo, thisName);
                            [sessionFileInfo, RFMapping, RFMappingMetadata, allCenters] = analyseRFMapping(sessionFileInfo, thisName);
                            % plotRFMapping(sessionFileInfo, RFMapping,RFMappingMetadata,false,false)

                        % elseif contains(thisName, 'DotMotion_RFMapping')
                        %     fprintf('  Running DotMotion RFMapping Analysis...\n');
                        %     [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                        %     [bonsaiData, sessionFileInfo] = getTrialGroups(sessionFileInfo, thisName);
                        %     [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                        %     % plotDotFieldRFs(sessionFileInfo, doSmooth)

                        elseif contains(thisName, 'DotMotion_SpeedTuning')
                            fprintf('  Running DotMotion Analysis...\n');
                            [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                            [bonsaiData, sessionFileInfo] = getTrialGroups(sessionFileInfo, thisName);
                            [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                            % plotSpeedTuning(sessionFileInfo, response)
                            %plotSpeedTuning_TemporNasalSpeedInc(sessionFileInfo, response, doSmooth)
                  

                        elseif contains(thisName, 'DirTuning')
                            fprintf('  Running Direction Tuning Analysis...\n');
                            
                            [~, sessionFileInfo] = get2PFramesByTrial(sessionFileInfo, thisName, true);
                            [~, sessionFileInfo] = getTrialGroups(sessionFileInfo, thisName);
                            [~, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, thisName, signalName);
                            % plotDirectionTuning(sessionFileInfo, response,tr, false)
                            disp('  Direction tuning processed.');
                        end

                        fprintf('  Completed: %s\n', thisName);
                    end % End of loop
                else
                    disp('RFMapping, DotFields, or dirTuning not found in otherVisualStim list.');
                end
            end

            % C. Process all VR Stimuli
            if ~isempty(vrStimNames)
                fprintf('Found %d VR stimulus file(s).\n', length(vrStimNames));

                % C.1 Preprocessing Looop (Generates individual response.mat files for ALL runs)
                for thisVRStim = 1:length(vrStimNames)
                    vrStimName = vrStimNames{thisVRStim};
                    fprintf('Processing VR Stim: %s\n', vrStimName);
                    try
                        % Preprocessing steps
                        
                        [~, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, vrStimName); % Reads and saves all tables; no additional computation done here
                        [~, sessionFileInfo] = findBonsaiPeripheralLag(sessionFileInfo, vrStimName, 1, interpRate);
                        [~, sessionFileInfo] = alignVRBonsaiToPeripheralData(sessionFileInfo,vrStimName);
                        [~, ~, ~, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo,interpRate,'TwoPFrameTime', vrStimName);
                        [response, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo, vrStimName);
                        [~, sessionFileInfo] = getVRTrialIndices(sessionFileInfo, vrStimName); %new
                        [~, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, vrStimName);
                        [~, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, vrStimName, zScoreProcessedSignals, applyTemporalSmoothing); %, prctlF, windowSize

                        % Calculates and saves LapPositionActivity (un-shuffled)
                        [~, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, vrStimName, useZScoredProcessedSignals);

                        [~, sessionFileInfo] = getRunningSpeedAcrossLaps(sessionFileInfo, vrStimName);
                        clear response;

                    catch
                        fprintf('WARNING: Could not preprocess VRStim name %s and will not be able to combine across VRRuns if appropriate..\n', vrStimName)
                    end
                end

                % Determine Path
                shouldProcessSeparately = any(ismember(vrStimNames, doNotCombine));


                % Path condition: TRUE if it's a single run OR if it's a multi-run exception.
                runAsIndividualLoop = (isscalar(vrStimNames)) || (length(vrStimNames) > 1 && shouldProcessSeparately);
                % Execute Shuffle and Analysis

                if runAsIndividualLoop
                    fprintf('Processing via INDIVIDUAL LOOP path.\n');

                    for thisVRStim = 1:length(vrStimNames)

                        vrStimName = vrStimNames{thisVRStim};
                        stimIdx = find(strcmp(vrStimName, {sessionFileInfo.stimFiles.name}), 1);
                        % load the individual response structure into the environment
                        response = load(sessionFileInfo.stimFiles(stimIdx).Response);
                        % Compute Shuffle Matrix (Stitch to itself)
                        fprintf(' -> Computing shuffle matrix for run: %s\n', vrStimName);
                        % change shuffle to include spikes 
                        %plotPopulationActivityAcrossConditions(sessionFileInfo, response, signalName, applyTemporalSmoothing)
                        %plotPopulationActivityAcrossConditions(sessionFileInfo, response, 'dFFNeuropilCorrected', applyTemporalSmoothing)
                        [response, sessionFileInfo] = computeShuffleMatrixForSession(sessionFileInfo,response,{vrStimName}, useZScoredProcessedSignals);
                       
                        %plotSortedPopulationResponse_OddEven(sessionFileInfo, response, signalName, true);

                        % plotAllNeuronSummariesToPDF_SparseNoiseAndVR(sessionFileInfo, response, applyTemporalSmoothing)
                        % plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing, signalName)                       % Run checks
                        % plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing, 'dFFNeuropilCorrected')
                        %plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing, 'dFF')

                        [~,~] = getRangeSignificance_fromShuffle(sessionFileInfo, response); % Run on final response
                        [~, ~] = getPeakSignificance_fromShuffle(sessionFileInfo, response); % Run on final response
                        [~,~,~]= computeVarianceAcrossPositionBins(sessionFileInfo, response); % Run on final response
                        [~,~,~] = checkOddEvenCorrelation(sessionFileInfo, response); % Run on final response
                        [~,~,~] = checkHalvesCorrelation(sessionFileInfo, response);

                        % For boutons specifically
                        if strcmpi(typeImaged, 'Boutons')
                            disp('Finding highly correlated boutons for this session')
                            [~,~,~,~] = findHightlyCorrelatedROIs(sessionFileInfo);
                        end
                    end

                else

                    % Combination Required (N>1 and NOT an exception) ---
                    fprintf('Processing via COMBINED PATH...This is to keep Sonali sane\n');
                    % Combine Response sturctures
                    fprintf('Checking for and combining Response across multiple VR runs...\n')

                    responseVRRuns = loadDataStructuresByKeyword(sessionFileInfo, responseName{thisMouse}, 'Response');
                    % TODOPRIORITY: check and append additional information from
                    % the ucl open vr into this
                    [response, sessionFileInfo] = combineResponseForVRRuns(responseVRRuns, sessionFileInfo);

                    % Compute Shuffle Matrix (Stitch ALL raw signals together)
                    fprintf('Computing shuffle matrix for COMBINED signals...\n')
               
                    %plotSortedPopulationResponse_OddEven(sessionFileInfo, response, signalName, true)
                    %plotPopulationActivityAcrossConditions(sessionFileInfo, response, signalName, applyTemporalSmoothing);
                    %plotPopulationActivityAcrossConditions(sessionFileInfo, response, 'dFFNeuropilCorrected', applyTemporalSmoothing);
                    % plotPopulationActivityAcrossConditions(sessionFileInfo, response, 'spks', applyTemporalSmoothing);
                    [response, sessionFileInfo] = computeShuffleMatrixForSession(sessionFileInfo, response, vrStimNames, useZScoredProcessedSignals);

                    % plotAllROIConditionSummaries(sessionFileInfo, response)
                    %plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing,'dFFNeuropilCorrected')
                    %plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing,signalName)
                    % plotAllNeuronSummariesToPDF_withThresholds(sessionFileInfo, response, true, signalName)
                    % Run checks
                    [~,~] = getRangeSignificance_fromShuffle(sessionFileInfo, response); % Run on final response
                    [~, ~] = getPeakSignificance_fromShuffle(sessionFileInfo, response); % Run on final response
                    [~,~,~]= computeVarianceAcrossPositionBins(sessionFileInfo, response); % Run on final response
                    [~,~,~] = checkOddEvenCorrelation(sessionFileInfo, response); % Run on final response
                    [~,~,~] = checkHalvesCorrelation(sessionFileInfo, response);
                    % plotPopulationActivityAcrossConditions_HalvesStableROIsOnly(sessionFileInfo, response, signalName, applyTemporalSmoothing);
                    % plotPopulationActivityAcrossConditions_HalvesStableROIsOnly(sessionFileInfo, response, 'dFFNeuropilCorrected', applyTemporalSmoothing);

                    % plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing, signalName);
                    % plotAllNeuronConditionsSummaries_VR_and_RF(sessionFileInfo, response, applyTemporalSmoothing,'dFFNeuropilCorrected');

                    % For boutons specifically
                    if strcmpi(typeImaged, 'Boutons')
                        disp('Finding highly correlated boutons for this session')
                        [~,~,~,~] = findHightlyCorrelatedROIs(sessionFileInfo); %TODO: flatten!!!
                    end

                end
            end

            fprintf('  Session %s processed successfully!\n', sessionName);

        catch ME

            % Still display the warning for immediate feedback
            warning('    Error processing %s, session %s: %s', mousenumber, sessionName, ME.message);
            fprintf(2, '    Error in function %s, line %d.\n', ME.stack(1).name, ME.stack(1).line);

            % Prepare the data for the log file
            timestamp = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
            mouseStr = string(mousenumber);
            sessionStr = string(sessionName);
            % Clean up the error message to remove newlines for better CSV formatting
            errorMessage = string(strrep(ME.message, newline, ' '));
            errorFunction = string(ME.stack(1).name);
            errorLine = ME.stack(1).line;

            % Create a table with the new error information
            errorData = table(timestamp, mouseStr, sessionStr, errorMessage, errorFunction, errorLine, ...
                'VariableNames', logHeaders);

            % Append the error data to the CSV file
            writetable(errorData, logFilePath, 'WriteMode', 'append', 'WriteVariableNames', false);
        end
    end
end
disp('All mice and sessions processed!');