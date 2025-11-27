function [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo, samplingRate, mainTimeToUse, VRStimName, plotFlag, trimNaNs)
%
%   Aligns and interpolates all VR experiment signals (Suite2P, Bonsai, peripheral)
%   to a unified 2P timebase using the Arduino clock.
%   This function supports flexible frame time alignment, bonsai lag correction, and visual
%   inspection of raw vs interpolated signals.
%
%   MODIFIED VERSION
%   Includes new 'trimNaNs' flag to remove NaN padding from the start/end
%   of all signals, ensuring a common, valid time window. This will only be
%   relavent for Soma recordings without z-motion correction. 
%
% Inputs:
%   sessionFileInfo : struct
%       Contains paths and metadata for the session (generated earlier in your pipeline).
%   samplingRate : (optional, default = 60)
%       Target sampling rate in Hz for all interpolated signals (e.g., 60 Hz).
%   mainTimeToUse : string (optional, default = 'TwoPFrameTime')
%       Timebase to align to; Use 'TwoPFrameTime' or 'ArduinoTime'
%   VRStimName : string
%       Name of the  visual stimulus file to process (e.g., 'RFStim_001').
%   plotFlag : logical (optional, default = true)
%       If true, generates sanity check plots.
%   trimNaNs : logical (optional, default = false) % 
%       If true, finds the common non-NaN time window across all signals
%       and trims all data to that window.
%
% Outputs:
%   processedTwoPData : struct
%       Contains resampled F, Fneu, spks, frame times, and ROI metadata.
%   bonsaiData : struct
%       Bonsai-tracked signals (e.g., mouse position, trial info, quadstate), all corrected for lag and resampled.
%   peripheralData : struct
%       Peripheral signals (e.g., photodiode, wheel), resampled to the same timebase.
%
% Example usage:
%   sessionFileInfo = getSessionInfo('Mouse01', 'Day1');
%   [processedTwoPData, bonsaiData, peripheralData, ~] = resamplAndAlignVR_BonsaiPeripheralSuite2P_trim(sessionFileInfo, 60, 'TwoPFrameTime', 'VRStim_001', true, true); % Trim NaNs
%
% Aman and Sonali - April 2025
%% Define default paramters and load appropriate data files
if nargin < 2, samplingRate = 60; end
if nargin < 3, mainTimeToUse = 'TwoPFrameTime'; end % This is the interrupt-arduino time from the Bonsai Arduino
if nargin < 5, plotFlag = true; end
if nargin < 6, trimNaNs = false; end % <<< NEW: Default for trimNaNs

% if nargin < 5, VR_number = 1; end
% Pick out VRCorr
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName,sessionFileInfo.stimFiles(iStim).name);
end
% stimList = find(bonsaiData.isVRstim==1); %%% have iStim (or) VR number as an input so you can run multiple VR files. And here just point to the first case if not specified
% iStim = stimList(VR_number);
iStim = find(bonsaiData.isVRstim==1);
% Load data files
if exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
    load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData')
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
else
    error('Merged Bonsai-Suite2P, BonsaiData and/or Peripheral data not found for this session.');
end
% Save output path for a new 2P data file
stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_processed2PData_' sessionFileInfo.stimFiles(iStim).name '.mat'];
sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
%% Create unified time base from the selected main time
% Concat the arduino time previously split into planes
raw2PTimes = vertcat(twoPData.(mainTimeToUse));
[unique2PTimes, mainTimeUniqueIdx] = unique(raw2PTimes);
% Range to interpolate all vectors
resample2PTimes = unique2PTimes(1):1/samplingRate:unique2PTimes(end);
sampleTimes = resample2PTimes;
% Define interpolation method for all and trialInfo
generalInterpMethod     = 'linear';
trialInfoInterpMethod   = 'nearest'; %%% is 'previous' more accurate
%% Interpolate: Two-photon time vectors
timeFields = {'TwoPFrameTime', 'ArduinoTime'}; % Excluding RenderFrameCount and LastSyncPulseTime
disp('Processing TwoP Frame Times')
% Arduino Times %%%% I don't think this is needed
for thisField = 1:numel(timeFields)
    fieldName = timeFields{thisField};
    concatenatedTimeVec = vertcat(twoPData.(fieldName));
    concatenatedTimeVec = concatenatedTimeVec(mainTimeUniqueIdx);
    processedTwoPData.(fieldName) = interp1(concatenatedTimeVec, concatenatedTimeVec, sampleTimes, generalInterpMethod)';
end
processedTwoPData.resample2PTimeUsed = mainTimeToUse; % For future use.
%% Interpolate: Suite2p data
disp('Processing Suite2P Data')
roiFields = {'F', 'Fneu', 'spks'};
%interpMethods = {'linear', 'linear', 'nearest'}; % Change if incorrect
for thisField = 1:numel(roiFields)
    processedTwoPData.(roiFields{thisField}) = [];
end
processedTwoPData.roiPlaneIdentity = [];
processedTwoPData.iscell = [];
processedTwoPData.redcell = [];
processedTwoPData.stat = {};
processedTwoPData.ops = {};
processedTwoPData.planeName = {};
for thisPlaneIdx = 1:numel(twoPData)
    % Using the chosen arduino 2p plane time
    rawArduinoPlaneTime = double(twoPData(thisPlaneIdx).(mainTimeToUse)); % why is double used here?
    for thisField = 1:numel(roiFields)
        fieldName = roiFields{thisField};
        signal = double(twoPData(thisPlaneIdx).(fieldName));
        interpolated = interp1(rawArduinoPlaneTime, signal', sampleTimes, generalInterpMethod)';
        processedTwoPData.(fieldName) = [processedTwoPData.(fieldName); interpolated];
    end
    processedTwoPData.iscell = [processedTwoPData.iscell; twoPData(thisPlaneIdx).iscell];
    processedTwoPData.redcell = [processedTwoPData.redcell; twoPData(thisPlaneIdx).redcell];
    processedTwoPData.roiPlaneIdentity = [processedTwoPData.roiPlaneIdentity; repmat(thisPlaneIdx-1, size(twoPData(thisPlaneIdx).F,1), 1)];
    processedTwoPData.stat = [processedTwoPData.stat, twoPData(thisPlaneIdx).stat];
    processedTwoPData.ops{end+1} = twoPData(thisPlaneIdx).ops;
    processedTwoPData.planeName{end+1} = twoPData(thisPlaneIdx).planeName;
end
%%  Interpolate: Peripheral - Wheel (no lag correction)
disp('Processing Peripheral Data: Wheel')
if isfield(peripheralData, 'Wheel')
    rawTime     = peripheralData.Wheel.rawArduinoTime;
    rawValue    = peripheralData.Wheel.rawValue;
    % @Aman Alternative naming for organising for interpolated vectors/
    % ResampledCorrArduinoTime? ResampledCorrValue?
    peripheralData.Wheel.Value      = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Wheel.sampleTimes = sampleTimes';
    %     peripheralData.Wheel.ArduinoTime= interp1(rawTime, rawTime, sampleTimes, generalInterpMethod, NaN)';
end
%%  Interpolate: Peripheral - Photodiode (no lag correction)
disp('Processing Peripheral Data: PD')
if isfield(peripheralData, 'Photodiode')
    rawTime = peripheralData.Photodiode.rawArduinoTime;
    rawValue = peripheralData.Photodiode.rawValue;
    peripheralData.Photodiode.Value = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Photodiode.sampleTimes = sampleTimes';
    %     peripheralData.Photodiode.ArduinoTime = interp1(rawTime, rawTime, sampleTimes, generalInterpMethod, NaN)';
end

%% Interpolate: Peripheral - Quadstate (lag corrected) / (changed to peripheral: 18.11.25)
% disp('Processing Bonsai Data: Quad')
% if isfield(bonsaiData, 'Quadstate')
%     rawValue = bonsaiData.Quadstate.rawValue;
%     rawTime = bonsaiData.Quadstate.rawArduinoTime;
%     % Lag corrected
%     lagCorrT = bonsaiData.Quadstate.rawCorrectedArduinoTime;
%     bonsaiData.Quadstate.Value = interp1(lagCorrT, rawValue, sampleTimes, generalInterpMethod, NaN)';
%     bonsaiData.Quadstate.sampleTimes = sampleTimes';
%     % Uncorrected
% %     bonsaiData.Quadstate.uncorrectedValue = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
% %     bonsaiData.Quadstate.uncorrectedArduinoTime = interp1(rawTime, rawTime, sampleTimes, generalInterpMethod, NaN)';
% end

disp('Processing Bonsai Data: Quad')
if isfield(peripheralData, 'Quadstate')
    rawValue = peripheralData.Quadstate.rawValue;
    rawTime = peripheralData.Quadstate.rawArduinoTime;
    % Lag corrected
    lagCorrT = peripheralData.Quadstate.rawCorrectedArduinoTime;
    peripheralData.Quadstate.Value = interp1(lagCorrT, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Quadstate.sampleTimes = sampleTimes';
    % Uncorrected
%     bonsaiData.Quadstate.uncorrectedValue = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
%     bonsaiData.Quadstate.uncorrectedArduinoTime = interp1(rawTime, rawTime, sampleTimes, generalInterpMethod, NaN)';
end
%% Interpolate: Bonsai - Mouse Position (lag corrected & uncorrected)
disp('Processing Bonsai Data: Mouse Position')
if isfield(bonsaiData, 'MousePos')
    rawValue = bonsaiData.MousePos.rawValue;
    rawTime = bonsaiData.MousePos.rawCorrectedArduinoTime;
    %     lagCorrT = bonsaiData.MousePos.rawCorrectedArduinoTime;
    % Lag corrected
    bonsaiData.MousePos.Value = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    bonsaiData.MousePos.sampleTimes = sampleTimes';
    %     bonsaiData.MousePos.ArduinoTime = interp1(lagCorrT, lagCorrT, sampleTimes, generalInterpMethod, NaN)';
    % Uncorrected %%%%% don't need this
%     bonsaiData.MousePos.uncorrectedValue = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
%     bonsaiData.MousePos.uncorrectedArduinoTime = interp1(rawTime, rawTime, sampleTimes, generalInterpMethod, NaN)';
end
%% Interpolate: Bonsai - TrialInfo @Aman - not sure if this is right
disp('Processing Bonsai Data: Trial Info')
if isfield(bonsaiData, 'TrialInfo')
    correctedStartTimeAll = bonsaiData.TrialInfo.rawCorrectedArduinoTime;
    uncorrectedStartTimeAll = bonsaiData.TrialInfo.rawArduinoTime;
    % Snap to nearest value in resampled 2P time vector
    bonsaiData.TrialInfo.StartTimeAll = interp1(sampleTimes, sampleTimes, correctedStartTimeAll, trialInfoInterpMethod);
    bonsaiData.TrialInfo.uncorrectedStartTimeAll = interp1(sampleTimes, sampleTimes, uncorrectedStartTimeAll, trialInfoInterpMethod);
end


%% NEW: Trim NaN padding if requested 
if trimNaNs
    disp('Trimming NaN padding from start/end of signals...');
    % Start with a full mask (sampleTimes is a row vector)
    combinedValidMask = true(size(sampleTimes));
    % 1. Check processedTwoPData (N_ROIs x N_Time)
    % A timepoint is invalid if *any* ROI is NaN at that time
    % We use all(isfinite(...)) to find columns with NO NaNs.
    combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.F), 1);
    combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.Fneu), 1);
    combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.spks), 1);

    % Check peripheralData (N_Time x 1 vectors)
    % Ensure they are row vectors for logical AND
    % if isfield(peripheralData, 'Wheel')
    %     combinedValidMask = combinedValidMask & isfinite(peripheralData.Wheel.Value)'; 
    % end
    % if isfield(peripheralData, 'Photodiode')
    %     combinedValidMask = combinedValidMask & isfinite(peripheralData.Photodiode.Value)';
    % end
    % 
    % if isfield(peripheralData, 'Quadstate')
    %     combinedValidMask = combinedValidMask & isfinite(peripheralData.Quadstate.Value)';
    % end
    % 
    % % Check bonsaiData (N_Time x 1 vectors)
    % if isfield(bonsaiData, 'MousePos')
    %     combinedValidMask = combinedValidMask & isfinite(bonsaiData.MousePos.Value)';
    % end


    
    % Check if we'd be removing everything
    if ~any(combinedValidMask)
        warning('Trimming NaNs would remove all data. Skipping trim.');
    else
        % Apply the mask to all time-series data
        processedTwoPData.trimmedNaNPadding = true; 
        processedTwoPData.trimmedMetaData = ['Trimming ' num2str(numel(sampleTimes)) ' samples down to ' num2str(sum(combinedValidMask)) ' samples.'];
        disp(['Trimming ' num2str(numel(sampleTimes)) ' samples down to ' num2str(sum(combinedValidMask)) ' samples.']);

        % Get the new, trimmed time vector (for filtering TrialInfo)
        newSampleTimes = sampleTimes(combinedValidMask);

        % Trim processedTwoPData
        processedTwoPData.F = processedTwoPData.F(:, combinedValidMask);
        processedTwoPData.Fneu = processedTwoPData.Fneu(:, combinedValidMask);
        processedTwoPData.spks = processedTwoPData.spks(:, combinedValidMask);
        processedTwoPData.TwoPFrameTime = processedTwoPData.TwoPFrameTime(combinedValidMask);
        processedTwoPData.ArduinoTime = processedTwoPData.ArduinoTime(combinedValidMask);

        % Trim peripheralData
        if isfield(peripheralData, 'Wheel')
            peripheralData.Wheel.Value = peripheralData.Wheel.Value(combinedValidMask);
            peripheralData.Wheel.sampleTimes = peripheralData.Wheel.sampleTimes(combinedValidMask);
        end
        if isfield(peripheralData, 'Photodiode')
            peripheralData.Photodiode.Value = peripheralData.Photodiode.Value(combinedValidMask);
            peripheralData.Photodiode.sampleTimes = peripheralData.Photodiode.sampleTimes(combinedValidMask);
        end
        if isfield(peripheralData, 'Quadstate')
            peripheralData.Quadstate.Value = peripheralData.Quadstate.Value(combinedValidMask);
            peripheralData.Quadstate.sampleTimes = peripheralData.Quadstate.sampleTimes(combinedValidMask);
        end

        % Trim bonsaiData
        if isfield(bonsaiData, 'MousePos')
            bonsaiData.MousePos.Value = bonsaiData.MousePos.Value(combinedValidMask);
            bonsaiData.MousePos.sampleTimes = bonsaiData.MousePos.sampleTimes(combinedValidMask);
        end

        % Filter TrialInfo data to be within the new time range
        if isfield(bonsaiData, 'TrialInfo')
            minTime = newSampleTimes(1);
            maxTime = newSampleTimes(end);

            % Find trials whose snapped start times fall within the valid range
            keepTrials = (bonsaiData.TrialInfo.StartTimeAll >= minTime) & ...
                         (bonsaiData.TrialInfo.StartTimeAll <= maxTime);

            bonsaiData.TrialInfo.StartTimeAll = bonsaiData.TrialInfo.StartTimeAll(keepTrials);

            % Also filter the uncorrected times using the same trial indices
            if isfield(bonsaiData.TrialInfo, 'uncorrectedStartTimeAll') && ...
               numel(bonsaiData.TrialInfo.uncorrectedStartTimeAll) == numel(keepTrials)
                bonsaiData.TrialInfo.uncorrectedStartTimeAll = bonsaiData.TrialInfo.uncorrectedStartTimeAll(keepTrials);
            end
            
            % NOTE: If TrialInfo contains other fields (e.g., trial outcomes)
            % that are indexed per-trial, they must also be filtered here
            % using the 'keepTrials' logical mask.
        end
    end
end
%% Sanity check plots
if plotFlag
    % 2P Times
    figure('Name', 'TwoP Arduino Frame Times');
    hold on;
    histogram(diff(raw2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Original');
    histogram(diff(unique2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Unique');
    histogram(diff(processedTwoPData.(mainTimeToUse)), 'BinWidth', 0.001, 'DisplayName', 'Resampled');
    xlabel('Time Diff (s)'); ylabel('Count');
    legend; title('2P Frame Time Distribution');
    xlim([0 0.2]);
    % Neuron trace comparison
    planeIndex = 1;
    fOrig = double(twoPData(planeIndex).F);
    originalTime = double(twoPData(planeIndex).(mainTimeToUse));
    nROIs = size(fOrig, 1);
    
    roiMask = processedTwoPData.roiPlaneIdentity == (planeIndex - 1);
    
    % Check if any ROIs remain for this plane after trimming
    if any(roiMask)
        fResampled = processedTwoPData.F(roiMask, :);
        
        % Ensure we don't try to plot more ROIs than exist
        if size(fResampled,1) < 5
             roiIndices = 1:size(fResampled,1);
        else
             roiIndices = randperm(size(fResampled,1), 5);
        end
        
        figure('Name', 'Calcium Trace Comparison');
        for idx = 1:numel(roiIndices)
            roi = roiIndices(idx);
            subplot(numel(roiIndices), 1, idx);
            hold on;
            % Note: fOrig uses 'roi', fResampled uses 'roi' as an index into the *resampled* matrix
            plot(originalTime, fOrig(roi, :), 'k','LineWidth', 1.2, 'DisplayName', 'Original');
            plot(processedTwoPData.(mainTimeToUse), fResampled(roi, :), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
            title(sprintf('ROI %d (Plane %d)', roi, planeIndex));
            ylabel('F');
            if idx == 1
                legend();
            end
        end
        xlabel('Time (s)');
        sgtitle('Neuron Trace Resampling: Original vs Interpolated');
    else
        disp('Skipping calcium trace plot: No ROIs found for plane 1 in processed data (possibly trimmed).');
    end

    figure('Name', 'Resampled Peripheral Signals');
    plotIdx = 1;
    if isfield(peripheralData, 'Photodiode')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Photodiode.sampleTimes, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('Photodiode: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end
    if isfield(peripheralData, 'Wheel')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(peripheralData.Wheel.rawArduinoTime, peripheralData.Wheel.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Wheel.sampleTimes, peripheralData.Wheel.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('Wheel: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end
    figure('Name', 'Resampled Bonsai Signals with Lag Correction');
    plotIdx = 1;
    if isfield(bonsaiData, 'MousePos')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.MousePos.sampleTimes, bonsaiData.MousePos.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('MousePos: Raw vs Lag-Corrected'); legend; ylabel('Position'); xlabel('Time');
    end
    if isfield(bonsaiData, 'Quadstate')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.Quadstate.rawArduinoTime, bonsaiData.Quadstate.rawValue, 'k', 'DisplayName', 'raw Time'); hold on;
        plot(bonsaiData.Quadstate.sampleTimes, bonsaiData.Quadstate.Value, 'r-', 'DisplayName', 'Corrected/Interpolated Time');
        title('Quadstate: Raw vs Lag-Corrected'); legend; ylabel('Quad'); xlabel('Time');
    end
end
%% Save
disp('Saving processed data files...');
save(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'processedTwoPData', '-v7.3');
disp('Saved processedTwoPData');
save(sessionFileInfo.stimFiles(iStim).BonsaiData, "bonsaiData");
disp('Saved bonsaiData');
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData");
disp('Saved peripheralData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Saved sessionFileInfo');
end