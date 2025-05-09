function [ES, bonsaiData, ES] = resamplAndAlign_BonsaiPeripheralSuite2P_VR(sessionFileInfo, samplingRate, mainTimeToUse, plotFlag)
%
%   Aligns and interpolates all VR experiment signals (Suite2P, Bonsai, peripheral)
%   to a unified 2P timebase using the Arduino clock.
%   This function supports flexible frame time alignment, bondai lag correction, and visual
%   inspection of raw vs interpolated signals.
%
% Inputs:
%   sessionFileInfo : struct
%       Contains paths and metadata for the session (generated earlier in your pipeline).
%   samplingRate : (optional, default = 60)
%       Target sampling rate in Hz for all interpolated signals (e.g., 60 Hz).
%   mainTimeToUse : string (optional, default = 'TwoPFrameTime')
%       Timebase to align to; Use 'TwoPFrameTime' or 'ArduinoTime'
%   plotFlag : logical (optional, default = true)
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
%   [processedTwoPData, bonsaiData, peripheralData] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo, 60, 'TwoPFrameTime', true);
%
% Aman and Sonali - April 2025
%% Define default paramters and load appropriate data files
if nargin < 2, samplingRate = 60; end
if nargin < 3, mainTimeToUse = 'TwoPFrameTime'; end % This is the interrupt-arduino time from the Bonsai Arduino
if nargin < 4, plotFlag = true; end
if nargin < 5, VR_number = 1; end

% Pick out VRCorr
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp('VRCorr',sessionFileInfo.stimFiles(iStim).name);
end

stimList = find(bonsaiData.isVRstim==1); %%% have iStim (or) VR number as an input so you can run multiple VR files. And here just point to the first case if not specified
iStim = stimList(VR_number);

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
stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_ES_' sessionFileInfo.stimFiles(iStim).name '.mat'];
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
    ES.(fieldName) = interp1(concatenatedTimeVec, concatenatedTimeVec, sampleTimes, generalInterpMethod)';
end

ES.resample2PTimeUsed = mainTimeToUse; % For future use.

%% Interpolate: Suite2p data
disp('Processing Suite2P Data')
roiFields = {'F', 'Fneu', 'spks'};
%interpMethods = {'linear', 'linear', 'nearest'}; % Change if incorrect

for thisField = 1:numel(roiFields)
    ES.(roiFields{thisField}) = [];
end

ES.roiPlaneIdentity = [];
ES.iscell = [];
ES.redcell = [];
ES.stat = {};
ES.ops = {};
ES.planeName = {};

for thisPlaneIdx = 1:numel(twoPData)
    % Using the chosen arduino 2p plane time
    rawArduinoPlaneTime = double(twoPData(thisPlaneIdx).(mainTimeToUse)); % why is double used here?

    for thisField = 1:numel(roiFields)
        fieldName = roiFields{thisField};
        signal = double(twoPData(thisPlaneIdx).(fieldName));
        interpolated = interp1(rawArduinoPlaneTime, signal', sampleTimes, generalInterpMethod)';
        ES.(fieldName) = [ES.(fieldName); interpolated];
    end

    ES.iscell = [ES.iscell; twoPData(thisPlaneIdx).iscell];
    ES.redcell = [ES.redcell; twoPData(thisPlaneIdx).redcell];
    ES.roiPlaneIdentity = [ES.roiPlaneIdentity; repmat(thisPlaneIdx-1, size(twoPData(thisPlaneIdx).F,1), 1)];
    ES.stat = [ES.stat, twoPData(thisPlaneIdx).stat];
    ES.ops{end+1} = twoPData(thisPlaneIdx).ops;
    ES.planeName{end+1} = twoPData(thisPlaneIdx).planeName;
end

%%  Interpolate: Peripheral - Wheel (no lag correction)
disp('Processing Peripheral Data: Wheel')
if isfield(peripheralData, 'Wheel')
    rawTime     = peripheralData.Wheel.rawArduinoTime;
    rawValue    = peripheralData.Wheel.rawValue;
    % @Aman Alternative naming for organising for interpolated vectors/
    % ResampledCorrArduinoTime? ResampledCorrValue?
    ES.wheelValue      = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    ES.wheelSampleTimes = sampleTimes;
    %     peripheralData.Wheel.ArduinoTime= interp1(rawTime, rawTime, sampleTimes, generalInterpMethod, NaN)';
end

%%  Interpolate: Peripheral - Photodiode (no lag correction)
disp('Processing Peripheral Data: PD')
if isfield(peripheralData, 'Photodiode')
    rawTime = peripheralData.Photodiode.rawArduinoTime;
    rawValue = peripheralData.Photodiode.rawValue;
    ES.photodiodeValue = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    ES.photodiodeSampleTimes = sampleTimes';

end

%% Interpolate: Bonsai - Mouse Position (lag corrected & uncorrected)
disp('Processing Bonsai Data: Mouse Position')
if isfield(bonsaiData, 'MousePos')
    rawValue = bonsaiData.MousePos.rawValue;
    rawCorrectedTime = bonsaiData.MousePos.rawCorrectedArduinoTime;
    % Lag corrected
     ES.mousePosValue = interp1(rawCorrectedTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
     ES.mousePosSampleTimes = sampleTimes';
end


%% Interpolate: Bonsai - TrialInfo @Aman - not sure if this is right
disp('Processing Bonsai Data: Trial Info')
if isfield(bonsaiData, 'TrialInfo')
    correctedStartTimeAll = bonsaiData.TrialInfo.rawCorrectedArduinoTime;
    uncorrectedStartTimeAll = bonsaiData.TrialInfo.rawArduinoTime;
    % Snap to nearest value in resampled 2P time vector % Try previous 
    ES.StartTimeAll = interp1(sampleTimes, sampleTimes, correctedStartTimeAll, trialInfoInterpMethod);
    ES.uncorrectedStartTimeAll = interp1(sampleTimes, sampleTimes, uncorrectedStartTimeAll, trialInfoInterpMethod);
end

%% Interpolate: Bonsai - Quadstate (lag corrected)
disp('Processing Bonsai Data: Quad')
if isfield(bonsaiData, 'Quadstate')
    rawValue = bonsaiData.Quadstate.rawValue;
    rawTime = bonsaiData.Quadstate.rawArduinoTime;
    % Lag corrected
    rawCorrectedTime = bonsaiData.Quadstate.rawCorrectedArduinoTime;
    ES.quadstateValue = interp1(rawCorrectedTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    ES.quadstateSampleTime = sampleTimes; 
    % Uncorrected
    ES.quadstateUncorrectedValue = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';

end

%% Sanity check plots
if plotFlag
    % 2P Times
    figure('Name', 'TwoP Arduino Frame Times');
    hold on;
    histogram(diff(raw2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Original');
    histogram(diff(unique2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Unique');
    histogram(diff(ES.(mainTimeToUse)), 'BinWidth', 0.001, 'DisplayName', 'Resampled');
    xlabel('Time Diff (s)'); ylabel('Count');
    legend; title('2P Frame Time Distribution');
    xlim([0 0.2]);

    % Neuron trace comparison
    planeIn dex = 1;
    fOrig = double(twoPData(planeIndex).F);
    originalTime = double(twoPData(planeIndex).(mainTimeToUse));
    nROIs = size(fOrig, 1);
    roiIndices = randperm(nROIs, 5);
    roiMask = ES.roiPlaneIdentity == (planeIndex - 1);
    fResampled = ES.F(roiMask, :);

    figure('Name', 'Calcium Trace Comparison');
    for idx = 1:numel(roiIndices)
        roi = roiIndices(idx);
        subplot(numel(roiIndices), 1, idx);
        hold on;
        plot(originalTime, fOrig(roi, :), 'k','LineWidth', 1.2, 'DisplayName', 'Original');
        plot(ES.(mainTimeToUse), fResampled(roi, :), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
        title(sprintf('ROI %d (Plane %d)', roi, planeIndex));
        ylabel('F');
        if idx == 1
            legend();
        end
    end
    xlabel('Time (s)');
    sgtitle('Neuron Trace Resampling: Original vs Interpolated');

    figure('Name', 'Resampled Peripheral Signals');
    plotIdx = 1;
    if isfield(ES, 'Photodiode')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(ES.photodiodeSampleTimes, ES.photodiodeValue, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('Photodiode: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end

    if isfield(ES, 'Wheel')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(ES.Wheel.rawArduinoTime, ES.Wheel.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(ES.Wheel.ArduinoTime, ES.Wheel.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('Wheel: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end

    figure('Name', 'Resampled Bonsai Signals with Lag Correction');
    plotIdx = 1;
    if isfield(bonsaiData, 'MousePos')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.MousePos.ArduinoTime, bonsaiData.MousePos.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('MousePos: Raw vs Lag-Corrected'); legend; ylabel('Position'); xlabel('Time');
    end

    if isfield(bonsaiData, 'Quadstate')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.Quadstate.rawArduinoTime, bonsaiData.Quadstate.rawValue, 'k', 'DisplayName', 'raw Time'); hold on;
        plot(bonsaiData.Quadstate.ArduinoTime, bonsaiData.Quadstate.Value, 'r-', 'DisplayName', 'Corrected/Interpolated Time');
        title('Quadstate: Raw vs Lag-Corrected'); legend; ylabel('Quad'); xlabel('Time');
    end
end

%% Save
disp('Saving processed data files...');
save(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'ES', '-v7.3');
disp('Saved processedTwoPData');
save(sessionFileInfo.stimFiles(iStim).BonsaiData, "bonsaiData");
disp('Saved bonsaiData');
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "ES");
disp('Saved peripheralData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Saved sessionFileInfo');
end