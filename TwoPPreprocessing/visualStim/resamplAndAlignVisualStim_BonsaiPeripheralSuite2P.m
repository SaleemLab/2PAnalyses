function [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVisualStim_BonsaiPeripheralSuite2P(sessionFileInfo, samplingRate, mainTimeToUse, StimName, plotFlag, trimNaNs, overwrite)
%   Aligns and interpolates all VR experiment signals (Suite2P, Bonsai, peripheral)
%   to a unified 2P timebase using the Arduino clock.
%
%   - Overwrites processedTwoPData file as a FLATTENED file (-struct).
%   - Appends/Updates bonsaiData and peripheralData as nested structs.

%% Define default parameters
if nargin < 2, samplingRate = 60; end
if nargin < 3, mainTimeToUse = 'TwoPFrameTime'; end 
if nargin < 5, plotFlag = false; end
if nargin < 6, trimNaNs = true; end
if nargin < 7, overwrite = true; end 

%% File Paths
bonsaiData = []; 
isStim = strcmp(StimName, {sessionFileInfo.stimFiles.name});
thisStim = find(isStim, 1);
if isempty(thisStim)
    error('Stimulus %s not found in sessionFileInfo', StimName);
end

% Define Input File Paths
fileMerged = sessionFileInfo.stimFiles(thisStim).mergedBonsai2PSuite2pData;
fileBonsai = sessionFileInfo.stimFiles(thisStim).BonsaiData;
filePeripheral = sessionFileInfo.stimFiles(thisStim).processedPeripheralData;

% Define Output File Path (processedTwoPData)
stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_processed2PData_' sessionFileInfo.stimFiles(thisStim).name '.mat'];
outputSavePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
sessionFileInfo.stimFiles(thisStim).processedMergedBonsaiSuite2pData = outputSavePath;

%% Overwrite Check (Updated for Flattened Loading)
if exist(outputSavePath, 'file') && ~overwrite
    disp(['Processed 2P Data for ' StimName ' already exists. Loading flattened file...']);
    
    % Reconstruct the struct from individual variables in the flat file
    processedTwoPData = load(outputSavePath);
    
    if exist(fileBonsai, 'file')
        load(fileBonsai, 'bonsaiData');
    else
        bonsaiData = [];
    end
    
    if exist(filePeripheral, 'file')
        load(filePeripheral, 'peripheralData');
    else
        peripheralData = [];
    end
    
    disp('Data loaded from file. Skipping processing.');
    return; 
end

%% Load raw data streams 
disp(['Processing ' StimName '...']);
if ~exist(fileMerged, 'file'), error('Merged data file not found: %s', fileMerged); end
if ~exist(filePeripheral, 'file'), error('Peripheral data file not found: %s', filePeripheral); end

load(fileMerged, 'twoPData');
load(filePeripheral, 'peripheralData');
if ~exist(fileBonsai, 'file')
    warning('Bonsai data file not found, continuing without it: %s', fileBonsai);
else
    load(fileBonsai, 'bonsaiData');
end

%% Create unified time base 
raw2PTimes = vertcat(twoPData.(mainTimeToUse));
[unique2PTimes, mainTimeUniqueIdx] = unique(raw2PTimes);
resample2PTimes = unique2PTimes(1):1/samplingRate:unique2PTimes(end);
sampleTimes = resample2PTimes;
generalInterpMethod     = 'linear';
trialInfoInterpMethod   = 'nearest'; 

%% Interpolate: Two-photon time vectors
timeFields = {'TwoPFrameTime', 'ArduinoTime'}; 
disp('Processing TwoP Frame Times')
for thisField = 1:numel(timeFields)
    fieldName = timeFields{thisField};
    concatenatedTimeVec = vertcat(twoPData.(fieldName));
    concatenatedTimeVec = concatenatedTimeVec(mainTimeUniqueIdx);
    processedTwoPData.(fieldName) = interp1(concatenatedTimeVec, concatenatedTimeVec, sampleTimes, generalInterpMethod)';
end
processedTwoPData.resample2PTimeUsed = mainTimeToUse; 

%% Interpolate: Suite2p data
disp('Processing Suite2P Data')
roiFields = {'F', 'Fneu', 'spks'};
for thisField = 1:numel(roiFields), processedTwoPData.(roiFields{thisField}) = []; end
processedTwoPData.roiPlaneIdentity = [];
processedTwoPData.iscell = [];
processedTwoPData.redcell = [];
processedTwoPData.stat = {};
processedTwoPData.ops = {};
processedTwoPData.planeName = {};

for thisPlaneIdx = 1:numel(twoPData)
    rawArduinoPlaneTime = double(twoPData(thisPlaneIdx).(mainTimeToUse)); 
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

%% Interpolate: Peripheral - Wheel
disp('Processing Peripheral Data: Wheel')
if isfield(peripheralData, 'Wheel')
    rawTime     = peripheralData.Wheel.rawArduinoTime;
    rawValue    = peripheralData.Wheel.rawValue;
    peripheralData.Wheel.Value      = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Wheel.sampleTimes = sampleTimes';
end

%% GrayScreen Check: Remove fields if necessary 
if contains(StimName, 'GrayScreen', 'IgnoreCase',true)
    if exist('peripheralData', 'var')
        peripheralfieldsToRemove = {'Quadstate', 'Photodiode'};
        existingFieldsPeripheral = intersect(fieldnames(peripheralData), peripheralfieldsToRemove);
        if ~isempty(existingFieldsPeripheral)
            peripheralData = rmfield(peripheralData, existingFieldsPeripheral);
            fprintf('GrayScreen detected. In peripheralData, removed fields: %s\n', strjoin(existingFieldsPeripheral, ', '));
        end
    end
    
    if ~isempty(bonsaiData)
        bonsaiFieldsToRemove = {'MousePos'};
        existingFieldsBonsai = intersect(fieldnames(bonsaiData), bonsaiFieldsToRemove);
        if ~isempty(existingFieldsBonsai)
            bonsaiData = rmfield(bonsaiData, existingFieldsBonsai);
            fprintf('GrayScreen detected. In bonsaiData, removed fields: %s\n', strjoin(existingFieldsBonsai, ', '));
        end
    end 
end

%% Interpolate: Peripheral - Photodiode
disp('Processing Peripheral Data: PD')
if isfield(peripheralData, 'Photodiode')
    rawTime = peripheralData.Photodiode.rawArduinoTime;
    rawValue = peripheralData.Photodiode.rawValue;
    peripheralData.Photodiode.Value = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Photodiode.sampleTimes = sampleTimes';
end

%% Interpolate: Peripheral - Quadstate 
disp('Processing Peripheral Data: Quad')
if isfield(peripheralData, 'Quadstate')
    rawValue = peripheralData.Quadstate.rawValue;
    rawTime = peripheralData.Quadstate.rawArduinoTime;
    [uniqueRawTime, idx] = unique(rawTime);
    uniqueRawValue = double(rawValue(idx));
    
    Value = interp1(uniqueRawTime, uniqueRawValue, sampleTimes, 'previous', NaN)'; 
    peripheralData.Quadstate.Value = Value;
    peripheralData.Quadstate.sampleTimes = sampleTimes';
end

%% Interpolate: Bonsai - TrialInfo
disp('Processing Bonsai Data: StimOnset Info')
fieldsToProcess = { 'ArduinoTimeRaw', 'ArduinoTime' };
for thisField = 1:size(fieldsToProcess, 1)
    rawField = fieldsToProcess{thisField, 1};
    targetField = fieldsToProcess{thisField, 2};
    if isfield(bonsaiData, rawField)
        rawTime = bonsaiData.(rawField);
        processedValues = interp1(sampleTimes, sampleTimes, rawTime, trialInfoInterpMethod); 
        bonsaiData.(targetField) = processedValues;
        fprintf('Processed %s into %s\n', rawField, targetField);
    end
end

%% Trim NaN padding
if trimNaNs
    disp('Trimming NaN padding...');
    combinedValidMask = true(size(sampleTimes));
    
    combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.F), 1);
    combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.Fneu), 1);
    combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.spks), 1);
    
    if ~any(combinedValidMask)
        warning('Trimming NaNs would remove all data. Skipping trim.');
    else
        processedTwoPData.trimmedNaNPadding = true;
        newSampleTimes = sampleTimes(combinedValidMask);
        
        processedTwoPData.F = processedTwoPData.F(:, combinedValidMask);
        processedTwoPData.Fneu = processedTwoPData.Fneu(:, combinedValidMask);
        processedTwoPData.spks = processedTwoPData.spks(:, combinedValidMask);
        processedTwoPData.TwoPFrameTime = processedTwoPData.TwoPFrameTime(combinedValidMask);
        processedTwoPData.ArduinoTime = processedTwoPData.ArduinoTime(combinedValidMask);
        
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
        
        % Filter TrialInfo data
        if isfield(bonsaiData, 'ArduinoTime') && ~isempty(bonsaiData.ArduinoTime)
            referenceTime = bonsaiData.ArduinoTime;
            minTime = newSampleTimes(1);
            maxTime = newSampleTimes(end);
            keepTrials = (referenceTime >= minTime) & (referenceTime <= maxTime);
            nTotalTrials = numel(referenceTime);
            
            infoFields = fieldnames(bonsaiData);
            for f = 1:numel(infoFields)
                if numel(bonsaiData.(infoFields{f})) == nTotalTrials
                    bonsaiData.(infoFields{f}) = bonsaiData.(infoFields{f})(keepTrials);
                end
            end
        end
    end
end

%% Sanity check plots
if plotFlag
    % 2P Times
    figure('Name', ['TwoP Arduino Frame Times: ' StimName]);
    hold on;
    histogram(diff(raw2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Original');
    histogram(diff(unique2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Unique');
    histogram(diff(processedTwoPData.(mainTimeToUse)), 'BinWidth', 0.001, 'DisplayName', 'Resampled');
    xlabel('Time Diff (s)'); ylabel('Count'); legend; title('2P Frame Time Distribution'); xlim([0 0.2]);
    
    % Trace Comparison
    planeIndex = 1;
    roiMask = processedTwoPData.roiPlaneIdentity == (planeIndex - 1);
    if any(roiMask)
        fResampled = processedTwoPData.F(roiMask, :);
        fOrig = double(twoPData(planeIndex).F);
        originalTime = double(twoPData(planeIndex).(mainTimeToUse));
        if size(fResampled,1) < 5, roiIndices = 1:size(fResampled,1); else, roiIndices = randperm(size(fResampled,1), 5); end
        figure('Name', 'Calcium Trace Comparison');
        for idx = 1:numel(roiIndices)
            roi = roiIndices(idx);
            subplot(numel(roiIndices), 1, idx); hold on;
            plot(originalTime, fOrig(roi, :), 'k','LineWidth', 1.2, 'DisplayName', 'Original');
            plot(processedTwoPData.(mainTimeToUse), fResampled(roi, :), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
            ylabel('F'); if idx==1, legend; end
        end
        sgtitle('Neuron Trace Resampling');
    end
    
    % Peripheral Plots
    figure('Name', ['Peripheral Signals: ' StimName]); t = tiledlayout('flow');
    if isfield(peripheralData, 'Photodiode')
        nexttile;
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Photodiode.sampleTimes, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Resampled');
        title('Photodiode'); legend;
    end
    if isfield(peripheralData, 'Wheel')
        nexttile;
        plot(peripheralData.Wheel.rawArduinoTime, peripheralData.Wheel.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Wheel.sampleTimes, peripheralData.Wheel.Value, 'r-', 'DisplayName', 'Resampled');
        title('Wheel'); legend;
    end
    if isfield(peripheralData, 'Quadstate')
        nexttile;
        plot(peripheralData.Quadstate.rawArduinoTime, peripheralData.Quadstate.rawValue, 'k', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Quadstate.sampleTimes, peripheralData.Quadstate.Value, 'r-', 'DisplayName', 'Resampled');
        title('Quadstate'); legend;
    end
end

%% Save Section (Flattened processedTwoPData only)
disp('Saving processed data files...');

% Save processedTwoPData: Unpack struct into individual variables
save(outputSavePath, '-struct', 'processedTwoPData', '-v7.3');
disp(['Saved Flattened processedTwoPData to: ' outputSavePath]);

% bonsaiData: Standard append as a struct
if exist(fileBonsai, 'file')
    save(fileBonsai, 'bonsaiData', '-append');
    disp('Updated bonsaiData (Appended Struct)');
end

% peripheralData: Standard append as a struct
if exist(filePeripheral, 'file')
    save(filePeripheral, 'peripheralData', '-append');
    disp('Updated peripheralData (Appended Struct)');
end

% Update session info
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end