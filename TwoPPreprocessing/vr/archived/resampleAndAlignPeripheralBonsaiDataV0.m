function [peripheralData, bonsaiData] = resampleAndAlignPeripheralBonsaiData(sessionFileInfo, correctLag, plotFlag)
% Resample and align peripheral and bonsai data to a shared (2P frame) timebase
%
% This function interpolates peripheral (e.g., wheel, photodiode, quad) and 
% Bonsai (e.g., mouse position, trial timestamps) signals to a common resampled 
% time vector, previously defined during TwoP preprocessing. Optionally applies 
% a lag correction based on Arduino–Bonsai alignment.
%
% Inputs:
%   - sessionFileInfo (struct):
%       Structure containing file paths and metadata for the session, including
%       locations of Bonsai, Suite2P, and processed TwoP files.
%
%   - correctLag (logical, optional):
%       Whether to apply lag correction to timing-sensitive signals (default: true).
%       Applies only to: Quadstate, Wheel, MousePos, and TrialInfo.
%
%   - plotFlag (logical, optional):
%       Whether to plot sanity-check comparisons of raw vs resampled signals
%       (default: true).
%
% Outputs:
%   - peripheralData (struct):
%       Updated peripheral signals structure with the following fields:
%         - corrected/uncorrectedArduinoTime: interpolated timebase (based on lag)
%         - Value: resampled signal values
% 
%   - bonsaiData (struct):
%       Updated Bonsai signals (TrialInfo and MousePos) aligned to the same resampled
%       base and:  
%         - corrected/uncorrectedArduinoTime: interpolated timebase (based on lag)
%
% Notes (to recall after PIPs): 
%   - The shared resampledTime is extracted from peripheralData.resamplingInfo.
%   - Interpolation is done with linear method. This will however change
%       for the other time vectors that will be interpolated in the future. 
%   - resampledTime is derived from uniqueChosen2PTime(1):1/samplingRate:uniqueChosen2PTime(end);
%   - Lag correction is applied as:
%         correctedTime = rawArduinoTime - (xcorrBestLag / samplingRate)
%   - Wheel data is saved with both corrected and uncorrected time for flexibility.
%   - TrialInfo only contains discrete event times, so is aligned by nearest frame.
%   - Other time vectors such as BonsaiTime, RenderFrameCount and
%       LastSyncPuleTime missing here. @Aman 
%       
%
% Usage Example:
%   [peripheralData, bonsaiData] = resampleAndAlignPeripheralBonsaiData(sessionFileInfo, true, true);
%
%
% Sonali and Aman - March 2025


%% Assign default parameters 
if nargin < 2
    correctLag = true;
end 
if nargin < 3
    plotFlag = true;
end     

%% Load data files 
for iStim = 1:length(sessionFileInfo.stimFiles) % Find VR Stimulus
    bonsaiData.isVRstim(iStim) = strcmp('VRCorr', sessionFileInfo.stimFiles(iStim).name);
end
iStim = find(bonsaiData.isVRstim==1); 

if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');     
else
    error('PeriphralData and/or BonsaiData missing for VR recording');
end

%% Extract previously defined variables
TwoPTime = peripheralData.resampledTimes.TwoPTime;
samplingRate  = peripheralData.resampledTimes.samplingRate;

lagShift = 0;
if correctLag && isfield(bonsaiData, 'LagInfo') && isfield(bonsaiData.LagInfo, 'xcorrBestLag')
    lagShift = bonsaiData.LagInfo.xcorrBestLag*(1 / samplingRate);
end

%% Signals to align and interpolate 
allSignalNames = {'Photodiode', 'Wheel', 'Quadstate', 'MousePos'};
FieldsToResample = {'rawArduinoTime', 'rawValue'}; 

for signalIdx = 1:numel(allSignalNames)
    signalName = allSignalNames{signalIdx};

    if strcmp(signalName, 'MousePos')
        if ~isfield(bonsaiData, signalName), warning('%s not found in bonsaiData. Skipping.', signalName); continue; end
        data = bonsaiData.(signalName);
        targetStruct = 'bonsai';
    else
        if ~isfield(peripheralData, signalName), warning('%s not found in peripheralData. Skipping.', signalName); continue; end
        data = peripheralData.(signalName);
        targetStruct = 'peripheral';
    end

    for fieldIdx = 1:numel(FieldsToResample)
        rawFieldName = FieldsToResample{fieldIdx};
        if ~isfield(data, rawFieldName), continue; end

        isTimeField = contains(rawFieldName, 'Time');
        rawData = data.(rawFieldName);

        % Determine output name based on type and correction
        if isTimeField
            fieldBase = erase(rawFieldName, 'raw');
            if strcmp(rawFieldName, 'rawArduinoTime') && any(strcmp(signalName, {'Quadstate', 'Wheel', 'MousePos'}))
                outputFieldName = ['corrected' fieldBase];
            elseif strcmp(rawFieldName, 'rawArduinoTime')
                outputFieldName = ['uncorrected' fieldBase];
            else
                outputFieldName = fieldBase;
            end
        else
            outputFieldName = erase(rawFieldName, 'raw');
        end

        % Time interpolation
        if isTimeField
            if strcmp(rawFieldName, 'rawArduinoTime')
                [uniqueTime, ~] = unique(rawData);
                correctedTime = uniqueTime;
                if any(strcmp(signalName, {'Quadstate', 'Wheel', 'MousePos'}))
                    correctedTime = uniqueTime - lagShift;
                end

                if strcmp(signalName, 'Wheel')
                    ArduinoTime_uncorrected = interp1(uniqueTime, uniqueTime, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN)';
                    peripheralData.(signalName).uncorrectedArduinoTime = ArduinoTime_uncorrected;
                end

                interpolated = interp1(uniqueTime, correctedTime, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN)';
            else
                [uniqueTime, ~] = unique(rawData);
                interpTime = uniqueTime(1):1/samplingRate:uniqueTime(end);
                interpolated = interp1(uniqueTime, uniqueTime, interpTime, 'linear', NaN)';
            end

            if strcmp(targetStruct, 'bonsai')
                bonsaiData.(signalName).(outputFieldName) = interpolated;
            else
                peripheralData.(signalName).(outputFieldName) = interpolated;
            end

        % Value interpolation
        elseif strcmp(rawFieldName, 'rawValue')
            if isfield(data, 'rawArduinoTime')
                timeRef = data.rawArduinoTime;
                [uniqueTime, uniqueIdx] = unique(timeRef);

                if any(strcmp(signalName, {'Quadstate', 'Wheel', 'MousePos'}))
                    uniqueTime = uniqueTime - lagShift;
                end

                rawData = rawData(uniqueIdx);
                interpolated = interp1(uniqueTime, rawData, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN)';

                if strcmp(targetStruct, 'bonsai')
                    bonsaiData.(signalName).(outputFieldName) = interpolated;
                else
                    peripheralData.(signalName).(outputFieldName) = interpolated;
                end
            end
        end
    end
end




%% TrialInfo (timestamp only, aligned via MousePos.correctedArduinoTime)
rawMousePositionTime = bonsaiData.MousePos.rawArduinoTime;

if (~isempty(bonsaiData.TrialInfo))
    disp('Processing lap information..');
    
    % Extract the start time of the lap from the ArduinoTime column;
    rawStartTimeAll = bonsaiData.TrialInfo.rawArduinoTime;
  
    % Find the index of the start time in the mouse_position_time vector; 
    for n = 1:length(rawStartTimeAll)
        % Check if identical lap start time is found in the mouse position
        % time vector 
        startInd = find(rawStartTimeAll(n) == rawMousePositionTime);
        % Find the index in quadstate_time that matches the current start time
 
        % Handle rare cases where the exact start time isn't found in mouse position time
        if isempty(startInd)
            for count = 1:10
                startInd = find(rawStartTimeAll(n)+count == rawMousePositionTime);
                if ~isempty(startInd)
                    break;
                end
            end
        end
        
        % If the lap starts after the last recorded mouse position time,
        % break the process; Incomplete lab before recorded ended. 
        if rawStartTimeAll(n) >= rawMousePositionTime(end)
            break;  % Exit the loop as there's no valid data left
        end

        % Peripherals.corrected_sglxTime()
        bonsaiData.TrialInfo.originalStartTimeAll(n, 1) = rawStartTimeAll(startInd);  % Used mouse position time instead of peripherals.Time
        % Masa's code aligns based on reward time; here interpolating based
        % corrected and interpolated mouseposition time. 
        bonsaiData.TrialInfo.newStartTimeAll(n, 1) = interp1(TwoPTime, TwoPTime, bonsaiData.TrialInfo.originalStartTimeAll(n,1), 'nearest');
%        
        
    end




%% Sanity check plots 
if plotFlag
    figure('Name','Sanity Check: Raw vs Resampled');
    plotIdx = 1;
    exampleSignals = {'Photodiode', 'MousePos', 'Quadstate'};

    for s = 1:numel(exampleSignals)
        signal = exampleSignals{s};

        % Dynamically assign source struct
        if isfield(peripheralData, signal)
            dataSource = peripheralData;
        elseif isfield(bonsaiData, signal)
            dataSource = bonsaiData;
        else
            warning('%s not found in either data struct. Skipping.', signal);
            continue;
        end

        if isfield(dataSource.(signal), 'rawArduinoTime') && ...
           isfield(dataSource.(signal), 'rawValue') && ...
           isfield(dataSource.(signal), 'Value')

            rawT = dataSource.(signal).rawArduinoTime;
            rawV = dataSource.(signal).rawValue;

            if isfield(dataSource.(signal), 'correctedArduinoTime')
                interpT = dataSource.(signal).correctedArduinoTime;
            elseif isfield(dataSource.(signal), 'uncorrectedArduinoTime')
                interpT = dataSource.(signal).uncorrectedArduinoTime;
            else
                interpT = TwoPTime;
            end

            interpV = dataSource.(signal).Value;

            subplot(numel(exampleSignals), 1, plotIdx);
            hold on;
            plot(rawT, rawV, 'k.', 'DisplayName', 'Raw');
            plot(interpT, interpV, 'r-', 'DisplayName', 'Interpolated');
            title(sprintf('%s: Raw vs Resampled', signal));
            xlabel('Time (s)');
            ylabel('Signal');
            legend;
            plotIdx = plotIdx + 1;
        end
    end
end

%% Saving 
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData")
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
