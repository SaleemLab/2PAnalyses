function [peripheralData, bonsaiData] = resampleAndAlignPeripheralBonsaiDataV2(sessionFileInfo, correctLag, plotFlag)
% TODO: PD on and offs missing; possibly not important but save anyway.
% Interpolate all other 'time' vectors saved for each Bonsai and
% Checking times 
% Peripheral field - currently only focused on Arduino Times. 
%
% Resample and align peripheral and Bonsai data to a shared 2P timebase
%
% This function interpolates peripheral (e.g., Wheel, Photodiode) and Bonsai
% (e.g., Mouse Position, Quadstate, TrialInfo) signals to a common TwoP timebase
% defined during Suite2P processing. Lag correction is applied only to Bonsai signals
% when specified. Trial start timestamps are mapped to the closest frame using the
% corrected TwoP time, and corresponding 2P frame indices are also saved.
%
% Inputs:
%   - sessionFileInfo (struct): contains session file paths and metadata.
%   - correctLag (logical, optional): apply lag correction to Bonsai signals (default: true).
%   - plotFlag (logical, optional): plot sanity checks of resampled signals (default: true).
%
% Outputs: (CHANGE)
%   - peripheralData (struct):
%       - .TwoPTime: reference time vector used for interpolation (uncorrected).
%       - .Wheel / .Photodiode:
%           - .Value: interpolated signal values.
%           - .ArduinoTime: interpolated time vector (same as TwoPTime).
% 
%   - bonsaiData (struct):
%       - .correctedTwoPTime: lag-corrected time vector used for interpolation.
%       - .MousePos / .Quadstate:
%           - .Value: lag-corrected interpolated values.
%           - .ArduinoTime: interpolated lag-corrected time vector.
%           - .correctedArduinoTime: lag-corrected raw time vector.
%           - .uncorrectedArduinoTime (Quadstate only): original raw time vector.
%       - .TrialInfo:
%           - .originalStartTimeAll: mapped start times from raw ArduinoTime to raw TwoPTime.
%           - .alignedStartTimeAll: snapped to nearest frame in correctedTwoPTime.
%           - .alignedStartFrameAll: index of the nearest frame in correctedTwoPTime.
%
% Notes:
%   - Bonsai signals are aligned using: correctedTime = rawTime - (xcorrBestLag / samplingRate)
%   - Signals are interpolated using: TwoPTime(1):1/samplingRate:TwoPTime(end)
%   - TrialInfo timestamps are mapped to nearest corrected frame using correctedTwoPTime.
%
% Usage:
%   [peripheralData, bonsaiData] = resampleAndAlignPeripheralBonsaiData(sessionFileInfo, true, true);
%
% Sonali and Aman - March 2025

if nargin < 2, correctLag = false; end
if nargin < 3, plotFlag = true; end

for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp('VRCorr', sessionFileInfo.stimFiles(iStim).name);
end

iStim = find(bonsaiData.isVRstim == 1);

if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
   exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
else
    error('Missing Bonsai or peripheral data file');
end

% Resampled two photon frame time to use for interpolation 
TwoPTime = peripheralData.resampledTimes.TwoPTime;
samplingRate = peripheralData.resampledTimes.samplingRate;
peripheralData.TwoPTime = TwoPTime;

% Lag shift for Bonsai correction
lagShift = 0;
if correctLag && isfield(bonsaiData, 'LagInfo') && isfield(bonsaiData.LagInfo, 'xcorrBestLag')
    lagShift = bonsaiData.LagInfo.xcorrBestLag * (1 / samplingRate);
    % Lag corrected & resampled 
    
    bonsai2PTime = bonsaiData.resampledTimes.TwoPTime - lagShift;
    bonsaiData.correctedTwoPTime = bonsai2PTime;
else
    % Resampled bonsai time 
    bonsai2PTime = bonsaiData.resampledTimes.TwoPTime;
    bonsaiData.TwoPTime = bonsai2PTime;
end



%% Peripheral: Wheel (no lag correction)
if isfield(peripheralData, 'Wheel')
    rawT = peripheralData.Wheel.rawArduinoTime;
    rawV = peripheralData.Wheel.rawValue;
    peripheralData.Wheel.Value = interp1(rawT, rawV, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
    peripheralData.Wheel.ArduinoTime = interp1(rawT, rawT, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
end

%% Peripheral: Photodiode (no lag correction)
if isfield(peripheralData, 'Photodiode')
    rawT = peripheralData.Photodiode.rawArduinoTime;
    rawV = peripheralData.Photodiode.rawValue;
    peripheralData.Photodiode.Value = interp1(rawT, rawV, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
    peripheralData.Photodiode.ArduinoTime = interp1(rawT, rawT, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
end

%% Mouse Position (lag corrected)
if isfield(bonsaiData, 'MousePos')
    rawV = bonsaiData.MousePos.rawValue;
    if correctLag 
        % Correct raw mouse position time for lag 
        Time = bonsaiData.MousePos.rawArduinoTime - lagShift;
    else
        % Uncorrected
        Time = bonsaiData.MousePos.rawArduinoTime; 
    end 
    % Interpolate using either the corrected or uncorrected '2PTime'
    % defined above 
    bonsaiData.MousePos.Value = interp1(Time, rawV, bonsai2PTime, 'linear', NaN);
    bonsaiData.MousePos.ArduinoTime = interp1(Time, Time, bonsai2PTime, 'linear', NaN);
end

%% Quadstate (lag corrected)
if isfield(bonsaiData, 'Quadstate')
    rawV = bonsaiData.Quadstate.rawValue;
    if correctLag 
        % Correct raw mouse position time for lag 
        Time = bonsaiData.Quadstate.rawArduinoTime - lagShift;
    else
        % Uncorrected
        Time = bonsaiData.Quadstate.rawArduinoTime; 
    end 
    bonsaiData.Quadstate.Value = interp1(Time, rawV, bonsai2PTime, 'linear', NaN);
    bonsaiData.Quadstate.ArduinoTime = interp1(Time, Time, bonsai2PTime, 'linear', NaN);
end

%% TrialInfo (align to correctedTwoPTime based on raw2PTime)
if isfield(bonsaiData, 'TrialInfo')
    rawStartTimeAll = bonsaiData.TrialInfo.rawArduinoTime;
    raw2PTime = bonsaiData.resampledTimes.rawMainUnique2PTime;
    for n = 1:length(rawStartTimeAll)
        [~, startInd] = min(abs(raw2PTime - rawStartTimeAll(n)));
        if rawStartTimeAll(n) >= raw2PTime(end)
            break;
        end
        bonsaiData.TrialInfo.originalStartTimeAll(n, 1) = raw2PTime(startInd);
        % Again, using the bonsai2PTime corrected or uncorrected time defined above 
        [~, frameIdx] = min(abs(bonsai2PTime - raw2PTime(startInd)));
        bonsaiData.TrialInfo.alignedStartTimeAll(n, 1) = bonsai2PTime(frameIdx);
        bonsaiData.TrialInfo.alignedStartFrameAll(n, 1) = frameIdx;
    end
end

%% Plot sanity checks
if plotFlag && correctLag
    figure('Name', 'Resampled Signals with Lag Correction');
    plotIdx = 1;

    if isfield(peripheralData, 'Photodiode')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Photodiode.ArduinoTime, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('Photodiode: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end

    if isfield(bonsaiData, 'MousePos')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.MousePos.ArduinoTime, bonsaiData.MousePos.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('MousePos: Raw vs Lag-Corrected'); legend; ylabel('Position'); xlabel('Time');
    end

    if isfield(bonsaiData, 'Quadstate')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.Quadstate.rawArduinoTime, bonsaiData.Quadstate.rawValue, 'k--', 'DisplayName', 'raw Time'); hold on;
        plot(bonsaiData.Quadstate.ArduinoTime, bonsaiData.Quadstate.Value, 'b-', 'DisplayName', 'Corrected/Interpolated Time');
%         plot(bonsaiData.correctedTwoPTime, bonsaiData.Quadstate.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Quadstate: Raw vs Lag-Corrected'); legend; ylabel('Quad'); xlabel('Time');
    end

else
    figure('Name', 'Resampled Signals without Lag Correction');
    plotIdx = 1;
    if isfield(peripheralData, 'Photodiode')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.Photodiode.ArduinoTime, bonsaiData.Photodiode.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Photodiode: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end

    if isfield(bonsaiData, 'MousePos')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.MousePos.ArduinoTime, bonsaiData.MousePos.Value, 'r-', 'DisplayName', 'Interpolated');
        title('MousePos: Raw vs Interpolated (No Lag Correction)'); legend; ylabel('Position'); xlabel('Time');
    end

    if isfield(bonsaiData, 'Quadstate')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.Quadstate.rawArduinoTime, bonsaiData.Quadstate.rawValue, 'k--', 'DisplayName', 'Uncorrected Time'); hold on;
        plot(bonsaiData.c, bonsaiData.Quadstate.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Quadstate: Raw vs Lag-Corrected vs Interpolated'); legend; ylabel('Quad'); xlabel('Time');
    end

end

%% Save
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData")
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
