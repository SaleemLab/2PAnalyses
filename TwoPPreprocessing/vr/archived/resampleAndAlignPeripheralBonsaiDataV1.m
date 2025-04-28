function [peripheralData, bonsaiData] = resampleAndAlignPeripheralBonsaiData(sessionFileInfo, correctLag, plotFlag)
% PD on and offs missing; possibly not important but save anyway.
% Interpolate all other 'time' vectors saved for each Bonsai and
% Peripheral field - currently only focused on Arduino Times. 
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
% Outputs:
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

if nargin < 2, correctLag = true; end
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

TwoPTime = peripheralData.resampledTimes.TwoPTime;
samplingRate = peripheralData.resampledTimes.samplingRate;
peripheralData.TwoPTime = TwoPTime;

% Lag shift for Bonsai correction
lagShift = 0;
if correctLag && isfield(bonsaiData, 'LagInfo') && isfield(bonsaiData.LagInfo, 'xcorrBestLag')
    lagShift = bonsaiData.LagInfo.xcorrBestLag * (1 / samplingRate);
end

bonsaiData.TwoPTime = TwoPTime;
bonsaiData.correctedTwoPTime = TwoPTime - lagShift;

%% Wheel (no lag correction)
if isfield(peripheralData, 'Wheel')
    rawT = peripheralData.Wheel.rawArduinoTime;
    rawV = peripheralData.Wheel.rawValue;
    peripheralData.Wheel.Value = interp1(rawT, rawV, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
    peripheralData.Wheel.ArduinoTime = interp1(rawT, rawT, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
end

%% Photodiode (no lag correction)
if isfield(peripheralData, 'Photodiode')
    rawT = peripheralData.Photodiode.rawArduinoTime;
    rawV = peripheralData.Photodiode.rawValue;
    peripheralData.Photodiode.Value = interp1(rawT, rawV, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
    peripheralData.Photodiode.ArduinoTime = interp1(rawT, rawT, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear', NaN);
end

%% Mouse Position (lag corrected)
if isfield(bonsaiData, 'MousePos')
    rawT = bonsaiData.MousePos.rawArduinoTime;
    rawV = bonsaiData.MousePos.rawValue;
    correctedT = rawT - lagShift;
    bonsaiData.MousePos.correctedArduinoTime = correctedT;
    bonsaiData.MousePos.Value = interp1(correctedT, rawV, bonsaiData.correctedTwoPTime, 'linear', NaN);
    bonsaiData.MousePos.ArduinoTime = interp1(correctedT, correctedT, bonsaiData.correctedTwoPTime, 'linear', NaN);
end

%% Quadstate (lag corrected)
if isfield(bonsaiData, 'Quadstate')
    rawT = bonsaiData.Quadstate.rawArduinoTime;
    rawV = bonsaiData.Quadstate.rawValue;
    correctedT = rawT - lagShift;
    bonsaiData.Quadstate.uncorrectedArduinoTime = rawT;
    bonsaiData.Quadstate.correctedArduinoTime = correctedT;
    bonsaiData.Quadstate.Value = interp1(correctedT, rawV, bonsaiData.correctedTwoPTime, 'linear', NaN);
    bonsaiData.Quadstate.ArduinoTime = interp1(correctedT, correctedT, bonsaiData.correctedTwoPTime, 'linear', NaN);
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
        [~, frameIdx] = min(abs(bonsaiData.correctedTwoPTime - raw2PTime(startInd)));
        bonsaiData.TrialInfo.alignedStartTimeAll(n, 1) = bonsaiData.correctedTwoPTime(frameIdx);
        bonsaiData.TrialInfo.alignedStartFrameAll(n, 1) = frameIdx;
    end
end

%% Plot sanity checks
if plotFlag
    figure('Name', 'Resampled Signals with/without Lag Correction');
    plotIdx = 1;

    if isfield(peripheralData, 'Photodiode')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(TwoPTime, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Photodiode: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end

    if isfield(bonsaiData, 'MousePos')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.correctedTwoPTime, bonsaiData.MousePos.Value, 'r-', 'DisplayName', 'Corrected/Interpolated');
        title('MousePos: Raw vs Lag-Corrected'); legend; ylabel('Position'); xlabel('Time');
    end

    if isfield(bonsaiData, 'Quadstate')
        subplot(4,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.Quadstate.uncorrectedArduinoTime, bonsaiData.Quadstate.rawValue, 'k--', 'DisplayName', 'Uncorrected Time'); hold on;
        plot(bonsaiData.Quadstate.correctedArduinoTime, bonsaiData.Quadstate.rawValue, 'b-', 'DisplayName', 'Corrected Time');
        plot(bonsaiData.correctedTwoPTime, bonsaiData.Quadstate.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Quadstate: Raw vs Lag-Corrected vs Interpolated'); legend; ylabel('Quad'); xlabel('Time');
    end
end

%% Save
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData")
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
