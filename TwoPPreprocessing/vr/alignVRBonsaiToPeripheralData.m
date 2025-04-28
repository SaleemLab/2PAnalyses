function [bonsaiData, sessionFileInfo] = alignVRBonsaiToPeripheralData(sessionFileInfo,VRStimName,plotFlag)
%
% Align Bonsai (ArduinoTime) to peripheal data using the the lag (currently only for
% method 1) 
%
% Inputs:
%   - sessionFileInfo (struct): contains session file paths and metadata.
%   - plotFlag (logical, optional): plot sanity checks of resampled signals (default: true).
%
% Outputs: % 
%   - bonsaiData (struct):
%       - .MousePos / .Quadstate / .TrialInfo:
%           - .correctedArduinoTime: lag-corrected raw time vector.
%       - .LagInfo
%           - .lagShift: Lag shift used in seconds. 
%
% Notes:
%   - Only the arduino Times are aligned using: correctedTime = rawTime - (xcorrBestLag / samplingRate)
%
% Usage:
%   [bonsaiData] = alignVRBonsaiToPeripheralData(sessionFileInfo);
%
% Aman and Sonali - March 2025

if nargin < 3, plotFlag = true; end

for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
end

iStim = find(bonsaiData.isVRstim == 1);

if exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') 
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
else
    error('Missing Bonsai or peripheral data file');
end

% Lag shift for Bonsai correction
lagShift = bonsaiData.LagInfo.xcorrBestLag * (1 / bonsaiData.LagInfo.samplingRate);
bonsaiData.LagInfo.lagShift = lagShift;

%% Mouse Position
if isfield(bonsaiData, 'MousePos')
    bonsaiData.MousePos.rawCorrectedArduinoTime = bonsaiData.MousePos.rawArduinoTime - lagShift; %%% substitute 'rawC' with 'c'
end

%% Quadstate
if isfield(bonsaiData, 'Quadstate')
    bonsaiData.Quadstate.rawCorrectedArduinoTime =  bonsaiData.Quadstate.rawArduinoTime - lagShift;%%% substitute 'rawC' with 'c'
end

%% TrialInfo 
if isfield(bonsaiData, 'TrialInfo')
    bonsaiData.TrialInfo.rawCorrectedArduinoTime = bonsaiData.TrialInfo.rawArduinoTime - lagShift;  %%% substitute 'rawC' with 'c'
end

%% Plot sanity checks
if plotFlag
    figure('Name', 'Lag Correction Sanity Plots');
    if isfield(bonsaiData, 'MousePos')
        subplot(3,1,1); 
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.MousePos.rawCorrectedArduinoTime, bonsaiData.MousePos.rawValue, 'r-', 'DisplayName', ' Corrected');
        title('MousePos: Raw vs Lag-Corrected'); legend; ylabel('Position'); xlabel('Time (s)');
    end

    if isfield(bonsaiData, 'Quadstate')
        subplot(3,1,2);
        plot(bonsaiData.Quadstate.rawArduinoTime, bonsaiData.Quadstate.rawValue, 'k--', 'DisplayName', 'raw Time'); hold on;
        plot(bonsaiData.Quadstate.rawCorrectedArduinoTime, bonsaiData.Quadstate.rawValue, 'r-', 'DisplayName', 'Corrected');
        title('Quadstate: Raw vs Lag-Corrected'); legend; ylabel('Quad'); xlabel('Time (s)');
    end

    if isfield(bonsaiData, 'TrialInfo')
        subplot(3,1,3);
        scatter(bonsaiData.TrialInfo.rawArduinoTime, repmat(6, size(bonsaiData.TrialInfo.rawArduinoTime)),'k', 'DisplayName', 'raw Time'); hold on;
        scatter(bonsaiData.TrialInfo.rawCorrectedArduinoTime, repmat(6, size(bonsaiData.TrialInfo.rawCorrectedArduinoTime)),'red', 'filled', 'DisplayName', 'Corrected'); hold on;
        title('TrialStart: Raw vs Lag-Corrected'); legend; ylabel('Lap Start'); xlabel('Time (s)');
    end

end

%% Save
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end
