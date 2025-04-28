function [processedTwoPData] = resampleMergedBonsaiSuite2PDataV1(sessionFileInfo, samplingRate, mainTimeToUse, plotFlag)

% TODO: (1) Rename resampledTimes structure to resampled2PTimes to make this
% more clear. 
%
% This function concatenates and resamples calcium imaging signals and associated
% time fields from merged Bonsai-Suite2P data. It aligns all Suite2P planes and
% time fields to a common uniform time base defined by a selected main time 
% variable (e.g., 'ArduinoTime' or 'TwoPFrameTime'), resampled at a fixed rate.
%
% The output `processedTwoPData` includes aligned time fields, interpolated ROI
% activity traces, and all supporting metadata for analysis.
%
% Inputs:
%   - sessionFileInfo (struct): 
%       Struct containing file paths, animal/session identifiers, and save locations.
%
%   - samplingRate (float, optional): 
%       The target sampling rate in Hz for interpolation (default: 60).
%
%   - mainTimeToUse (char, optional): 
%       The field to use as the primary time base ('ArduinoTime' or 'TwoPFrameTime').
%       All other time fields and calcium traces are aligned to this reference.
%
%   - plotFlag (logical or char, optional): 
%       Whether to display summary plots of time interpolation and calcium traces
%       (default: true).
%
% Outputs:
%   - processedTwoPData (struct):
%       Contains:
%         - F, Fneu, spks           : Interpolated calcium traces
%         - Time fields (ArduinoTime, TwoPFrameTime, etc.)
%         - roiPlaneIdentity        : Per-ROI plane index
%         - iscell, redcell, stat   : ROI classification & metadata
%         - ops, planeName          : Per-plane metadata
%         - resamplingInfo          : Struct with sampling info and interpolation basis
%
% Dependencies:
%   - Requires preprocessed merged Suite2P + Bonsai data for the selected stimulus
%   - Should be called after `mergeBonsaiSuite2pFiles` and `processPeripheralFiles`
%
% Example usage:
%   processedTwoPData = resampleMergedBonsaiSuite2PData(sessionFileInfo, 60, 'ArduinoTime', true);
%
% Sonali and Aman - March 2025


if nargin < 2, samplingRate = 60; end
if nargin < 3, mainTimeToUse = 'ArduinoTime'; end
if nargin < 4, plotFlag = true; end

for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp('VRCorr',sessionFileInfo.stimFiles(iStim).name);
end

iStim = find(bonsaiData.isVRstim==1);

if exist(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file')
    load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData')
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
else
    error('Merged Bonsai-Suite2P data not found for the specified session.');
end

% Save output path
stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_processedBonsai2PData_' sessionFileInfo.stimFiles(iStim).name '.mat'];
sessionFileInfo.stimFiles(iStim).processedBonsai2PSuite2pData = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);

%% Create unified time base from the selected main time
mainTime = vertcat(twoPData.(mainTimeToUse));
[uniqueTime, uniqueIdx] = unique(mainTime);
TwoPTime = interp1(uniqueTime, uniqueTime, uniqueTime(1):1/samplingRate:uniqueTime(end), 'linear')';

%% Resample time-related fields
timeFields = {'TwoPFrameTime', 'ArduinoTime'}; % Excluding RenderFrameCount and LastSyncPulseTime
unalignedTimeFields = {'BonsaiTime', 'RenderFrameCount', 'LastSyncPulseTime'}; % Bonsai clock and Arduino clocks wont work
fprintf('Processing TwoP Frame Times')
% Arduino Times 
for thisField = 1:numel(timeFields)
    fieldName = timeFields{thisField};

    if strcmp(fieldName, mainTimeToUse)
        % Skip interpolation for the main reference time; Interpolated
        % above. 
        processedTwoPData.(fieldName) = TwoPTime;
        continue
    end
    concatenated = vertcat(twoPData.(fieldName));
    concatenated = concatenated(uniqueIdx); % Use the mainTimeToUse index since these data were saved concurrently. 
    % Interpolate the other time-stamps with itself to match the size of
    % the twoPTime. @Aman 
    processedTwoPData.(fieldName) = interp1(concatenated, concatenated, TwoPTime(1):1/samplingRate:TwoPTime(end), 'linear')';
end

% Simply resample other unaligned time fields to 60Hz to have the same
% length? Should this match the lenght of the twoPTime? Or have this
% interpolated to 60Hz with its on start and end times. @Aman 
for thisOtherField = 1:numel(unalignedTimeFields)
    otherFieldName = unalignedTimeFields{thisOtherField};
    otherConcatenated = vertcat(twoPData.(otherFieldName));
%     [~,otherUniqueIdx] = unique(otherConcatenated);
    otherConcatenated = otherConcatenated(uniqueIdx); 
    processedTwoPData.(otherFieldName) = interp1(uniqueTime, otherConcatenated, uniqueTime(1):1/samplingRate:uniqueTime(end), 'linear')';
end

%% Resample Suite2p traces
fprintf('Processing Suite2P Data')
roiFields = {'F', 'Fneu', 'spks'};
interpMethods = {'linear', 'linear', 'nearest'};
for thisField = 1:numel(roiFields)
    processedTwoPData.(roiFields{thisField}) = [];
end
processedTwoPData.roiPlaneIdentity = [];
processedTwoPData.iscell = [];
processedTwoPData.redcell = [];
processedTwoPData.stat = {}; 
processedTwoPData.ops = {};
processedTwoPData.planeName = {};

for p = 1:numel(twoPData)
    % Using the mainTimeToUse as the time vector for the signal
    originalTime = double(twoPData(p).(mainTimeToUse));

    for thisField = 1:numel(roiFields)
        field = roiFields{thisField};
        signal = double(twoPData(p).(field));
        interpolated = interp1(originalTime, signal', TwoPTime(1):1/samplingRate:TwoPTime(end), interpMethods{thisField})';
        processedTwoPData.(field) = [processedTwoPData.(field); interpolated];
    end

    processedTwoPData.iscell = [processedTwoPData.iscell; twoPData(p).iscell];
    processedTwoPData.redcell = [processedTwoPData.redcell; twoPData(p).redcell];
    processedTwoPData.roiPlaneIdentity = [processedTwoPData.roiPlaneIdentity; repmat(p-1, size(twoPData(p).F,1), 1)];
    processedTwoPData.stat = [processedTwoPData.stat, twoPData(p).stat];
    processedTwoPData.ops{end+1} = twoPData(p).ops;
    processedTwoPData.planeName{end+1} = twoPData(p).planeName;
end

%% Save resampling info
processedTwoPData.resampled2PTimes = struct(...
    'samplingRate', samplingRate, ...
    'mainTimeToUse', mainTimeToUse, ...
    'rawMain2PTime', mainTime, ...
    'rawMainUnique2PTime', uniqueTime, ...
    'uniqueIdx', uniqueIdx, ...
    'TwoPTime', TwoPTime);

bonsaiData.resampled2PTimes = processedTwoPData.resampled2PTimes;
peripheralData.resampled2PTimes = processedTwoPData.resampled2PTimes;

%% Optional plotting
if plotFlag
    figure;
    hold on;
    histogram(diff(mainTime), 'BinWidth', 0.001, 'DisplayName', 'Original');
    histogram(diff(uniqueTime), 'BinWidth', 0.001, 'DisplayName', 'Unique');
    histogram(diff(TwoPTime), 'BinWidth', 0.001, 'DisplayName', 'Resampled');
    xlabel('Time Diff (s)'); ylabel('Count');
    legend; title('Time Distribution');
    xlim([0 0.2]);

    
   % Neuron trace comparison
    planeIndex = 1;
    fOrig = double(twoPData(planeIndex).F);
    originalTime = double(twoPData(planeIndex).(mainTimeToUse));
    nROIs = size(fOrig, 1);
    roiIndices = randperm(nROIs, 5);
    roiMask = processedTwoPData.roiPlaneIdentity == (planeIndex - 1);
    fResampled = processedTwoPData.F(roiMask, :);
    
    figure('Name', 'Calcium Trace Comparison');
    for idx = 1:numel(roiIndices)
        roi = roiIndices(idx);
        subplot(numel(roiIndices), 1, idx);
        hold on;
        plot(originalTime, fOrig(roi, :), 'k','LineWidth', 1.2, 'DisplayName', 'Original');
        plot(processedTwoPData.(mainTimeToUse), fResampled(roi, :), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
        title(sprintf('ROI %d (Plane %d)', roi, planeIndex));
        ylabel('F');
        if idx == 1
            legend();
        end
    end
    xlabel('Time (s)');
    sgtitle('Neuron Trace Resampling: Original vs Interpolated');
end


fprintf('Saving processed data files...');
tic; 
save(sessionFileInfo.stimFiles(iStim).processedBonsai2PSuite2pData, 'processedTwoPData', '-v7.3');
fprintf('Saved processedTwoPData\n');
save(sessionFileInfo.stimFiles(iStim).BonsaiData, "bonsaiData");
fprintf('Saved bonsaiData\n');
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData");
fprintf('Saved peripheralData\n');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
fprintf('Saved sessionFileInfo\n');
elapsedTime = toc; 
fprintf('All files saved in %.2f seconds.\n', elapsedTime);
end
