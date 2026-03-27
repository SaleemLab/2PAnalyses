%%This version excludes bad frames 
function [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVisualStim_BonsaiPeripheralSuite2P(sessionFileInfo, samplingRate, mainTimeToUse, StimName, plotFlag, trimNaNs, overwrite)
%   Aligns and interpolates all VR experiment signals (Suite2P, Bonsai, peripheral)
%   to a unified 2P timebase using the Arduino clock.
%
%   - Overwrites processedTwoPData file as a FLATTENED file (-struct).
%   - Appends/Updates bonsaiData and peripheralData as nested structs.
%   - Excludes bad frames from ops. 
%   - Pupil tracking included (03/26)

%% Define default parameters
if nargin < 2, samplingRate = 60; end
if nargin < 3, mainTimeToUse = 'TwoPFrameTime'; end 
if nargin < 5, plotFlag = true; end
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
if isfield(sessionFileInfo.stimFiles(thisStim), 'BonsaiData')
    fileBonsai = sessionFileInfo.stimFiles(thisStim).BonsaiData;
else 
    fileBonsai = [];
end 
filePeripheral = sessionFileInfo.stimFiles(thisStim).processedPeripheralData;

% Define Output File Path (processedTwoPData)
stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_processed2PData_' sessionFileInfo.stimFiles(thisStim).name '.mat'];
outputSavePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
sessionFileInfo.stimFiles(thisStim).processedMergedBonsaiSuite2pData = outputSavePath;

%% Overwrite Check (Updated for Flattened Loading)
if exist(outputSavePath, 'file') && ~overwrite
    disp(['Processed 2P Data for ' StimName ' already exists. Loading flattened file...']);
    processedTwoPData = load(outputSavePath);
    if exist(fileBonsai, 'file'), load(fileBonsai, 'bonsaiData'); else, bonsaiData = []; end
    if exist(filePeripheral, 'file'), load(filePeripheral, 'peripheralData'); else, peripheralData = []; end
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
processedTwoPData.badFrames = {}; % Initialized as cell array
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
    
    % Interpolate badFrames mask (0/1)
    rawMask = double(twoPData(thisPlaneIdx).badframes);
    interpMask = interp1(rawArduinoPlaneTime, rawMask, sampleTimes, 'linear', 0);
    
    % Force to logical row vector to prevent memory errors later
    processedTwoPData.badFrames{thisPlaneIdx} = logical(interpMask(:)' > 0)';

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
if isfield(peripheralData, 'Photodiode') && ~isempty(peripheralData.Photodiode)
    rawTime = peripheralData.Photodiode.rawArduinoTime;
    rawValue = peripheralData.Photodiode.rawValue;
    peripheralData.Photodiode.Value = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Photodiode.sampleTimes = sampleTimes';
end

%% Interpolate: Peripheral - Quadstate 
disp('Processing Peripheral Data: Quad')
if isfield(peripheralData, 'Quadstate') && ~isempty(peripheralData.Quadstate)
    rawValue = peripheralData.Quadstate.rawValue;
    rawTime = peripheralData.Quadstate.rawArduinoTime;
    [uniqueRawTime, idx] = unique(rawTime);
    uniqueRawValue = double(rawValue(idx));
    
    Value = interp1(uniqueRawTime, uniqueRawValue, sampleTimes, 'previous', NaN)'; 
    peripheralData.Quadstate.Value = Value;
    peripheralData.Quadstate.sampleTimes = sampleTimes';
end

%%  Interpolate: Peripheral - Pupil (with clock alignment)  @ Written based on SGS alignEyeData() 02/26
disp('Processing Peripheral Data: Pupil')
if isfield(peripheralData, 'Pupil')
    
    % Sync pupil arduino clock to the 2P arduino clock
    uSyncEye = unique(peripheralData.Pupil.raw.LastSyncPulseTime);
    
    % We'll use the first plane's sync pulses as the reference; will all be
    % the same
    uSyncTwoP = unique(twoPData(1).LastSyncPulseTime);
    
    % Align the pupil timestamps to the 2P timestamps
    newPupilArduinoTime = align2PSyncPulses(uSyncEye, uSyncTwoP, peripheralData.Pupil.raw.ArduinoTime);
    
    % Interpolate Pupil features to sample times 
    % Note: We interpolate from the previously calculated .int fields  
    
    pupilFields = {'CentroidX', 'CentroidY', 'Area', 'MajorAxisLength', 'MinorAxisLength'};
    
    for thisfld = 1:numel(pupilFields)
        fName = pupilFields{thisfld};
        rawValue = peripheralData.Pupil.int.(fName);
        
        % values; keep consistent
        peripheralData.Pupil.Value.(fName) = interp1(newPupilArduinoTime, rawValue, sampleTimes, generalInterpMethod, 'extrap')'; 
    end
    
    peripheralData.Pupil.sampleTimes = sampleTimes';

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
        
        % Trim badFrames cell array
        for p = 1:numel(processedTwoPData.badFrames)
            processedTwoPData.badFrames{p} = processedTwoPData.badFrames{p}(combinedValidMask);
        end

        if isfield(peripheralData, 'Wheel')
            peripheralData.Wheel.Value = peripheralData.Wheel.Value(combinedValidMask);
            peripheralData.Wheel.sampleTimes = peripheralData.Wheel.sampleTimes(combinedValidMask);
        end
        if isfield(peripheralData, 'Photodiode') && ~isempty(peripheralData.Photodiode)
            peripheralData.Photodiode.Value = peripheralData.Photodiode.Value(combinedValidMask);
            peripheralData.Photodiode.sampleTimes = peripheralData.Photodiode.sampleTimes(combinedValidMask);
        end
        if isfield(peripheralData, 'Quadstate') && ~isempty(peripheralData.Photodiode)
            peripheralData.Quadstate.Value = peripheralData.Quadstate.Value(combinedValidMask);
            peripheralData.Quadstate.sampleTimes = peripheralData.Quadstate.sampleTimes(combinedValidMask);
        end

        % New
        if isfield(peripheralData, 'Pupil')
            pFields = fieldnames(peripheralData.Pupil.Value);
            for f = 1:numel(pFields)
                peripheralData.Pupil.Value.(pFields{f}) = peripheralData.Pupil.Value.(pFields{f})(combinedValidMask);
            end
            peripheralData.Pupil.sampleTimes = peripheralData.Pupil.sampleTimes(combinedValidMask);
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
        sampleTimes = newSampleTimes; % Update global sampleTimes for plotting
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
    
    % Bad Frame Masking Effect Example Plot
    figure('Name', ['Example: Bad Frame Masking Effect: ' StimName]);
    p_idx = 1; 
    roi_to_plot = find(processedTwoPData.roiPlaneIdentity == (p_idx-1), 1);
    if ~isempty(roi_to_plot)
        subplot(2,1,1);
        originalInterpolated = processedTwoPData.F(roi_to_plot, :);
        tempNanTrace = originalInterpolated;
        tempNanTrace(processedTwoPData.badFrames{p_idx}) = NaN;
        plot(sampleTimes, originalInterpolated, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'F (Interpolated)');
        hold on;
        plot(sampleTimes, tempNanTrace, 'r', 'LineWidth', 1.5, 'DisplayName', 'F (If NaNs applied)');
        title(sprintf('ROI %d Plane %d: Bad Frame Removal visualization', roi_to_plot, p_idx)); ylabel('F'); legend;
        subplot(2,1,2);
        plot(sampleTimes, processedTwoPData.badFrames{p_idx}, 'k');
        title('Bad Frame Mask (1 = Bad)'); ylabel('Status'); xlabel('Time (s)'); ylim([-0.1 1.1]);
    end

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
    if isfield(peripheralData, 'Photodiode') && ~isempty(peripheralData.Photodiode)
        nexttile;
        plot(peripheralData.Photodiode.sampleTimes, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Resampled');
        title('Photodiode'); legend;
    end
    if isfield(peripheralData, 'Wheel')
        nexttile;
        plot(peripheralData.Wheel.sampleTimes, peripheralData.Wheel.Value, 'r-', 'DisplayName', 'Resampled');
        title('Wheel'); legend;
    end

    %% pupil plot 
    if isfield(peripheralData, 'Pupil') && isfield(peripheralData, 'Wheel')
        figure('Name', 'Alignment Check: Pupil vs Wheel');
        tickToCmConversion = 3.1415 * 20 / 1024;  % Wheel radius 20 cm, 1024 ticks per revolution
        displacement = [0; diff(peripheralData.Wheel.Value * tickToCmConversion)];

        displacement(displacement < -100) = 0;  % Negative large jumps
        displacement(displacement > 100) = 0;   % Positive large jumps


        wheelSpeed = displacement ./ [0; diff(peripheralData.Wheel.sampleTimes)]; 
        smoothWin = 15; 
        
        pClean = movmean(peripheralData.Pupil.Value.Area, smoothWin, 'omitnan');
        wClean = movmean(wheelSpeed, smoothWin, 'omitnan');


        zP = (pClean - mean(pClean, 'omitnan')) / std(pClean, 'omitnan');
        zW = (wClean - mean(wClean, 'omitnan')) / std(wClean, 'omitnan');
        hold on;

        % Plot Wheel (Black) and Pupil (Red)
        plot(zW, 'k', 'LineWidth', 1.5, 'DisplayName', 'Running Speed');
        plot(zP, 'r', 'LineWidth', 2, 'DisplayName', 'Pupil Area');

        totalSamples = length(zP);
        fiveMinSamples = 60 * 5 * 60; % (60Hz * 60sec * 5min) = 18000

        startSample = 2000; 

        if totalSamples > (startSample + fiveMinSamples)
        
            endSample = startSample + fiveMinSamples;
        else
   
            endSample = totalSamples;
           
        end

        % 3. Apply to Plot
        legend('Location', 'northwest');
        xlabel('Samples');
        ylabel('Activity (Z-Score)');
        title('Pupil vs Running Speed');
        xlim([startSample, endSample]);
        ylim([-2 5]);
        
    end 
end

%% Save Section (Flattened processedTwoPData only)
disp('Saving processed data files...');
save(outputSavePath, '-struct', 'processedTwoPData', '-v7.3');
if exist(fileBonsai, 'file'), save(fileBonsai, 'bonsaiData', '-append'); end
if exist(filePeripheral, 'file'), save(filePeripheral, 'peripheralData', '-append'); end
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Done.');
end

% function [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVisualStim_BonsaiPeripheralSuite2P(sessionFileInfo, samplingRate, mainTimeToUse, StimName, plotFlag, trimNaNs, overwrite)
% %   Aligns and interpolates all VR experiment signals (Suite2P, Bonsai, peripheral)
% %   to a unified 2P timebase using the Arduino clock.
% %
% %   - Overwrites processedTwoPData file as a FLATTENED file (-struct).
% %   - Appends/Updates bonsaiData and peripheralData as nested structs.
% 
% %% Define default parameters
% if nargin < 2, samplingRate = 60; end
% if nargin < 3, mainTimeToUse = 'TwoPFrameTime'; end 
% if nargin < 5, plotFlag = false; end
% if nargin < 6, trimNaNs = true; end
% if nargin < 7, overwrite = true; end 
% 
% %% File Paths
% bonsaiData = []; 
% isStim = strcmp(StimName, {sessionFileInfo.stimFiles.name});
% thisStim = find(isStim, 1);
% if isempty(thisStim)
%     error('Stimulus %s not found in sessionFileInfo', StimName);
% end
% 
% % Define Input File Paths
% fileMerged = sessionFileInfo.stimFiles(thisStim).mergedBonsai2PSuite2pData;
% fileBonsai = sessionFileInfo.stimFiles(thisStim).BonsaiData;
% filePeripheral = sessionFileInfo.stimFiles(thisStim).processedPeripheralData;
% 
% % Define Output File Path (processedTwoPData)
% stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_processed2PData_' sessionFileInfo.stimFiles(thisStim).name '.mat'];
% outputSavePath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% sessionFileInfo.stimFiles(thisStim).processedMergedBonsaiSuite2pData = outputSavePath;
% 
% %% Overwrite Check (Updated for Flattened Loading)
% if exist(outputSavePath, 'file') && ~overwrite
%     disp(['Processed 2P Data for ' StimName ' already exists. Loading flattened file...']);
% 
%     % Reconstruct the struct from individual variables in the flat file
%     processedTwoPData = load(outputSavePath);
% 
%     if exist(fileBonsai, 'file')
%         load(fileBonsai, 'bonsaiData');
%     else
%         bonsaiData = [];
%     end
% 
%     if exist(filePeripheral, 'file')
%         load(filePeripheral, 'peripheralData');
%     else
%         peripheralData = [];
%     end
% 
%     disp('Data loaded from file. Skipping processing.');
%     return; 
% end
% 
% %% Load raw data streams 
% disp(['Processing ' StimName '...']);
% if ~exist(fileMerged, 'file'), error('Merged data file not found: %s', fileMerged); end
% if ~exist(filePeripheral, 'file'), error('Peripheral data file not found: %s', filePeripheral); end
% 
% load(fileMerged, 'twoPData');
% load(filePeripheral, 'peripheralData');
% if ~exist(fileBonsai, 'file')
%     warning('Bonsai data file not found, continuing without it: %s', fileBonsai);
% else
%     load(fileBonsai, 'bonsaiData');
% end
% 
% %% Create unified time base 
% raw2PTimes = vertcat(twoPData.(mainTimeToUse));
% [unique2PTimes, mainTimeUniqueIdx] = unique(raw2PTimes);
% resample2PTimes = unique2PTimes(1):1/samplingRate:unique2PTimes(end);
% sampleTimes = resample2PTimes;
% generalInterpMethod     = 'linear';
% trialInfoInterpMethod   = 'nearest'; 
% 
% %% Interpolate: Two-photon time vectors
% timeFields = {'TwoPFrameTime', 'ArduinoTime'}; 
% disp('Processing TwoP Frame Times')
% for thisField = 1:numel(timeFields)
%     fieldName = timeFields{thisField};
%     concatenatedTimeVec = vertcat(twoPData.(fieldName));
%     concatenatedTimeVec = concatenatedTimeVec(mainTimeUniqueIdx);
%     processedTwoPData.(fieldName) = interp1(concatenatedTimeVec, concatenatedTimeVec, sampleTimes, generalInterpMethod)';
% end
% processedTwoPData.resample2PTimeUsed = mainTimeToUse; 
% 
% %% Interpolate: Suite2p data
% disp('Processing Suite2P Data')
% roiFields = {'F', 'Fneu', 'spks'};
% for thisField = 1:numel(roiFields), processedTwoPData.(roiFields{thisField}) = []; end
% processedTwoPData.roiPlaneIdentity = [];
% processedTwoPData.iscell = [];
% processedTwoPData.redcell = [];
% processedTwoPData.stat = {};
% processedTwoPData.ops = {};
% processedTwoPData.planeName = {};
% 
% for thisPlaneIdx = 1:numel(twoPData)
%     rawArduinoPlaneTime = double(twoPData(thisPlaneIdx).(mainTimeToUse)); 
%     for thisField = 1:numel(roiFields)
%         fieldName = roiFields{thisField};
%         signal = double(twoPData(thisPlaneIdx).(fieldName));
%         interpolated = interp1(rawArduinoPlaneTime, signal', sampleTimes, generalInterpMethod)';
%         processedTwoPData.(fieldName) = [processedTwoPData.(fieldName); interpolated];
%     end
%     processedTwoPData.iscell = [processedTwoPData.iscell; twoPData(thisPlaneIdx).iscell];
%     processedTwoPData.redcell = [processedTwoPData.redcell; twoPData(thisPlaneIdx).redcell];
%     processedTwoPData.roiPlaneIdentity = [processedTwoPData.roiPlaneIdentity; repmat(thisPlaneIdx-1, size(twoPData(thisPlaneIdx).F,1), 1)];
%     processedTwoPData.stat = [processedTwoPData.stat, twoPData(thisPlaneIdx).stat];
%     processedTwoPData.ops{end+1} = twoPData(thisPlaneIdx).ops;
%     processedTwoPData.planeName{end+1} = twoPData(thisPlaneIdx).planeName;
% end
% 
% %% Interpolate: Peripheral - Wheel
% disp('Processing Peripheral Data: Wheel')
% if isfield(peripheralData, 'Wheel')
%     rawTime     = peripheralData.Wheel.rawArduinoTime;
%     rawValue    = peripheralData.Wheel.rawValue;
%     peripheralData.Wheel.Value      = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
%     peripheralData.Wheel.sampleTimes = sampleTimes';
% end
% 
% %% GrayScreen Check: Remove fields if necessary 
% if contains(StimName, 'GrayScreen', 'IgnoreCase',true)
%     if exist('peripheralData', 'var')
%         peripheralfieldsToRemove = {'Quadstate', 'Photodiode'};
%         existingFieldsPeripheral = intersect(fieldnames(peripheralData), peripheralfieldsToRemove);
%         if ~isempty(existingFieldsPeripheral)
%             peripheralData = rmfield(peripheralData, existingFieldsPeripheral);
%             fprintf('GrayScreen detected. In peripheralData, removed fields: %s\n', strjoin(existingFieldsPeripheral, ', '));
%         end
%     end
% 
%     if ~isempty(bonsaiData)
%         bonsaiFieldsToRemove = {'MousePos'};
%         existingFieldsBonsai = intersect(fieldnames(bonsaiData), bonsaiFieldsToRemove);
%         if ~isempty(existingFieldsBonsai)
%             bonsaiData = rmfield(bonsaiData, existingFieldsBonsai);
%             fprintf('GrayScreen detected. In bonsaiData, removed fields: %s\n', strjoin(existingFieldsBonsai, ', '));
%         end
%     end 
% end
% 
% %% Interpolate: Peripheral - Photodiode
% disp('Processing Peripheral Data: PD')
% if isfield(peripheralData, 'Photodiode')
%     rawTime = peripheralData.Photodiode.rawArduinoTime;
%     rawValue = peripheralData.Photodiode.rawValue;
%     peripheralData.Photodiode.Value = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
%     peripheralData.Photodiode.sampleTimes = sampleTimes';
% end
% 
% %% Interpolate: Peripheral - Quadstate 
% disp('Processing Peripheral Data: Quad')
% if isfield(peripheralData, 'Quadstate')
%     rawValue = peripheralData.Quadstate.rawValue;
%     rawTime = peripheralData.Quadstate.rawArduinoTime;
%     [uniqueRawTime, idx] = unique(rawTime);
%     uniqueRawValue = double(rawValue(idx));
% 
%     Value = interp1(uniqueRawTime, uniqueRawValue, sampleTimes, 'previous', NaN)'; 
%     peripheralData.Quadstate.Value = Value;
%     peripheralData.Quadstate.sampleTimes = sampleTimes';
% end
% 
% %% Interpolate: Bonsai - TrialInfo
% disp('Processing Bonsai Data: StimOnset Info')
% fieldsToProcess = { 'ArduinoTimeRaw', 'ArduinoTime' };
% for thisField = 1:size(fieldsToProcess, 1)
%     rawField = fieldsToProcess{thisField, 1};
%     targetField = fieldsToProcess{thisField, 2};
%     if isfield(bonsaiData, rawField)
%         rawTime = bonsaiData.(rawField);
%         processedValues = interp1(sampleTimes, sampleTimes, rawTime, trialInfoInterpMethod); 
%         bonsaiData.(targetField) = processedValues;
%         fprintf('Processed %s into %s\n', rawField, targetField);
%     end
% end
% 
% %% Trim NaN padding
% if trimNaNs
%     disp('Trimming NaN padding...');
%     combinedValidMask = true(size(sampleTimes));
% 
%     combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.F), 1);
%     combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.Fneu), 1);
%     combinedValidMask = combinedValidMask & all(isfinite(processedTwoPData.spks), 1);
% 
%     if ~any(combinedValidMask)
%         warning('Trimming NaNs would remove all data. Skipping trim.');
%     else
%         processedTwoPData.trimmedNaNPadding = true;
%         newSampleTimes = sampleTimes(combinedValidMask);
% 
%         processedTwoPData.F = processedTwoPData.F(:, combinedValidMask);
%         processedTwoPData.Fneu = processedTwoPData.Fneu(:, combinedValidMask);
%         processedTwoPData.spks = processedTwoPData.spks(:, combinedValidMask);
%         processedTwoPData.TwoPFrameTime = processedTwoPData.TwoPFrameTime(combinedValidMask);
%         processedTwoPData.ArduinoTime = processedTwoPData.ArduinoTime(combinedValidMask);
% 
%         if isfield(peripheralData, 'Wheel')
%             peripheralData.Wheel.Value = peripheralData.Wheel.Value(combinedValidMask);
%             peripheralData.Wheel.sampleTimes = peripheralData.Wheel.sampleTimes(combinedValidMask);
%         end
%         if isfield(peripheralData, 'Photodiode')
%             peripheralData.Photodiode.Value = peripheralData.Photodiode.Value(combinedValidMask);
%             peripheralData.Photodiode.sampleTimes = peripheralData.Photodiode.sampleTimes(combinedValidMask);
%         end
%         if isfield(peripheralData, 'Quadstate')
%             peripheralData.Quadstate.Value = peripheralData.Quadstate.Value(combinedValidMask);
%             peripheralData.Quadstate.sampleTimes = peripheralData.Quadstate.sampleTimes(combinedValidMask);
%         end
% 
%         % Filter TrialInfo data
%         if isfield(bonsaiData, 'ArduinoTime') && ~isempty(bonsaiData.ArduinoTime)
%             referenceTime = bonsaiData.ArduinoTime;
%             minTime = newSampleTimes(1);
%             maxTime = newSampleTimes(end);
%             keepTrials = (referenceTime >= minTime) & (referenceTime <= maxTime);
%             nTotalTrials = numel(referenceTime);
% 
%             infoFields = fieldnames(bonsaiData);
%             for f = 1:numel(infoFields)
%                 if numel(bonsaiData.(infoFields{f})) == nTotalTrials
%                     bonsaiData.(infoFields{f}) = bonsaiData.(infoFields{f})(keepTrials);
%                 end
%             end
%         end
%     end
% end
% 
% %% Sanity check plots
% if plotFlag
%     % 2P Times
%     figure('Name', ['TwoP Arduino Frame Times: ' StimName]);
%     hold on;
%     histogram(diff(raw2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Original');
%     histogram(diff(unique2PTimes), 'BinWidth', 0.001, 'DisplayName', 'Unique');
%     histogram(diff(processedTwoPData.(mainTimeToUse)), 'BinWidth', 0.001, 'DisplayName', 'Resampled');
%     xlabel('Time Diff (s)'); ylabel('Count'); legend; title('2P Frame Time Distribution'); xlim([0 0.2]);
% 
%     % Trace Comparison
%     planeIndex = 1;
%     roiMask = processedTwoPData.roiPlaneIdentity == (planeIndex - 1);
%     if any(roiMask)
%         fResampled = processedTwoPData.F(roiMask, :);
%         fOrig = double(twoPData(planeIndex).F);
%         originalTime = double(twoPData(planeIndex).(mainTimeToUse));
%         if size(fResampled,1) < 5, roiIndices = 1:size(fResampled,1); else, roiIndices = randperm(size(fResampled,1), 5); end
%         figure('Name', 'Calcium Trace Comparison');
%         for idx = 1:numel(roiIndices)
%             roi = roiIndices(idx);
%             subplot(numel(roiIndices), 1, idx); hold on;
%             plot(originalTime, fOrig(roi, :), 'k','LineWidth', 1.2, 'DisplayName', 'Original');
%             plot(processedTwoPData.(mainTimeToUse), fResampled(roi, :), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
%             ylabel('F'); if idx==1, legend; end
%         end
%         sgtitle('Neuron Trace Resampling');
%     end
% 
%     % Peripheral Plots
%     figure('Name', ['Peripheral Signals: ' StimName]); t = tiledlayout('flow');
%     if isfield(peripheralData, 'Photodiode')
%         nexttile;
%         plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
%         plot(peripheralData.Photodiode.sampleTimes, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Resampled');
%         title('Photodiode'); legend;
%     end
%     if isfield(peripheralData, 'Wheel')
%         nexttile;
%         plot(peripheralData.Wheel.rawArduinoTime, peripheralData.Wheel.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
%         plot(peripheralData.Wheel.sampleTimes, peripheralData.Wheel.Value, 'r-', 'DisplayName', 'Resampled');
%         title('Wheel'); legend;
%     end
%     if isfield(peripheralData, 'Quadstate')
%         nexttile;
%         plot(peripheralData.Quadstate.rawArduinoTime, peripheralData.Quadstate.rawValue, 'k', 'DisplayName', 'Raw'); hold on;
%         plot(peripheralData.Quadstate.sampleTimes, peripheralData.Quadstate.Value, 'r-', 'DisplayName', 'Resampled');
%         title('Quadstate'); legend;
%     end
% end
% 
% %% Save Section (Flattened processedTwoPData only)
% disp('Saving processed data files...');
% 
% % Save processedTwoPData: Unpack struct into individual variables
% save(outputSavePath, '-struct', 'processedTwoPData', '-v7.3');
% disp(['Saved Flattened processedTwoPData to: ' outputSavePath]);
% 
% % bonsaiData: Standard append as a struct
% if exist(fileBonsai, 'file')
%     save(fileBonsai, 'bonsaiData', '-append');
%     disp('Updated bonsaiData (Appended Struct)');
% end
% 
% % peripheralData: Standard append as a struct
% if exist(filePeripheral, 'file')
%     save(filePeripheral, 'peripheralData', '-append');
%     disp('Updated peripheralData (Appended Struct)');
% end
% 
% % Update session info
% save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
% disp('Done.');
% end