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
%   Modified Nov: Creates and saves sessionROIData 
%
% Aman and Sonali - April 2025

%% Define default paramters and load appropriate data files
if nargin < 2, samplingRate = 60; end
if nargin < 3, mainTimeToUse = 'TwoPFrameTime'; end 
if nargin < 5, plotFlag = false; end
if nargin < 6, trimNaNs = true; end 

% Pick out VRCorr
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName,sessionFileInfo.stimFiles(iStim).name);
end
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

%% Create unified time base
raw2PTimes = vertcat(twoPData.(mainTimeToUse));
[unique2PTimes, mainTimeUniqueIdx] = unique(raw2PTimes);
sampleTimes = unique2PTimes(1):1/samplingRate:unique2PTimes(end);
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
for thisField = 1:numel(roiFields)
    processedTwoPData.(roiFields{thisField}) = [];
end
processedTwoPData.badFrames = {}; 
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

    % Force to logical to prevent memory overflow (logical is 8x smaller than double)
    processedTwoPData.badFrames{thisPlaneIdx} = logical(interpMask > 0);

    processedTwoPData.iscell = [processedTwoPData.iscell; twoPData(thisPlaneIdx).iscell];
    processedTwoPData.redcell = [processedTwoPData.redcell; twoPData(thisPlaneIdx).redcell];
    processedTwoPData.roiPlaneIdentity = [processedTwoPData.roiPlaneIdentity; repmat(thisPlaneIdx-1, size(twoPData(thisPlaneIdx).F,1), 1)];
    processedTwoPData.stat = [processedTwoPData.stat, twoPData(thisPlaneIdx).stat];
    processedTwoPData.ops{end+1} = twoPData(thisPlaneIdx).ops;
    processedTwoPData.planeName{end+1} = twoPData(thisPlaneIdx).planeName;
end

%% Interpolate: Peripheral and Bonsai Data
disp('Processing Peripheral Data: Wheel')
if isfield(peripheralData, 'Wheel')
    rawTime     = peripheralData.Wheel.rawArduinoTime;
    rawValue    = peripheralData.Wheel.rawValue;
    peripheralData.Wheel.Value      = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Wheel.sampleTimes = sampleTimes';
end

disp('Processing Peripheral Data: PD')
if isfield(peripheralData, 'Photodiode')
    rawTime = peripheralData.Photodiode.rawArduinoTime;
    rawValue = peripheralData.Photodiode.rawValue;
    peripheralData.Photodiode.Value = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Photodiode.sampleTimes = sampleTimes';
end

disp('Processing Bonsai Data: Quad')
if isfield(peripheralData, 'Quadstate')
    rawValue = peripheralData.Quadstate.rawValue;
    lagCorrT = peripheralData.Quadstate.rawCorrectedArduinoTime;
    peripheralData.Quadstate.Value = interp1(lagCorrT, rawValue, sampleTimes, generalInterpMethod, NaN)';
    peripheralData.Quadstate.sampleTimes = sampleTimes';
end


%%  Interpolate: Peripheral - Pupil (with clock alignment)  @ Written based on SGS alignEyeData() 02/26
disp('Processing Peripheral Data: Pupil')
if isfield(peripheralData, 'Pupil')
    % Sync pupil arduino clock to the 2P arduino clock
    uSyncEye = unique(peripheralData.Pupil.raw.LastSyncPulseTime);
    % We'll use the first plane's sync pulses as the reference; will all be
    % the same if multi-plane recording..wait no it wont @SONALI TODO 
    uSyncTwoP = unique(twoPData(1).LastSyncPulseTime);
    % Align the pupil timestamps to the 2P timeline
    % This handles any offset or drift between the two Arduinos 
    newPupilArduinoTime = align2PSyncPulses(uSyncEye, uSyncTwoP, peripheralData.Pupil.raw.ArduinoTime);
    % Interpolate Pupil features to sample times 
    % Note: We interpolate from the previously calculated .int fields 
    pupilFields = {'CentroidX', 'CentroidY', 'Area', 'MajorAxisLength', 'MinorAxisLength'};
    for thisfld = 1:numel(pupilFields)
        fName = pupilFields{thisfld};
        rawValue = peripheralData.Pupil.int.(fName);
        % values; keep consistent
        % what does extrap do to nans? are there even nans in pupil area?
        % @Sonali TODO
        peripheralData.Pupil.Value.(fName) = interp1(newPupilArduinoTime, rawValue, sampleTimes, generalInterpMethod, 'extrap')';
    end

    peripheralData.Pupil.sampleTimes = sampleTimes';

end

disp('Processing Bonsai Data: Mouse Position')
if isfield(bonsaiData, 'MousePos')
    rawValue = bonsaiData.MousePos.rawValue;
    rawTime = bonsaiData.MousePos.rawCorrectedArduinoTime;
    correctedStartTimeAll = bonsaiData.TrialInfo.rawCorrectedArduinoTime;
    newValue = interp1(rawTime, rawValue, sampleTimes, generalInterpMethod, NaN)';
    newTime = sampleTimes';
    if sum(diff(rawTime)>2) > 3 
        tidx = find([diff(rawTime) ;0]>2); 
        if sum(rawValue(tidx) < 0)>0 
             for i = 1:length(correctedStartTimeAll)
                 candidateFrame = find(rawTime(tidx) < correctedStartTimeAll(i));
                 if isempty(candidateFrame)
                     candidateFrame = 1;
                 end
                 [~,idx]= min(abs(correctedStartTimeAll(i)-rawTime(tidx(candidateFrame))));
                 [~,idx1]= min(abs(newTime-rawTime(tidx(idx))));
                 [~,idx2]= min(abs(newTime-correctedStartTimeAll(i)));
                 newValue(idx1:idx2) = newValue(idx1);
             end
        end
    end
    bonsaiData.MousePos.Value = newValue;
    bonsaiData.MousePos.sampleTimes = newTime;
end

disp('Processing Bonsai Data: Trial Info')
if isfield(bonsaiData, 'TrialInfo')
    correctedStartTimeAll = bonsaiData.TrialInfo.rawCorrectedArduinoTime;
    uncorrectedStartTimeAll = bonsaiData.TrialInfo.rawArduinoTime;
    bonsaiData.TrialInfo.StartTimeAll = interp1(sampleTimes, sampleTimes, correctedStartTimeAll, trialInfoInterpMethod);
    bonsaiData.TrialInfo.uncorrectedStartTimeAll = interp1(sampleTimes, sampleTimes, uncorrectedStartTimeAll, trialInfoInterpMethod);
end

%% Trim NaN padding if requested 
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
        processedTwoPData.trimmedMetaData = ['Trimming ' num2str(numel(sampleTimes)) ' samples down to ' num2str(sum(combinedValidMask)) ' samples.'];
        newSampleTimes = sampleTimes(combinedValidMask);
        
        % suite2p 
        processedTwoPData.F = processedTwoPData.F(:, combinedValidMask);
        processedTwoPData.Fneu = processedTwoPData.Fneu(:, combinedValidMask);
        processedTwoPData.spks = processedTwoPData.spks(:, combinedValidMask);
        % trim and used later down the pipeline 
        for p = 1:numel(processedTwoPData.badFrames)
            processedTwoPData.badFrames{p} = processedTwoPData.badFrames{p}(combinedValidMask);
        end
        

        processedTwoPData.TwoPFrameTime = processedTwoPData.TwoPFrameTime(combinedValidMask);
        processedTwoPData.ArduinoTime = processedTwoPData.ArduinoTime(combinedValidMask);
        sampleTimes = newSampleTimes;
        
        % New pupile 
        if isfield(peripheralData, 'Pupil')
            pFields = fieldnames(peripheralData.Pupil.Value);
            for f = 1:numel(pFields)
                peripheralData.Pupil.Value.(pFields{f}) = peripheralData.Pupil.Value.(pFields{f})(combinedValidMask);
            end
            peripheralData.Pupil.sampleTimes = peripheralData.Pupil.sampleTimes(combinedValidMask);
        end
        % wheel & times 
        if isfield(peripheralData, 'Wheel')
            peripheralData.Wheel.Value = peripheralData.Wheel.Value(combinedValidMask); 
            peripheralData.Wheel.sampleTimes = peripheralData.Wheel.sampleTimes(combinedValidMask); 
        end
        % pd & times 
        if isfield(peripheralData, 'Photodiode')
            peripheralData.Photodiode.Value = peripheralData.Photodiode.Value(combinedValidMask);
            peripheralData.Photodiode.sampleTimes = peripheralData.Photodiode.sampleTimes(combinedValidMask);
        end
        % mouse position & times 
        if isfield(bonsaiData, 'MousePos')
            bonsaiData.MousePos.Value = bonsaiData.MousePos.Value(combinedValidMask);
            bonsaiData.MousePos.sampleTimes = bonsaiData.MousePos.sampleTimes(combinedValidMask); 
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

    %  EXAMPLE: Show what NaNs would look like (Temporary visualize ONLY)
    figure('Name', 'Example: Bad Frame Masking Effect');
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
        title(sprintf('ROI %d: Visualization of Bad Frames', roi_to_plot));
        ylabel('Signal'); legend;
        subplot(2,1,2);
        plot(sampleTimes, processedTwoPData.badFrames{p_idx}, 'k');
        title('Bad Frame Mask (1 = Bad)'); ylabel('Status'); xlabel('Time (s)');
        ylim([-0.1 1.1]);
    end

    % Neuron trace comparison
    planeIndex = 1;
    fOrig = double(twoPData(planeIndex).F);
    originalTime = double(twoPData(planeIndex).(mainTimeToUse));
    roiMask = processedTwoPData.roiPlaneIdentity == (planeIndex - 1);
    if any(roiMask)
        fResampled = processedTwoPData.F(roiMask, :);
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
            plot(originalTime, fOrig(roi, :), 'k','LineWidth', 1.2, 'DisplayName', 'Original');
            plot(processedTwoPData.(mainTimeToUse), fResampled(roi, :), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Resampled');
            title(sprintf('ROI %d (Plane %d)', roi, planeIndex));
            ylabel('F'); if idx == 1; legend(); end
        end
        xlabel('Time (s)'); sgtitle('Neuron Trace Resampling: Original vs Interpolated');
    else
        disp('Skipping calcium trace plot: No ROIs found for plane 1.');
    end

    figure('Name', 'Resampled Peripheral Signals');
    plotIdx = 1;
    if isfield(peripheralData, 'Photodiode')
        subplot(3,1,plotIdx); 
        plot(peripheralData.Photodiode.rawArduinoTime, peripheralData.Photodiode.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Photodiode.sampleTimes, peripheralData.Photodiode.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Photodiode: Raw vs Interpolated'); legend; ylabel('PD'); xlabel('Time');
    end
    if isfield(peripheralData, 'Wheel')
        plotIdx = plotIdx + 1;
        subplot(3,1,plotIdx);
        plot(peripheralData.Wheel.rawArduinoTime, peripheralData.Wheel.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(peripheralData.Wheel.sampleTimes, peripheralData.Wheel.Value, 'r-', 'DisplayName', 'Interpolated');
        title('Wheel: Raw vs Interpolated'); legend; ylabel('Value'); xlabel('Time');
    end

    if isfield(peripheralData, 'Pupil') && isfield(peripheralData, 'Wheel')
        plotIdx = plotIdx + 1;
        subplot(3, 1, plotIdx); 
        figure('Name', 'Alignment Check: Pupil vs Running Speed');
        tickToCmConversion = 3.1415 * 20 / 1024;  % Wheel radius 20 cm, 1024 ticks per revolution
        displacement = [0; diff(peripheralData.Wheel.Value * tickToCmConversion)];

        displacement(displacement < -100) = 0;  % Negative large jumps
        displacement(displacement > 100) = 0;   % Positive large jumps
        wheelSpeed = displacement ./ [0; diff(peripheralData.Wheel.sampleTimes)];
        smoothWin = 15;
        %similar to pd
        pClean = movmean(peripheralData.Pupil.Value.Area, smoothWin, 'omitnan');
        wClean = movmean(wheelSpeed, smoothWin, 'omitnan');
        %z-scoring to compare; unsure if this is correct @Sonali TODO
        % zP = (pClean - mean(pClean, 'omitnan')) / std(pClean, 'omitnan');
        % zW = (wClean - mean(wClean, 'omitnan')) / std(wClean, 'omitnan');
        hold on;

        % Plot Wheel (Black) and Pupil (Red)
        yyaxis left 
        plot(wClean, 'k', 'LineWidth', 1.5, 'DisplayName', 'Running Speed');
        yyaxis right 
        plot(pClean, 'r', 'LineWidth', 2, 'DisplayName', 'Pupil Area');
        % some recordings (gray screen/darkness can be shorter than 5mins)
        totalSamples = length(pClean);
        fiveMinSamples = 60 * 5 * 60; % (60Hz * 60sec * 5min) = 18000
        startSample = 2000;
        if totalSamples > (startSample + fiveMinSamples)
            endSample = startSample + fiveMinSamples;
        else
            endSample = totalSamples;
        end
        legend('Location', 'northwest');
        xlabel('Samples');
        ylabel('Signals');
        title('Pupil vs Running Speed');
        xlim([startSample, endSample]);

    end
    figure('Name', 'Resampled Bonsai Signals with Lag Correction');
    plotIdx = 1;
    if isfield(bonsaiData, 'MousePos')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.MousePos.rawArduinoTime, bonsaiData.MousePos.rawValue, 'k.', 'DisplayName', 'Raw'); hold on;
        plot(bonsaiData.MousePos.sampleTimes, bonsaiData.MousePos.Value, 'r-', 'DisplayName', 'Interpolated');
        title('MousePos: Raw vs Lag-Corrected'); legend; ylabel('Position'); xlabel('Time');
    end
    if isfield(bonsaiData, 'Quadstate')
        subplot(2,1,plotIdx); plotIdx = plotIdx + 1;
        plot(bonsaiData.Quadstate.rawArduinoTime, bonsaiData.Quadstate.rawValue, 'k', 'DisplayName', 'raw Time'); hold on;
        plot(bonsaiData.Quadstate.sampleTimes, bonsaiData.Quadstate.Value, 'r-', 'DisplayName', 'Interpolated Time');
        title('Quadstate: Raw vs Lag-Corrected'); legend; ylabel('Quad'); xlabel('Time');
    end
end

%% Save
save(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData,'-struct', 'processedTwoPData', '-v7.3');
save(sessionFileInfo.stimFiles(iStim).BonsaiData, "bonsaiData");
save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData");
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Saved processed data with badFrames mask included.');
end
