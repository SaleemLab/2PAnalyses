function [alignedResult] = getSpeedStratifiedTuning(sessionFileInfo, VRStimName, targetNeuron)
    % Load response using your shared logic
    stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}), 1);
    if isempty(stimIdx), error('Stimulus not found'); end
    
    load(sessionFileInfo.stimFiles(stimIdx).Response, 'lapRunningSpeed', 'lapPositionActivity', 'lapPositionRunningSpeed');
    

    % Access your specific activity matrix [Neurons x Laps x Position]
    % We take [targetNeuron, all laps, all positions]
    ROISignal = squeeze(lapPositionActivity.dFFNeuropilCorrected(targetNeuron, :, :)); % [57 x 200]
    
    % Access the speed matrix (assuming it's binned the same way)
    % If not, we calculate mean speed per lap from response.lapRunningSpeed
    lapSpeeds = cellfun(@mean, lapRunningSpeed); % [57 x 1]
    
    % Define Speed Thresholds (Paper: 2, 10, 30)
    lowIdx  = lapSpeeds > 2  & lapSpeeds < 10;
    medIdx  = lapSpeeds >= 10 & lapSpeeds < 20;
    highIdx = lapSpeeds >= 20;
    
    % Calculate Mean Tuning for each speed bin
    tuningCurve(1, :) = nanmean(ROISignal(lowIdx, :), 1);
    tuningCurve(2, :) = nanmean(ROISignal(medIdx, :), 1);
    tuningCurve(3, :) = nanmean(ROISignal(highIdx, :), 1);
    
    % --- Alignment Logic ---
    % 1. Find the peak position using Medium speed
    [~, peakBin] = max(tuningCurve(2, :));
    
    % 2. Shift so peak is at the center (Bin 100)
    numBins = 200;
    shiftAmount = 100 - peakBin;
    alignedResult.data = circshift(tuningCurve, shiftAmount, 2);
    
    % 3. Set up X-axis (assuming 1cm bins for 200cm corridor)
    alignedResult.xAxis = (1:numBins) - 100; 
    alignedResult.speeds = [mean(lapSpeeds(lowIdx)), mean(lapSpeeds(medIdx)), mean(lapSpeeds(highIdx))];
end