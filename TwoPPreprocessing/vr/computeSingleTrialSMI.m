function sessionMatrix = computeSingleTrialSMI(sessionMatrix, fveThresh)
    if nargin < 2, fveThresh = 0.05; end % Default back to the paper's 5%

    offsetBins = 80;   
    excludeStart = 15; 
    excludeEnd = 15;
    
    for s = 1:length(sessionMatrix)
        if ~isfield(sessionMatrix(s).ConditionData, 'Baseline'), continue; end
        
        data = sessionMatrix(s).ConditionData.Baseline.LapActivity;
        [numROIs, numLaps, numBins] = size(data);
        validRange = (excludeStart + 1) : (numBins - excludeEnd);
        
        % Identify Reliable ROIs using FVE
        fveMask = sessionMatrix(s).FVE >= fveThresh;
        validIdx = find(fveMask);
        
        %  Find "Ground Truth" peak from Mean ODD
        meanOdd = squeeze(mean(data(:, 1:2:end, :), 2, 'omitnan'));
        
        evenIndices = 2:2:numLaps;
        numEven = length(evenIndices);
        
        % Matrix: [numValidROIs x numEven]
        allSingleTrialSMI = NaN(length(validIdx), numEven);
        
        for i = 1:length(validIdx)
            thsROI = validIdx(i);
            
            % Determine preferred bin from Odd average
            [peakVal, relIdx] = max(meanOdd(thsROI, validRange));
            if isnan(peakVal), continue; end
            
            prefBin = relIdx + validRange(1) - 1;
            
            % Identify partner bins
            potentialPartners = prefBin + [-offsetBins, offsetBins];
            matchingBins = potentialPartners(potentialPartners >= 1 & potentialPartners <= numBins);
            if isempty(matchingBins), continue; end
            
            % Calculate SMI for EACH EVEN TRIAL
            for t = 1:numEven
                currentLapIdx = evenIndices(t);
                
                Rp = max(0, data(thsROI, currentLapIdx, prefBin));
                Rn = max(0, mean(data(thsROI, currentLapIdx, matchingBins), 'omitnan'));
                
                if (Rp + Rn) > 1e-9
                    allSingleTrialSMI(i, t) = (Rp - Rn) / (Rp + Rn);
                end
            end
        end
        
        % Store results
        sessionMatrix(s).ConditionData.Baseline.SingleTrialSMI = allSingleTrialSMI;
        
        validValues = allSingleTrialSMI(~isnan(allSingleTrialSMI));
        fprintf('Session %d (FVE > %.2f): Median Single-Trial SMI: %.3f (n=%d trials)\n', ...
            s, fveThresh, median(validValues), length(validValues));
    end
end