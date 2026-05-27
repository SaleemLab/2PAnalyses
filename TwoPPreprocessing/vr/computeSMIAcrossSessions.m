function sessionMatrix = computeSMIAcrossSessions(sessionMatrix)
    % TODO: Includes a position range to detect peaks and then looks for the
    % -+80cm partner. 
    % landmark centres are 40, 80, 120 and 160
    
    partnerLandmarkPosition = 80;   
    landmarkCentres = [40, 80, 120, 160];
    peakWindow = 12;
    % exclude 30% of each lap 
    excludeStart = 30; 
    excludeEnd = 30;
    
    for s = 1:length(sessionMatrix)
        if ~isfield(sessionMatrix(s), 'ConditionData') || ...
           ~isfield(sessionMatrix(s).ConditionData, 'Baseline') || ...
           isempty(sessionMatrix(s).ConditionData.Baseline)
            continue; 
        end
        
        data = sessionMatrix(s).ConditionData.Baseline.LapActivity;
        [numROIs, ~, numBins] = size(data);
        validLapRange = (excludeStart + 1) : (numBins - excludeEnd);
        
        meanOdd = squeeze(mean(data(:, 1:2:end, :), 2, 'omitnan'));
        meanEven = squeeze(mean(data(:, 2:2:end, :), 2, 'omitnan'));
        
        smiValues = NaN(numROIs, 1);
        rpValues  = NaN(numROIs, 1);
        rnValues  = NaN(numROIs, 1);
        peakBins  = NaN(numROIs, 1);
        
        for thsROI = 1:numROIs
            landmarkBins = [];
            for iLand = 1:length(landmarkCentres)
                lowBound = max(1, landmarkCentres(iLand) - peakWindow);
                highBound = min(numBins, landmarkCentres(iLand) + peakWindow);
                landmarkBins = [landmarkBins, lowBound:highBound];
            end
            
            % Intersect landmark windows with the central 30% exclusion window
            validPeakBins = intersect(unique(landmarkBins), validLapRange);
            if isempty(validPeakBins)
                continue;
            end
            
            % Determine peak position from Odd
            [~, maxRelIdx] = max(meanOdd(thsROI, validPeakBins));
            prefBin = validPeakBins(maxRelIdx);
            peakBins(thsROI) = prefBin;
            
            % Measure in even; Response preferred (Rp)
            Rp = max(0, meanEven(thsROI, prefBin));
            rpValues(thsROI) = Rp;
            
            % Split point rule: Midpoint of active track (100cm) dictates direction
            if prefBin <= 100
                targetPartner = prefBin + partnerLandmarkPosition; % Look forward (+80cm)
            else
                targetPartner = prefBin - partnerLandmarkPosition; % Look backward (-80cm)
            end
            
            % Ensure the target partner falls cleanly within the total track boundaries
            if targetPartner >= 1 && targetPartner <= numBins
                matchingBin = targetPartner;
                Rn = max(0, meanEven(thsROI, matchingBin));
                rnValues(thsROI) = Rn;
            else
                continue; 
            end
            
            % Compute SMI
            if (Rp + Rn) > 1e-9
                smiValues(thsROI) = (Rp - Rn) / (Rp + Rn);
            else
                smiValues(thsROI) = NaN; % Assign NaN if SMI cannot be computed
            end
        end
        
        %
        sessionMatrix(s).ConditionData.Baseline.SMI = smiValues;
        
        % 
        sessionMatrix(s).ConditionData.Baseline.SMI_Metrics.SMI = smiValues;
        sessionMatrix(s).ConditionData.Baseline.SMI_Metrics.Rp = rpValues;
        sessionMatrix(s).ConditionData.Baseline.SMI_Metrics.Rn = rnValues;
        sessionMatrix(s).ConditionData.Baseline.SMI_Metrics.PeakBinPosition = peakBins;
        
        numKept = sum(~isnan(smiValues));
        fprintf('Session %d: %d/%d ROIs evaluated successfully.\n', ...
            s, numKept, numROIs);
    end
end