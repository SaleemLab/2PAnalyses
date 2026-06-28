function updatedMetrics = appendFilteredROIs(sessionMetrics, varargin)
% APPENDFILTEREDROIS Takes the output structure from getTuningDataByCondition,
% appends a new field 'FilteredROIs', and computes global/session ROI filter metrics.
%
% Example usage:
%   RSPData = getTuningDataByCondition(ExpRSPSessions);
%   RSPData = appendFilteredROIs(RSPData, 'RhoHalvesThreshold', 0.6, 'FilterEdgeSMI', true);

    p = inputParser;
    addRequired(p, 'sessionMetrics', @isstruct);
    addParameter(p, 'UseHalves', false, @islogical)
    addParameter(p, 'RhoHalvesThreshold', 0.8, @isnumeric); 
    addParameter(p, 'UseOddEven', false, @islogical)
    addParameter(p, 'RhoOddEvenThreshold', 0.8, @isnumeric);
    addParameter(p, 'UseExpVar_SigNullDist', true, @islogical);
    addParameter(p, 'ExpVarSigThreshold', 0.01, @isnumeric);
    addParameter(p, 'UseExpVar', false, @islogical);
    addParameter(p, 'cvExpvarThreshold', 0.1, @isnumeric); 
    addParameter(p, 'FilterDuplicateBoutons', true, @islogical);
    addParameter(p, 'FilterSomasByRF', false, @islogical);
    addParameter(p, 'RestrictSMIToLandmarkBoundaries', false, @islogical);
    addParameter(p, 'LandmarkBoundaries', 15, @isnumeric); 
    
    % filtering of cells with dominant peaks at track ends
    addParameter(p, 'FilterEdgeSMI', false, @islogical); 
    
    parse(p, sessionMetrics, varargin{:});
    params = p.Results;
    
    updatedMetrics = sessionMetrics;
    numSessions = length(sessionMetrics);
    
    % Define landmark bins to include based on the input boundary 
    landmarkCentres = [40, 80, 120, 160];
    landmarkBinsToCheck = [];
    for iLand = 1:length(landmarkCentres)
        thislandPos = landmarkCentres(iLand);
        landmarkBins = thislandPos-params.LandmarkBoundaries:thislandPos+params.LandmarkBoundaries;
        landmarkBinsToCheck = [landmarkBinsToCheck, landmarkBins];
    end
    TotalROIs = 0;
    KeptROIs = 0;
    FilteredOutROIs = 0;
    
    for thisSess = 1:numSessions
        thisSession = sessionMetrics(thisSess);
        
        if ~isfield(thisSession, 'ConditionData') || isempty(thisSession.ConditionData)
            warning('Session %d (Mouse %s) has no ConditionData. Skipping.', thisSess, thisSession.MouseID);
            updatedMetrics(thisSess).FilteredROIs = [];
            continue; 
        end
        
        condNames = fieldnames(thisSession.ConditionData);
        sampleData = thisSession.ConditionData.(condNames{1}).LapActivity;
        numTotalROIs = size(sampleData, 1); 
        
        roisToKeepIdx = 1:numTotalROIs; 
        
        % filter 1: cross-validated explained variance (cvExpVar) 
        if params.UseExpVar_SigNullDist
            if isfield(thisSession, 'cvExpVar') && ~isempty(thisSession.cvExpVar) 
                realVarIdx = find(thisSession.cvExpVar.pValues <= params.ExpVarSigThreshold); 
                roisToKeepIdx = intersect(roisToKeepIdx, realVarIdx);
            end
        end 
        
        if params.UseExpVar
            if isfield(thisSession, 'cvExpVar') && ~isempty(thisSession.cvExpVar) 
                realVarIdx = find( thisSession.cvExpVar.medianExpVar > params.cvExpvarThreshold);
                roisToKeepIdx = intersect(roisToKeepIdx, realVarIdx);
            end
        end
        
        % filter 2: Halves Rho 
        if params.UseHalves
            if isfield(thisSession, 'Rho_Halves') && ~isempty(thisSession.Rho_Halves)
                stableIdx = find(thisSession.Rho_Halves >= params.RhoHalvesThreshold);
                roisToKeepIdx = intersect(roisToKeepIdx, stableIdx);
            else
                warning('Mouse %s Session %s is missing ''Rho_Halves'' data.', thisSession.MouseID, thisSession.Session);
            end
        end 
        
        % Option Filter 2: Odd Even Rho 
        if params.UseOddEven
            if isfield(thisSession, 'Rho_OddEven') && ~isempty(thisSession.Rho_OddEven)
                stableIdx = find(thisSession.Rho_OddEven >= params.RhoOddEvenThreshold);
                roisToKeepIdx = intersect(roisToKeepIdx, stableIdx);
            else
                warning('Mouse %s Session %s is missing ''Rho_OddEven'' data.', thisSession.MouseID, thisSession.Session);
            end
        end 
        
        % Filter 3: Bouton-specific unique index filtering
        if params.FilterDuplicateBoutons && strcmpi(thisSession.TypeImaged, 'Boutons')
            if isfield(thisSession, 'uniqueBoutonIdx') && ~isempty(thisSession.uniqueBoutonIdx)
                roisToKeepIdx = intersect(roisToKeepIdx, thisSession.uniqueBoutonIdx);
            end
        end
        
        % Filter 4: Receptive Field Location Filter (-70 to -40 Azimuth)
        if params.FilterSomasByRF && strcmpi(thisSession.TypeImaged, 'Somas')
            if isfield(thisSession, 'sparseNoiseMatrix') && ~isempty(thisSession.sparseNoiseMatrix)
                roisWithinWindow = [];
                for neuronIdx = 1:numTotalROIs
                    rfRaw = flipud(thisSession.sparseNoiseMatrix.initMap{neuronIdx}(:, :, end, 4));
                    rfSmoothed = imgaussfilt(rfRaw, 1);
                    azimuthCoords = linspace(-70, 20, size(rfRaw, 2));
                    [~, maxIdx] = max(rfSmoothed(:));
                    [~, maxCol] = ind2sub(size(rfSmoothed), maxIdx);
                    peakAzimuth = azimuthCoords(maxCol);
                    if peakAzimuth >= -70 && peakAzimuth <= -40
                        roisWithinWindow(end+1) = neuronIdx;
                    end
                end
                roisToKeepIdx = intersect(roisToKeepIdx, roisWithinWindow);
            end
        end
        
        % Filter 5: Spatial Modulation Boundary Restrictions
        if params.RestrictSMIToLandmarkBoundaries
            if isfield(thisSession.SMI, 'PeakBinPosition') && ~isempty(thisSession.SMI.PeakBinPosition)
                ROIsWithinLandmarkWindows = find(ismember(thisSession.SMI.PeakBinPosition, landmarkBinsToCheck));
                roisToKeepIdx = intersect(roisToKeepIdx, ROIsWithinLandmarkWindows);
            end 
        end 
        
        % NEW FILTER: Exclude ROIs whose global maximum peak lies in the first/last 30cm
        if params.FilterEdgeSMI
            % Check inside the parsed SMI field/struct for the saved exclusion flags
            if isfield(thisSession, 'SMI') && isfield(thisSession.SMI, 'ExcludeEdgePeakCells') && ~isempty(thisSession.SMI.ExcludeEdgePeakCells)
                % Find indexes of cells that are NOT flagged for edge exclusions
                cleanSMIIdx = find(~thisSession.SMI.ExcludeEdgePeakCells);
                roisToKeepIdx = intersect(roisToKeepIdx, cleanSMIIdx);
            else
                warning('Mouse %s Session %s is missing ''SMI.ExcludeEdgePeakCells'' data.', thisSession.MouseID, thisSession.Session);
            end
        end
        
        updatedMetrics(thisSess).FilteredROIs = roisToKeepIdx;
        
        numKept = length(roisToKeepIdx);
        numFilteredOut = numTotalROIs - numKept;
        sessFilterPct = (numFilteredOut / numTotalROIs) * 100;
        
        TotalROIs = TotalROIs + numTotalROIs;
        KeptROIs = KeptROIs + numKept;
        FilteredOutROIs = FilteredOutROIs + numFilteredOut;
        
        fprintf('Mouse %s Session %s: Kept %d/%d ROIs (%d excluded, %.1f%% filtered out).\n', ...
            thisSession.MouseID, thisSession.Session, numKept, numTotalROIs, numFilteredOut, sessFilterPct);
    end
    
    fprintf('               Filter summary             \n');
    if TotalROIs > 0
        FilterPct = (FilteredOutROIs / TotalROIs) * 100;
        fprintf('Total Raw ROIs Processed: %d\n', TotalROIs);
        fprintf('Total ROIs Retained:      %d\n', KeptROIs);
        fprintf('Total ROIs Filtered Out:  %d\n', FilteredOutROIs);
        fprintf('Overall Exclusion Rate:   %.2f%%\n', FilterPct);
    else
        fprintf('No valid sessions were evaluated for counts.\n');
    end
end