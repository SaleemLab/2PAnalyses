function sessionMetrics = getTuningData(sessionTable, varargin)
  
    % getTuningData: Loads VR data, computes MEAN tuning curves 
    % and Modulation Index for each session in the input table.
    % Sonali Dec 2025
    
    p = inputParser;
    addRequired(p, 'sessionTable', @istable); 
    addParameter(p, 'signalToUse', 'dFF', @ischar); 
    addParameter(p, 'applySmoothing', true, @islogical);
    parse(p, sessionTable, varargin{:});
    params = p.Results;
    
    numSessions = height(sessionTable);
    sessionMetrics = struct('MouseID', {}, 'Day', {}, 'Session', {}, 'Type', {}, ...
                         'InjectionSite', {}, 'TargetArea', {}, 'OddMean', {}, 'EvenMean', {}, ...
                         'MeanTuning', {}, 'Modulation', {}, 'NumCells', {}, 'NumLaps', {}, ...
                         ... % Fields from sessionROIData
                         'TypeImaged_ROI', {}, 'TargetArea_ROI', {}, ...
                         'highlyCorrBoutons', {}, ...
                         ... % Fields from lapCorr_Halves
                         'lapCorr_HalvesRho', {}, 'lapCorr_HalvesP', {}, 'lapCorr_HalvesStableThreshold', {}, ...
                         'lapCorr_HalvesStableIdx', {}, 'lapCorr_HalvesFirstHalfLapsCount', {}, 'lapCorr_HalvesSecondHalfLapsCount', {}, ...
                         ... % Fields from lapCorr_OddEven
                         'lapCorr_OddEvenRho', {}, 'lapCorr_OddEvenP', {}, 'lapCorr_OddEvenStableThreshold', {}, ...
                         'lapCorr_OddEvenStableIdx', {}, ...
                         ... % Fields from nullDist_PeakTuningMetric
                         'isSignificantByPeakShuffling', {} ,'realPeakPercentileRank', {}, 'targetPercentile_Peak', {}, ...
                         ... % Fields from nullDist_RangeTuningMetric
                         'isSignificantByRange', {}, 'realRangePercentileRank', {}, 'targetPercentile_Range', {}, ...
                         ... % Fields from tuningCurveVariance
                         'ratioVarToTuningVar', {}, 'ratioVarToTuningRange', {} ...
                         );
    
    % wb = waitbar(0, 'Loading and processing sessions...');
    
    for thisSessions = 1:numSessions
        % waitbar(thisSessions/numSessions, wb, sprintf('Processing session %d/%d...', thisSessions, numSessions));
        
        row = sessionTable(thisSessions,:);
        currentSessionData = struct(); 
        
        currentSessionData.MouseID = row.MouseID{1};
        currentSessionData.Day = row.DayOfExperience;
        currentSessionData.Session = char(row.Session);
        
        % Type identification
        if ismember('TypeImaged', row.Properties.VariableNames)
             currentSessionData.Type = row.TypeImaged{1};
        elseif ismember('ImagedType', row.Properties.VariableNames)
             currentSessionData.Type = row.ImagedType{1};
        else
             currentSessionData.Type = 'Unknown';
        end
        
        if ismember('GCaMPInjectionSite', row.Properties.VariableNames)
            currentSessionData.InjectionSite = row.GCaMPInjectionSite{1};
            currentSessionData.TargetArea = row.TargetArea{1};
        else
             currentSessionData.InjectionSite = '';
             currentSessionData.TargetArea = '';
        end
        
        % Load Activity Data
        [lapActivity, numLaps, sfi, errorMessage] = loadVRLapSignal(currentSessionData.MouseID, currentSessionData.Session, params.signalToUse);
        [sessionROIData] = loadSessionROIData(sfi); 
        
        if isempty(lapActivity)
            warning('Skipping %s Session %s Day %d: %s', currentSessionData.MouseID, currentSessionData.Session, currentSessionData.Day, errorMessage);
            continue;
        end
        
        % Smooth data 
        if params.applySmoothing
            smoothingKernel = gausswin(9); 
            smoothingKernel = smoothingKernel / sum(smoothingKernel);
            lapActivity = smoothActivityTraces(lapActivity, smoothingKernel);
        end
        
        % Tuning curves
        meanOddLaps = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
        meanEvenLaps = squeeze(mean(lapActivity(:, 2:2:end, :), 2, 'omitnan'));
        meanAllLaps = squeeze(mean(lapActivity, 2, 'omitnan'));
        
        % Modulation Index
        tuningCurveMax = max(meanAllLaps, [], 2, 'omitnan'); 
        tuningCurveMin = min(meanAllLaps, [], 2, 'omitnan');
        tuningCurveMean = mean(meanAllLaps, 2, 'omitnan');
        modulationIndex = (tuningCurveMax - tuningCurveMin) ./ tuningCurveMean;
        modulationIndex(tuningCurveMean == 0) = NaN; 
        
        currentSessionData.OddMean = meanOddLaps;
        currentSessionData.EvenMean = meanEvenLaps;
        currentSessionData.MeanTuning = meanAllLaps;
        currentSessionData.Modulation = modulationIndex;
        currentSessionData.NumCells = size(lapActivity, 1);
        currentSessionData.NumLaps = numLaps;
        
        % Safe extraction from sessionROIData
        % Metadata
        currentSessionData.TypeImaged_ROI = safeGet(sessionROIData, 'typeImaged');
        currentSessionData.TargetArea_ROI = safeGet(sessionROIData, 'targetArea');
        
        % Highly Correlated Boutons
        if isfield(sessionROIData, 'highlyCorrBoutons')
            currentSessionData.highlyCorrBoutons = sessionROIData.highlyCorrBoutons.roisToKeep;
        else
            currentSessionData.highlyCorrBoutons = [];
        end
  
        % lapCorr_Halves
        if isfield(sessionROIData, 'lapCorr_Halves')
            currentSessionData.lapCorr_HalvesRho = sessionROIData.lapCorr_Halves.rho;
            currentSessionData.lapCorr_HalvesP = sessionROIData.lapCorr_Halves.p;
            currentSessionData.lapCorr_HalvesStableThreshold = sessionROIData.lapCorr_Halves.stableThreshold;
            currentSessionData.lapCorr_HalvesStableIdx = sessionROIData.lapCorr_Halves.stableIdx;
            currentSessionData.lapCorr_HalvesFirstHalfLapsCount = sessionROIData.lapCorr_Halves.firstHalfLapsCount;
            currentSessionData.lapCorr_HalvesSecondHalfLapsCount = sessionROIData.lapCorr_Halves.secondHalfLapsCount;
        else
            currentSessionData.lapCorr_HalvesRho = []; currentSessionData.lapCorr_HalvesP = [];
            currentSessionData.lapCorr_HalvesStableThreshold = []; currentSessionData.lapCorr_HalvesStableIdx = [];
            currentSessionData.lapCorr_HalvesFirstHalfLapsCount = []; currentSessionData.lapCorr_HalvesSecondHalfLapsCount = [];
        end
        
        % lapCorr_OddEven
        if isfield(sessionROIData, 'lapCorr_OddEven')
            currentSessionData.lapCorr_OddEvenRho = sessionROIData.lapCorr_OddEven.rho;
            currentSessionData.lapCorr_OddEvenP = sessionROIData.lapCorr_OddEven.p;
            currentSessionData.lapCorr_OddEvenStableThreshold = sessionROIData.lapCorr_OddEven.stableThreshold;
            currentSessionData.lapCorr_OddEvenStableIdx = sessionROIData.lapCorr_OddEven.stableIdx;
        else
            currentSessionData.lapCorr_OddEvenRho = []; currentSessionData.lapCorr_OddEvenP = [];
            currentSessionData.lapCorr_OddEvenStableThreshold = []; currentSessionData.lapCorr_OddEvenStableIdx = [];
        end
        
        % nullDist_PeakTuningMetric
        if isfield(sessionROIData, 'nullDist_PeakTuningMetric')
            currentSessionData.isSignificantByPeakShuffling = safeGet(sessionROIData.nullDist_PeakTuningMetric.isSignificantByPeakShuffling, params.signalToUse);
            currentSessionData.realPeakPercentileRank = safeGet(sessionROIData.nullDist_PeakTuningMetric.realPeakPercentileRank, params.signalToUse);
            currentSessionData.targetPercentile_Peak = sessionROIData.nullDist_PeakTuningMetric.targetPercentile;
        else
            currentSessionData.isSignificantByPeakShuffling = []; currentSessionData.realPeakPercentileRank = []; currentSessionData.targetPercentile_Peak = [];
        end
        
        % nullDist_RangeTuningMetric
        if isfield(sessionROIData, 'nullDist_RangeTuningMetric')
            currentSessionData.isSignificantByRange = safeGet(sessionROIData.nullDist_RangeTuningMetric.isSignificantByRange, params.signalToUse);
            currentSessionData.realRangePercentileRank = safeGet(sessionROIData.nullDist_RangeTuningMetric.realRangePercentileRank, params.signalToUse);
            currentSessionData.targetPercentile_Range = sessionROIData.nullDist_RangeTuningMetric.targetPercentile; 
        else
            currentSessionData.isSignificantByRange = []; currentSessionData.realRangePercentileRank = []; currentSessionData.targetPercentile_Range = [];
        end
        
        % tuningCurveVariance
        if isfield(sessionROIData, 'tuningCurveVariance')
            currentSessionData.ratioVarToTuningVar = sessionROIData.tuningCurveVariance.ratioVarToTuningVar;
            currentSessionData.ratioVarToTuningRange = sessionROIData.tuningCurveVariance.ratioVarToTuningRange;
        else
            currentSessionData.ratioVarToTuningVar = []; currentSessionData.ratioVarToTuningRange = [];
        end
        
        sessionMetrics(end+1) = orderfields(currentSessionData, sessionMetrics); 
    end
    
    % close(wb);
    fprintf('Data loaded for %d sessions.\n', length(sessionMetrics));
end

function val = safeGet(s, field)
    if isstruct(s) && isfield(s, field)
        val = s.(field);
    else
        val = [];
    end
end



%% Load sessionROIData (fixed but double check dec 2025 SSS)
function [sessionROIData] = loadSessionROIData(sessionFileInfo)
% FIX: Initialize output to an empty structure. This ensures the output 
% is assigned even if an error occurs, preventing the main function crash.
sessionROIData = struct(); 

try
    sessionROIDataFilePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
    
    if exist(sessionROIDataFilePath, 'file')
        sessionROIData = load(sessionROIDataFilePath);
    else
        warning('loadSessionROIData:MissingFile', ...
            'SessionROIData.mat is missing at path: %s', ...
            sessionROIDataFilePath);
    end
    
catch ME

    sessionROIData = struct(); 
    warning('loadSessionROIData:Error', ...
        'An error occurred while trying to load SessionROIData. Check if "otherSessFilePaths" exists in sessionFileInfo. Error: %s', ...
        ME.message);
end
end


%% local helper function (loadVRLapSignal)
function [lapActivity, numLaps, sfi, errorMessage] = loadVRLapSignal(subjectID, sessionID, signalName)
    lapActivity = [];
    numLaps = 0;
    errorMessage = '';
    sfi = struct(); %sfi just in case of a crash before assignment
    
    try
        % Find the path to the session's info file
        sessionInfoPath = findSessionFileInfoFilePath(subjectID, sessionID);
        if ~isfile(sessionInfoPath)
            errorMessage = 'InfoMissing';
            return;
        end
        
        % Load the session file info
        D = load(sessionInfoPath, 'sessionFileInfo');
        sfi = D.sessionFileInfo; % 
        
        % Find the correct VR stimulus file
        stimFileNames = string({sfi.stimFiles.name});
        vrFileIndex = find(contains(stimFileNames, "BaselineCorridor") & contains(stimFileNames, "CombinedRuns"), 1);
        
        if isempty(vrFileIndex) % Fallback: Find last run if no 'CombinedRuns' Need to find a better fix here 
            vrFileIndex = find(contains(stimFileNames, "BaselineCorridor") & ~contains(stimFileNames, "CombinedRuns"), 1, 'first'); 
        end
        
        if isempty(vrFileIndex)
            errorMessage = 'NoVRCorr';
            return; 
        end
        
        % Get the response file path
        responseFilePath = sfi.stimFiles(vrFileIndex).Response;
        if isempty(responseFilePath) || ~isfile(responseFilePath)
            errorMessage = 'RespMissing'; 
            return; 
        end
        
        % Load the response data
        R = load(responseFilePath);
        
        % Check if the signal exists in the lap data structure
        if ~isfield(R, 'lapPositionActivity') || ~isfield(R.lapPositionActivity, signalName)
            errorMessage = 'SignalMissing'; 
            return;
        end
        
        % Extract the signal data
        lapActivity = R.lapPositionActivity.(signalName);
        numLaps = size(lapActivity, 2);
        
        % Ensure there is more than one lap to be useful
        if numLaps < 2
            lapActivity = [];
            errorMessage = '<2Laps';
            return; 
        end
        
    catch ME
        % Added ME to the catch block to debug if needed
        errorMessage = 'LoadError';
        warning('loadVRLapSignal:LoadError', 'Error: %s', ME.message);
    end
end


%% Smooth activity traces
function smoothedData = smoothActivityTraces(activityData, filterKernel)
    smoothedData = activityData;
    
    % Loop over cells/channels (1st dimension)
    for cellIdx = 1:size(activityData, 1)
        % Loop over laps (2nd dimension)
        for lapIdx = 1:size(activityData, 2)
            
            % Get the trace for this specific cell and lap
            trace = squeeze(activityData(cellIdx, lapIdx, :));
            
            % Find NaN values to restore them after filtering
            nanMask = isnan(trace);
            if all(nanMask)
                continue; % Skip if the entire trace is NaN
            end
            
  
            trace(nanMask) = 0;
            

            filteredTrace = filtfilt(filterKernel, 1, trace);
            
            % Restore NaN values
            filteredTrace(nanMask) = NaN;
            
            % Assign the smoothed trace back to the output array
            smoothedData(cellIdx, lapIdx, :) = filteredTrace;
        end
    end
end
