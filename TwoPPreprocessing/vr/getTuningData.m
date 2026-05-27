function sessionMetrics = getTuningData(sessionTable, varargin)
    % getTuningData: Loads VR data, computes MEAN tuning curves 
    % and Modulation Index.
    
    p = inputParser;
    addRequired(p, 'sessionTable', @istable); 
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar); 
    addParameter(p, 'applySmoothing', false, @islogical);
    parse(p, sessionTable, varargin{:});
    params = p.Results;
    
    numSessions = height(sessionTable);
    
    % Define field names once to ensure consistency
    fields = {'MouseID', 'Day', 'Session', 'Type', 'InjectionSite', 'TargetArea', ...
              'OddMean', 'EvenMean', 'MeanTuning', 'Modulation', 'NumCells', 'NumLaps', ...
              'TypeImaged_ROI', 'TargetArea_ROI', 'highlyCorrBoutons', ...
              'lapCorr_HalvesRho', 'lapCorr_HalvesP', 'lapCorr_HalvesStableThreshold', ...
              'lapCorr_HalvesStableIdx', 'lapCorr_HalvesFirstHalfLapsCount', 'lapCorr_HalvesSecondHalfLapsCount', ...
              'lapCorr_OddEvenRho', 'lapCorr_OddEvenP', 'lapCorr_OddEvenStableThreshold', ...
              'lapCorr_OddEvenStableIdx', 'isSignificantByPeakShuffling', 'realPeakPercentileRank', ...
              'targetPercentile_Peak', 'isSignificantByRange', 'realRangePercentileRank', ...
              'targetPercentile_Range', 'ratioVarToTuningVar', 'ratioVarToTuningRange'};

    % Initialize master struct array
    sessionMetrics = cell2struct(cell(length(fields), numSessions), fields, 1);
    
    for thisSessions = 1:numSessions
        row = sessionTable(thisSessions,:);
        
        % Create a fresh struct for this session based on the master fields
        currentSessionData = struct();
        for f = 1:length(fields), currentSessionData.(fields{f}) = []; end
        
        % Basic Info
        currentSessionData.MouseID = char(row.MouseID);
        currentSessionData.Day = row.DayOfExperience;
        currentSessionData.Session = char(row.Session);
        
        fprintf('Processing Mouse %s | Day: %d | Session: %s \n', ...
            currentSessionData.MouseID, currentSessionData.Day, currentSessionData.Session)
        
        % Mapping Table Metadata
        vars = row.Properties.VariableNames;
        if ismember('TypeImaged', vars), currentSessionData.Type = row.TypeImaged{1};
        elseif ismember('ImagedType', vars), currentSessionData.Type = row.ImagedType{1}; end
        
        if ismember('GCaMPInjectionSite', vars), currentSessionData.InjectionSite = row.GCaMPInjectionSite{1}; end
        if ismember('TargetArea', vars), currentSessionData.TargetArea = row.TargetArea{1}; end
        
        % Load Activity Data
        [lapActivity, numLaps, sfi, ~] = loadVRLapSignal(currentSessionData.MouseID, currentSessionData.Session, params.signalToUse);
        [sessionROIData] = loadSessionROIData(sfi); 
        
        if isempty(lapActivity)
            % Still assign to maintain array indexing, then continue
            sessionMetrics(thisSessions) = currentSessionData;
            continue;
        end
        
        % Smooth data 
        if params.applySmoothing
            smoothingKernel = gausswin(15);
            smoothingKernel = smoothingKernel / sum(smoothingKernel);
            lapActivity = smoothActivityTraces(lapActivity, smoothingKernel);
        end
        
        % Tuning curves (using 'omitnan' to handle missing spatial bins)
        % Dimensions: [Cells x Laps x Bins] -> Mean across Laps (Dim 2)
        meanOddLaps = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
        meanEvenLaps = squeeze(mean(lapActivity(:, 2:2:end, :), 2, 'omitnan'));
        meanAllLaps = squeeze(mean(lapActivity, 2, 'omitnan'));
        
        % Fix for single-cell sessions (squeeze might over-flatten)
        if size(lapActivity, 1) == 1
            meanOddLaps = meanOddLaps'; meanEvenLaps = meanEvenLaps'; meanAllLaps = meanAllLaps';
        end
        
        % Modulation Index
        tuningCurveMax = max(meanAllLaps, [], 2, 'omitnan'); 
        tuningCurveMin = min(meanAllLaps, [], 2, 'omitnan');
        modulationIndex = (tuningCurveMax - tuningCurveMin);
        
        currentSessionData.OddMean = meanOddLaps;
        currentSessionData.EvenMean = meanEvenLaps;
        currentSessionData.MeanTuning = meanAllLaps;
        currentSessionData.Modulation = modulationIndex;
        currentSessionData.NumCells = size(lapActivity, 1);
        currentSessionData.NumLaps = numLaps;
        
        % ROI Metadata extraction
        currentSessionData.TypeImaged_ROI = safeGet(sessionROIData, 'typeImaged');
        currentSessionData.TargetArea_ROI = safeGet(sessionROIData, 'targetArea');
        currentSessionData.highlyCorrBoutons = safeGetField(sessionROIData, 'highlyCorrBoutons', 'roisToKeep', []);
  
        % Stability Metrics
        currentSessionData.lapCorr_HalvesRho = safeGetField(sessionROIData, 'lapCorr_Halves', 'rho', []);
        currentSessionData.lapCorr_HalvesP = safeGetField(sessionROIData, 'lapCorr_Halves', 'p', []);
        currentSessionData.lapCorr_HalvesStableThreshold = safeGetField(sessionROIData, 'lapCorr_Halves', 'stableThreshold', []);
        currentSessionData.lapCorr_HalvesStableIdx = safeGetField(sessionROIData, 'lapCorr_Halves', 'stableIdx', []);
        currentSessionData.lapCorr_HalvesFirstHalfLapsCount = safeGetField(sessionROIData, 'lapCorr_Halves', 'firstHalfLapsCount', []);
        currentSessionData.lapCorr_HalvesSecondHalfLapsCount = safeGetField(sessionROIData, 'lapCorr_Halves', 'secondHalfLapsCount', []);
        
        currentSessionData.lapCorr_OddEvenRho = safeGetField(sessionROIData, 'lapCorr_OddEven', 'rho', []);
        currentSessionData.lapCorr_OddEvenP = safeGetField(sessionROIData, 'lapCorr_OddEven', 'p', []);
        currentSessionData.lapCorr_OddEvenStableThreshold = safeGetField(sessionROIData, 'lapCorr_OddEven', 'stableThreshold', []);
        currentSessionData.lapCorr_OddEvenStableIdx = safeGetField(sessionROIData, 'lapCorr_OddEven', 'stableIdx', []);
        
        % Null Distribution Metrics
        if isfield(sessionROIData, 'nullDist_PeakTuningMetric')
            currentSessionData.isSignificantByPeakShuffling = safeGet(sessionROIData.nullDist_PeakTuningMetric.isSignificantByPeakShuffling, params.signalToUse);
            currentSessionData.realPeakPercentileRank = safeGet(sessionROIData.nullDist_PeakTuningMetric.realPeakPercentileRank, params.signalToUse);
            currentSessionData.targetPercentile_Peak = sessionROIData.nullDist_PeakTuningMetric.targetPercentile;
        end
        
        if isfield(sessionROIData, 'nullDist_RangeTuningMetric')
            currentSessionData.isSignificantByRange = safeGet(sessionROIData.nullDist_RangeTuningMetric.isSignificantByRange, params.signalToUse);
            currentSessionData.realRangePercentileRank = safeGet(sessionROIData.nullDist_RangeTuningMetric.realRangePercentileRank, params.signalToUse);
            currentSessionData.targetPercentile_Range = sessionROIData.nullDist_RangeTuningMetric.targetPercentile; 
        end
        
        currentSessionData.ratioVarToTuningVar = safeGetField(sessionROIData, 'tuningCurveVariance', 'ratioVarToTuningVar', []);
        currentSessionData.ratioVarToTuningRange = safeGetField(sessionROIData, 'tuningCurveVariance', 'ratioVarToTuningRange', []);
        
        % Final Assignment
        sessionMetrics(thisSessions) = currentSessionData; 
    end
    sessionMetrics = orderfields(sessionMetrics); 
end

%% --- Helper Functions ---
function val = safeGet(s, field)
    if isstruct(s) && isfield(s, field), val = s.(field); else, val = []; end
end

function val = safeGetField(s, sub, field, default)
    if isfield(s, sub) && isfield(s.(sub), field), val = s.(sub).(field); else, val = default; end
end

function [sessionROIData] = loadSessionROIData(sessionFileInfo)
    sessionROIData = struct(); 
    try
        path = sessionFileInfo.otherSessFilePaths.sessionROIData;
        if exist(path, 'file')
            D = load(path); 
            % Handle if loaded file is wrapped in a 'sessionROIData' variable
            if isfield(D, 'sessionROIData'), sessionROIData = D.sessionROIData; else, sessionROIData = D; end
        end
    catch
        sessionROIData = struct(); 
    end
end

function [lapActivity, numLaps, sfi, errorMessage] = loadVRLapSignal(subjectID, sessionID, signalName)
    lapActivity = []; numLaps = 0; errorMessage = ''; sfi = struct();
    try
        path = findSessionFileInfoFilePath(subjectID, sessionID);
        D = load(path, 'sessionFileInfo'); 
        sfi = D.sessionFileInfo;
        stimNames = string({sfi.stimFiles.name});
        
        % Identify potential baseline corridors (added 'VRCorr')
        baseIndices = find(contains(stimNames, "BaselineCorridor", 'IgnoreCase', true) | ...
                           contains(stimNames, "VRCorr", 'IgnoreCase', true));
        
        if isempty(baseIndices), errorMessage = 'NoBaselineFound'; return; end
        
        % Selection Priority Logic
        selectedIdx = [];
        combinedIdx = baseIndices(contains(stimNames(baseIndices), "CombinedRuns"));
        
        if ~isempty(combinedIdx)
            selectedIdx = combinedIdx(1);
        elseif isscalar(baseIndices)
            selectedIdx = baseIndices;    
        else
            idx00002 = baseIndices(contains(stimNames(baseIndices), "00002"));
            if ~isempty(idx00002), selectedIdx = idx00002(1); else, selectedIdx = baseIndices(1); end
        end
        
        if isempty(selectedIdx), errorMessage = 'SelectionFailed'; return; end
        
        R = load(sfi.stimFiles(selectedIdx).Response);
        if isfield(R.lapPositionActivity, signalName)
            lapActivity = R.lapPositionActivity.(signalName);
            numLaps = size(lapActivity, 2);
        else
            errorMessage = 'SignalNotFound';
        end
    catch
        errorMessage = 'LoadError'; 
    end
end

function smoothedData = smoothActivityTraces(activityData, kernel)
    smoothedData = activityData;
    for c = 1:size(activityData, 1)
        for l = 1:size(activityData, 2)
            trace = squeeze(activityData(c, l, :));
            mask = isnan(trace);
            if all(mask), continue; end
            % Fill NaNs for filtfilt then re-mask
            trace_filled = trace; trace_filled(mask) = 0;
            s = filtfilt(kernel, 1, double(trace_filled));
            s(mask) = NaN;
            smoothedData(c, l, :) = s;
        end
    end
end