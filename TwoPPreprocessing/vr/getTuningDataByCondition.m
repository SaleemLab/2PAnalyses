function sessionMetrics = getTuningDataByCondition(sessionTable, varargin)
    % getTuningDataByCondition: Processes a table of sessions, finds the best
    % stimulus file (Landmark > Baseline), and splits full lap activity by condition.
    
    p = inputParser;
    addRequired(p, 'sessionTable', @istable); 
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar); 
    addParameter(p, 'applySmoothing', true, @islogical);
    parse(p, sessionTable, varargin{:});
    params = p.Results;
    
    numSessions = height(sessionTable);
    sessionMetrics = struct(); 
    
    for thisSess = 1:numSessions
        row = sessionTable(thisSess,:);
        currentSessionData = struct(); 
        
        currentSessionData.MouseID = row.MouseID{1};
        currentSessionData.Day = row.DayOfExperience;
        currentSessionData.Session = char(row.Session);
        
        % Load & Split
        [condActivity, sfi, errorMessage] = loadVRLapSignalByCondition(currentSessionData.MouseID, currentSessionData.Session, params.signalToUse); 
        
        if isempty(fieldnames(condActivity))
            warning('Skipping %s: %s', currentSessionData.Session, errorMessage);
            continue;
        end
        
        condNames = fieldnames(condActivity);
        for c = 1:length(condNames)
            name = condNames{c};
            data = condActivity.(name); 
            
            if params.applySmoothing
                smoothingKernel = gausswin(15); 
                smoothingKernel = smoothingKernel / sum(smoothingKernel);
                % --- THIS IS THE FUNCTION THAT WAS MISSING ---
                data = smoothActivityTraces(data, smoothingKernel);
            end
            
            currentSessionData.ConditionData.(name).LapActivity = data;
            currentSessionData.ConditionData.(name).NumLaps = size(data, 2);
            currentSessionData.ConditionData.(name).MeanTuning = squeeze(mean(data, 2, 'omitnan'));
        end
        
        [sessionROIData] = loadSessionROIData(sfi); 
        currentSessionData.NumCells = size(condActivity.(condNames{1}), 1);
        currentSessionData.TargetArea_ROI = safeGet(sessionROIData, 'targetArea');
        
        if isempty(fieldnames(sessionMetrics))
            sessionMetrics = currentSessionData;
        else
            sessionMetrics(end+1) = currentSessionData;
        end
    end
end

%% --- HELPER: SMOOTHING ---
function smoothedData = smoothActivityTraces(activityData, filterKernel)
    smoothedData = activityData;
    % activityData is [Cells x Laps x Position]
    for cellIdx = 1:size(activityData, 1)
        for lapIdx = 1:size(activityData, 2)
            trace = squeeze(activityData(cellIdx, lapIdx, :));
            nanMask = isnan(trace);
            if all(nanMask), continue; end
            
            trace(nanMask) = 0;
            filteredTrace = filtfilt(filterKernel, 1, trace);
            filteredTrace(nanMask) = NaN;
            smoothedData(cellIdx, lapIdx, :) = filteredTrace;
        end
    end
end

%% --- HELPER: LOADING ---
function [condActivity, sfi, errorMessage] = loadVRLapSignalByCondition(subjectID, sessionID, signalName)
    condActivity = struct(); errorMessage = ''; sfi = struct();
    try
        sessionInfoPath = findSessionFileInfoFilePath(subjectID, sessionID);
        if ~isfile(sessionInfoPath), errorMessage = 'InfoMissing'; return; end
        D = load(sessionInfoPath, 'sessionFileInfo');
        sfi = D.sessionFileInfo; 
        
        stimNames = string({sfi.stimFiles.name});
        isLandmark = contains(stimNames, "LandManipCorridor");
        isBaseline = contains(stimNames, "BaselineCorridor");
        isCombined = contains(stimNames, "CombinedRuns");
        
        targetIdx = find(isLandmark & isCombined, 1);
        if isempty(targetIdx), targetIdx = find(isBaseline & isCombined, 1); end
        if isempty(targetIdx), targetIdx = find(isLandmark | isBaseline, 1, 'last'); end

        if isempty(targetIdx), errorMessage = 'NoValidStimFound'; return; end
        
        R = load(sfi.stimFiles(targetIdx).Response);
        fullActivity = R.lapPositionActivity.(signalName); 
        
        if isfield(R, 'trialIndicesByCondition')
            fn = fieldnames(R.trialIndicesByCondition);
            for i = 1:length(fn)
                trialIdx = R.trialIndicesByCondition.(fn{i});
                validLaps = trialIdx(trialIdx > 0 & trialIdx <= size(fullActivity, 2));
                if ~isempty(validLaps)
                    condActivity.(fn{i}) = fullActivity(:, validLaps, :);
                end
            end
        else
            condActivity.Default = fullActivity;
        end
    catch ME
        errorMessage = ME.message;
    end
end

%% --- UTILITIES ---
function [sessionROIData] = loadSessionROIData(sessionFileInfo)
    sessionROIData = struct(); 
    try
        path = sessionFileInfo.otherSessFilePaths.sessionROIData;
        if exist(path, 'file'), sessionROIData = load(path); end
    catch
        sessionROIData = struct(); 
    end
end

function val = safeGet(s, field)
    if isstruct(s) && isfield(s, field), val = s.(field); else, val = []; end
end



% % % Create a blank template based on the last entry
% 
% [sessionFileInfo.stimFiles.Response] = deal([]);
% 
% 
% newEntry = sessionFileInfo.stimFiles(end);
% 
% 
% newEntry.name = 'M26003_LandManipCorridor_20260325_CombinedRuns';
% newEntry.bonsai_filepaths = {}; 
% newEntry.eyetracking_filepaths = {}; 
% newEntry.tif_filepaths = {};
% newEntry.TwoPMetaData = '';
% newEntry.processedPeripheralData = '';
% newEntry.mergedBonsai2PSuite2pData = '';
% newEntry.Response = 'Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat';
% 
% 
% sessionFileInfo.stimFiles(end+1) = newEntry;
% 
% 
% save('\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_sessionFileInfo.mat', 'sessionFileInfo');