function sessionMetrics = getTuningData(sessionTable, varargin)
    % Original: sessionData = getTuningData(filteredTable, varargin)
    %
    % getTuningData: Loads VR data, computes MEAN tuning curves 
    % and Modulation Index for each session in the input table.
    

    p = inputParser;
    addRequired(p, 'sessionTable', @istable); % Renamed 'filteredTable'
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    parse(p, sessionTable, varargin{:});
    params = p.Results;
    

    numSessions = height(sessionTable);

    sessionMetrics = struct('MouseID', {}, 'Day', {}, 'Session', {}, 'Type', {}, ...
                         'InjectionSite', {}, 'TargetArea', {}, 'OddMean', {}, 'EvenMean', {}, ...
                         'MeanTuning', {}, 'Modulation', {}, 'NumCells', {}, 'NumLaps', {});
    
    wb = waitbar(0, 'Loading and processing sessions...');
    
    % Session Loop
    for thisSessions = 1:numSessions
        waitbar(thisSessions/numSessions, wb, sprintf('Processing session %d/%d...', thisSessions, numSessions));
        
        row = sessionTable(thisSessions,:);
        clear currentSessionData; % Renamed 'sData'
        
        % Copy Metadata
        currentSessionData.MouseID = row.MouseID{1};
        currentSessionData.Day = row.DayOfExperience;
        currentSessionData.Session = char(row.Session);
        
        % Handle different possible column names for 'Type'; fixed now.. 
        if ismember('TypeImaged', row.Properties.VariableNames)
             currentSessionData.Type = row.TypeImaged{1};
        elseif ismember('ImagedType', row.Properties.VariableNames)
             currentSessionData.Type = row.ImagedType{1};
        else
             currentSessionData.Type = 'Unknown';
        end
        
        % Handle optional column for 'Area'
        if ismember('GCaMPInjectionSite', row.Properties.VariableNames)
            currentSessionData.InjectionSite = row.GCaMPInjectionSite{1};
            currentSessionData.TargetArea = row.TargetArea{1};

        else
             currentSessionData.InjectionSite = '';
             currentSessionData.TargetArea = '';
        end
        
        % Load Activity Data
        [lapActivity, numLaps, errorMessage] = loadVRLapSignal(currentSessionData.MouseID, currentSessionData.Session, params.signalToUse);
        
        if isempty(lapActivity)
            warning('Skipping %s Session %s Day %d: %s', currentSessionData.MouseID, currentSessionData.Session, currentSessionData.Day, errorMessage);
            continue;
        end
        
        % Smooth data (CHECK) 
        if params.applySmoothing
            smoothingKernel = gausswin(9); 
            smoothingKernel = smoothingKernel / sum(smoothingKernel);
            
            lapActivity = smoothActivityTraces(lapActivity, smoothingKernel);
        end
        
        % 4. Tuning Curves
        % (Mean activity across position bins, separated by lap type)
        meanOddLaps = squeeze(mean(lapActivity(:, 1:2:end, :), 2, 'omitnan'));
        meanEvenLaps = squeeze(mean(lapActivity(:, 2:2:end, :), 2, 'omitnan'));
        meanAllLaps = squeeze(mean(lapActivity, 2, 'omitnan'));
        
        % Handle case with only one cell, where squeeze might remove too
        % much; Exclude these sessions instead?
        % if ismatrix(lapActivity) && size(lapActivity,1) == 1
        %      meanOddLaps = reshape(meanOddLaps, 1, []);
        %      meanEvenLaps = reshape(meanEvenLaps, 1, []);
        %      meanAllLaps = reshape(meanAllLaps, 1, []);
        % end
        
        % Modulation (acoss space) Index
        % Modulation Index: (Max - Min) / Mean
        tuningCurveMax = max(meanAllLaps, 2, 'omitnan'); % Find max within each row i.e., position bin 
        tuningCurveMin = min(meanAllLaps, 2, 'omitnan');
        tuningCurveMean = mean(meanAllLaps, 2, 'omitnan');
        
        modulationIndex = (tuningCurveMax - tuningCurveMin) ./ tuningCurveMean;
        
        % % Clean invalid MI values (e.g., from non-active cells where mean=0)
        % invalidMIFilter = (tuningCurveMean == 0 | isnan(tuningCurveMean) | isinf(modulationIndex));
        % modulationIndex(invalidMIFilter) = NaN;
        
        % Store Data for this Session
        currentSessionData.OddMean = meanOddLaps;
        currentSessionData.EvenMean = meanEvenLaps;
        currentSessionData.MeanTuning = meanAllLaps;
        currentSessionData.Modulation = modulationIndex;
        currentSessionData.NumCells = size(lapActivity, 1);
        currentSessionData.NumLaps = numLaps;
        
        % Append the struct for this session to the main array
        sessionMetrics(end+1) = orderfields(currentSessionData, sessionMetrics); 
    end
    
    close(wb);
    fprintf('Data loaded for %d sessions (with Modulation metrics).\n', length(sessionMetrics));
end

%% local helper function 
function [lapActivity, numLaps, errorMessage] = loadVRLapSignal(subjectID, sessionID, signalName)

    lapActivity = [];
    numLaps = 0;
    errorMessage = '';
    
    try
        % Find the path to the session's info file
        sessionInfoPath = findSessionFileInfoFilePath(subjectID, sessionID);
        if ~isfile(sessionInfoPath)
            errorMessage = 'InfoMissing';
            return;
        end
        
        % Load the session file info
        D = load(sessionInfoPath, 'sessionFileInfo');
        sfi = D.sessionFileInfo;
        
        % Find the correct VR stimulus file
        stimFileNames = string({sfi.stimFiles.name});
        vrFileIndex = find(contains(stimFileNames, "VRCorr") & contains(stimFileNames, "CombinedRuns"), 1);
        
        if isempty(vrFileIndex) % Fallback: Find last run if no 'CombinedRuns'
            vrFileIndex = find(contains(stimFileNames, "VRCorr") & ~contains(stimFileNames, "CombinedRuns"), 1, 'last'); 
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
        R = load(responseFilePath, 'response');
        
        % Check if the signal exists in the lap data structure
        if ~isfield(R.response, 'lapPositionActivity') || ~isfield(R.response.lapPositionActivity, signalName)
            errorMessage = 'SignalMissing'; 
            return;
        end
        
        % Extract the requested signal data
        lapActivity = R.response.lapPositionActivity.(signalName);
        numLaps = size(lapActivity, 2);
        
        % Ensure there is more than one lap to be useful
        if numLaps < 2
            lapActivity = [];
            errorMessage = '<2Laps';
            return; 
        end
        
    catch
        errorMessage = 'LoadError';
    end
end

% Move to utilities!! 
%% 
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
            
            % Temporarily set NaNs to 0 for filtering
            trace(nanMask) = 0;
            
            % Apply zero-phase digital filter
            filteredTrace = filtfilt(filterKernel, 1, trace);
            
            % Restore NaN values
            filteredTrace(nanMask) = NaN;
            
            % Assign the smoothed trace back to the output array
            smoothedData(cellIdx, lapIdx, :) = filteredTrace;
        end
    end
end