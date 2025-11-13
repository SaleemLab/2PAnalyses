function [lapActivity, numLaps, msg] = loadVRData(mouseID, sessionRow, signalToUse)
% loadVRData Loads lap activity for a specific session row.
% Returns empty activity and a message if something fails.

    lapActivity = []; numLaps = 0; msg = '';
    
    % Handle potential duplicate grid entries (take first)
    if height(sessionRow) > 1
        warning('Multiple entries found for %s. Using first. Please check and possibly combined Response files', char(sessionRow.Session(1)));
        sessionRow = sessionRow(1, :);
    end
    sessionStr = char(sessionRow.Session);

    try
        % 1. Locate Session Info
        infoPath = findSessionFileInfoFilePath(mouseID, sessionStr);
        if ~isfile(infoPath)
             msg = 'Info Missing'; return; 
        end
        loadedInfo = load(infoPath, 'sessionFileInfo');
        sfi = loadedInfo.sessionFileInfo;

        % 2. Find correct VRCorr file (Combined preferred, then last VR only)
        stimNames = string({sfi.stimFiles.name});
        idx = find(contains(stimNames, "VRCorr") & contains(stimNames, "CombinedRuns"), 1);
        if isempty(idx)
             idx = find(contains(stimNames, "VRCorr") & ~contains(stimNames, "CombinedRuns"), 1, 'last');
        end

        if isempty(idx) || ~isfield(sfi.stimFiles(idx), 'Response') || isempty(sfi.stimFiles(idx).Response)
             msg = 'No VRCorr/Response'; return;
        end

        % 3. Load Response
        respPath = sfi.stimFiles(idx).Response;
        if ~isfile(respPath)
            msg = 'Resp File Missing'; return; 
        end
        loadedResp = load(respPath, 'response');
        
        if ~isfield(loadedResp, 'response') || ...
           ~isfield(loadedResp.response, 'lapPositionActivity') || ...
           ~isfield(loadedResp.response.lapPositionActivity, signalToUse)
             msg = 'Signal Missing'; return;
        end
        
        % 4. Extract Data
        lapActivity = loadedResp.response.lapPositionActivity.(signalToUse);
        numLaps = size(lapActivity, 2);
        
        if numLaps < 2
             lapActivity = []; msg = '<2 Laps';
        end

    catch ME
        % msg = ME.message; % Use this for detailed debugging if needed
        msg = 'Load Error'; 
    end
end