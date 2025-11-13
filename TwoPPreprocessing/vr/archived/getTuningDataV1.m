function sessionData = getTuningDataV1(filteredTable, varargin)
% getTuningData Loads VR data and computes MEAN tuning curves.
% Now uses MEAN across laps instead of Median.

    p = inputParser;
    addRequired(p, 'filteredTable', @istable);
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar);
    addParameter(p, 'applySmoothing', true, @islogical);
    parse(p, filteredTable, varargin{:});
    params = p.Results;

    numSessions = height(filteredTable);
    
    % Updated pre-allocation with 'Mean' field names 
    sessionData = struct('MouseID', {}, 'Day', {}, 'Session', {}, 'Type', {}, ...
                         'InjectionSite', {}, 'TargetArea', {},  'OddMean', {}, 'EvenMean', {}, ...
                         'MeanTuning', {}, 'NumCells', {}, 'NumLaps', {});

    wb = waitbar(0, 'Loading and processing sessions...');
    
    for i = 1:numSessions
        waitbar(i/numSessions, wb, sprintf('Processing session %d/%d...', i, numSessions));
        row = filteredTable(i,:);
        clear sData; % Vital: clear previous iteration's data
        
        % --- Metadata ---
        sData.MouseID = row.MouseID{1};
        sData.Day = row.DayOfExperience;
        sData.Session = char(row.Session);
        
        if ismember('TypeImaged', row.Properties.VariableNames)
             sData.Type = row.TypeImaged{1};
        else
             sData.Type = 'Unknown';
        end
        
        if ismember('GCaMPInjectionSite', row.Properties.VariableNames)
            sData.InjectionSite = row.GCaMPInjectionSite{1};
        else
             sData.InjectionSite = '';
        end

        if ismember('TargetArea', row.Properties.VariableNames)
            sData.TargetArea = row.TargetArea{1};
        else
            sData.TargetArea = '';
        end

        % --- Load ---
        [lapActivity, numLaps, err] = internalLoadVR(sData.MouseID, sData.Session, params.signalToUse);
        if isempty(lapActivity)
            warning('Skipping %s Day %d: %s', sData.MouseID, sData.Day, err);
            continue;
        end

        % --- Process ---
        if params.applySmoothing
            w = gausswin(6); w = w/sum(w);
            lapActivity = smoothTrace(lapActivity, w);
        end

        oddLaps = lapActivity(:, 1:2:end, :);
        evenLaps = lapActivity(:, 2:2:end, :);
        
        % --- CHANGE IS HERE: Changed median() to mean() ---
        meanOdd = squeeze(mean(oddLaps, 2, 'omitnan'));
        meanEven = squeeze(mean(evenLaps, 2, 'omitnan'));
        meanAll = squeeze(mean(lapActivity, 2, 'omitnan'));

        if ismatrix(lapActivity) && size(lapActivity,1) == 1
             meanOdd = reshape(meanOdd, 1, []);
             meanEven = reshape(meanEven, 1, []);
             meanAll = reshape(meanAll, 1, []);
        end

        % --- Updated Field Names ---
        sData.OddMean = meanOdd;
        sData.EvenMean = meanEven;
        sData.MeanTuning = meanAll;
        sData.NumCells = size(lapActivity, 1);
        sData.NumLaps = numLaps;

        sessionData(end+1) = orderfields(sData, sessionData); 
    end
    close(wb);
    fprintf('Data loaded for %d sessions (using MEAN across laps).\n', length(sessionData));
end

%% Local Helper Functions
function [lapAct, nLaps, msg] = internalLoadVR(mID, sStr, sig)
    lapAct = []; nLaps = 0; msg = '';
    try
        infoP = findSessionFileInfoFilePath(mID, sStr);
        if ~isfile(infoP), msg='InfoMissing';
            return; 
        end
        D = load(infoP, 'sessionFileInfo'); 
        sfi = D.sessionFileInfo;
        stims = string({sfi.stimFiles.name});
        idx = find(contains(stims,"VRCorr") & contains(stims,"CombinedRuns"), 1);
        if isempty(idx), 
            idx = find(contains(stims,"VRCorr") & ~contains(stims,"CombinedRuns"),1,'last'); 
        end
        if isempty(idx), msg='NoVRCorr'; 
            return; 
        end
        respP = sfi.stimFiles(idx).Response;
        if isempty(respP) || ~isfile(respP), 
            msg='RespMissing'; 
            return; 
        end
        R = load(respP, 'response');
        if ~isfield(R.response, 'lapPositionActivity') || ~isfield(R.response.lapPositionActivity, sig)
            msg='SignalMissing'; 
            return;
        end
        lapAct = R.response.lapPositionActivity.(sig);
        nLaps = size(lapAct, 2);
        if nLaps < 2, lapAct=[]; msg='<2Laps'; return; end
    catch
        msg='LoadError';
    end
end

function smoothed = smoothTrace(act, w)
    smoothed = act;
    for c=1:size(act,1)
        for l=1:size(act,2)
            t = squeeze(act(c,l,:));
            mask = isnan(t);
            if all(mask), continue; end
            t(mask) = 0;
            f = filtfilt(w,1,t);
            f(mask) = NaN;
            smoothed(c,l,:) = f;
        end
    end
end