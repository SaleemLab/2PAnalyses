function sessionMetrics = getTuningDataByCondition(sessionTable, varargin)
% load data based on trial conditions; 
% saves various parameters to filter rois 
% temporarity computes fraction of variance explained [will add this to
% each sessions soon] 

    p = inputParser;
    addRequired(p, 'sessionTable', @istable); 
    addParameter(p, 'signalToUse', 'dFFNeuropilCorrected', @ischar); 
    % addParameter(p, 'applySmoothing', true, @islogical);
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
        currentSessionData.TypeImaged = char(row.TypeImaged);
        if contains(currentSessionData.TypeImaged, 'Somas')
            currentSessionData.TargetArea = 'V1';
        end
        fprintf('Processing Mouse %s Session %s \n', currentSessionData.MouseID, currentSessionData.Session)
        
        [condActivity, sfi, movementVisualGain, SMI, flaggedLaps, ~] = loadVRLapSignalByCondition(currentSessionData.MouseID, currentSessionData.Session, params.signalToUse); 
        if isempty(fieldnames(condActivity)), continue; end
        % new: saves the gain 
        
        
        currentSessionData.movementVisualGain=movementVisualGain; 

        currentSessionData.SMI = SMI;

        currentSessionData.flaggedLaps = flaggedLaps; 
        
        if length(currentSessionData.movementVisualGain)>1
            if isscalar(unique(currentSessionData.movementVisualGain))
                currentSessionData.movementVisualGain = movementVisualGain(1);
            end
        end 
        [sessionROIData] = loadSessionROIData(sfi); 
        
        currentSessionData.Rho_OddEven = safeGetField(sessionROIData, 'lapCorr_OddEven', 'rho', []);
        currentSessionData.Rho_Halves  = safeGetField(sessionROIData, 'lapCorr_Halves', 'rho', []);
        currentSessionData.TuningVarMetrics = safeGet(sessionROIData, 'tuningCurveVariance');
        currentSessionData.NullDist_Range = safeGet(sessionROIData, 'nullDist_RangeTuningMetric');
        currentSessionData.NullDist_Peak = safeGet(sessionROIData, 'nullDist_PeakTuningMetric');
        currentSessionData.cvExpVar = safeGetField(sessionROIData, 'crossValExpVar', 'cvExpVar'); 
        currentSessionData.cvExpvar_nullDist_pVales = safeGetField(sessionROIData, 'crossValExpVar', 'pValues');
        

        if  strcmpi(currentSessionData.TypeImaged, 'Somas')
            currentSessionData.sparseNoiseMatrix = safeGet(sessionROIData, 'sparseNoiseRF');
        end
      
     
        if  strcmpi(currentSessionData.TypeImaged, 'Boutons')
            currentSessionData.uniqueBoutonIdx = safeGetField(sessionROIData, 'highlyCorrBoutons','roisToKeep');
            currentSessionData.highCorrBoutonIdx = safeGetField(sessionROIData, 'highlyCorrBoutons','roisToDiscard');
            % this seems to be missing for somas; have imported above from
            % data table 
            currentSessionData.TargetArea = safeGet(sessionROIData, 'targetArea');
        end

        condNames = fieldnames(condActivity);
        for c = 1:length(condNames)
            name = condNames{c};
            data = condActivity.(name); 
            % if params.applySmoothing
            %     data = smoothActivityTraces(data, gausswin(15)/sum(gausswin(15)));
            % end
            currentSessionData.ConditionData.(name).LapActivity = data;
            currentSessionData.ConditionData.(name).NumLaps = size(data, 2);
        end
        
        if isempty(fieldnames(sessionMetrics))
            sessionMetrics = currentSessionData;
        else
            sessionMetrics(end+1) = currentSessionData;
        end
    end
end

% smooth before plotting; skipping this here.. 
% function smoothedData = smoothActivityTraces(activityData, filterKernel)
%     smoothedData = activityData;
%     for cellIdx = 1:size(activityData, 1)
%         for lapIdx = 1:size(activityData, 2)
%             trace = squeeze(activityData(cellIdx, lapIdx, :));
%             nanMask = isnan(trace);
%             if all(nanMask), continue; end
%             trace(nanMask) = 0;
%             filteredTrace = filtfilt(filterKernel, 1, double(trace));
%             filteredTrace(nanMask) = NaN;
%             smoothedData(cellIdx, lapIdx, :) = filteredTrace;
%         end
%     end
% end

function [condActivity, sfi, movementVisualGain, SMI, flaggedLaps, errorMessage] = loadVRLapSignalByCondition(subjectID, sessionID, signalName)
    condActivity = struct(); errorMessage = ''; sfi = struct();
    try
        sessionInfoPath = findSessionFileInfoFilePath(subjectID, sessionID);
        D = load(sessionInfoPath, 'sessionFileInfo');
        sfi = D.sessionFileInfo; 
        stimNames = string({sfi.stimFiles.name});
        
        targetIdx = find(contains(stimNames, "LandManipCorridor") & contains(stimNames, "CombinedRuns"), 1);
        
        if isempty(targetIdx)
            targetIdx = find(contains(stimNames, "LandManipCorridor"), 1);
        end
        
        if isempty(targetIdx)
            targetIdx = find(contains(stimNames, "BaselineCorridor") & contains(stimNames, "CombinedRuns"), 1);
        end
        
        if isempty(targetIdx)
            targetIdx = find(contains(stimNames, "BaselineCorridor"), 1);
        end

        if isempty(targetIdx)
            error('No valid LandManip or Baseline corridor found.');
        end
        
        responseFilePath = sfi.stimFiles(targetIdx).Response;
    
        response = load(responseFilePath, 'lapPositionActivity', 'movementVisualGain', 'trialIndicesByCondition', 'SMI_Metrics', 'flaggedLaps', 'stimName');
        
        
        % response = load(sfi.stimFiles(targetIdx).Response);        
        if ~isfield(response.lapPositionActivity, signalName)
            error('Signal %s missing from data.', signalName);
         end
        movementVisualGain = [];
        if isfield(response, 'movementVisualGain')
            movementVisualGain = response.movementVisualGain; 
        end 
        fullActivity = response.lapPositionActivity.(signalName); 

        if isfield(response, 'trialIndicesByCondition')
            fn = fieldnames(response.trialIndicesByCondition);
            for i = 1:length(fn)
                trialIdx = response.trialIndicesByCondition.(fn{i});
                validLaps = trialIdx(trialIdx > 0 & trialIdx <= size(fullActivity, 2));
                if ~isempty(validLaps)
                    condActivity.(fn{i}) = fullActivity(:, validLaps, :); 
                end
            end
        else
            condActivity.Default = fullActivity;
        end

        if isfield(response, 'SMI_Metrics') && ~isempty(response.SMI_Metrics)
           SMI = response.SMI_Metrics.(signalName);       
        end

        if isfield(response, 'flaggedLaps') %&& ~isempty(response.flaggedLaps)
           flaggedLaps = response.flaggedLaps;       
        end
        
    catch ME
        errorMessage = ME.message;
    end
end

function [sessionROIData] = loadSessionROIData(sessionFileInfo)
    sessionROIData = struct(); 
    try
        path = sessionFileInfo.otherSessFilePaths.sessionROIData;
        if exist(path, 'file')
            D = load(path); 
            if isfield(D, 'sessionROIData'), sessionROIData = D.sessionROIData; else, sessionROIData = D; end
        end
    catch
        sessionROIData = struct(); 
    end
end

function val = safeGet(s, field)
    if isstruct(s) && isfield(s, field), val = s.(field); else, val = []; end
end

function val = safeGetField(s, sub, field, default)
    if isstruct(s) && isfield(s, sub) && isfield(s.(sub), field), val = s.(sub).(field); else, val = default; end
end