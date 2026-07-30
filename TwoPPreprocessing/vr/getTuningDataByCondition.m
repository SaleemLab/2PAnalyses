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
        
        [condActivity, sfi, movementVisualGain, SMI, SMR, flaggedLaps, ~] = loadVRLapSignalByCondition(currentSessionData.MouseID, currentSessionData.Session, params.signalToUse); 
        if isempty(fieldnames(condActivity)), continue; end
        % new: saves the gain 
        
        
        currentSessionData.movementVisualGain=movementVisualGain; 

        currentSessionData.SMI = SMI;
        currentSessionData.SMR = SMR;

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
        currentSessionData.cvExpVar = safeGetField(sessionROIData, 'crossValExpVar', params.signalToUse); 

        

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

function [condActivity, sfi, movementVisualGain, SMI, SMR, flaggedLaps, errorMessage] = loadVRLapSignalByCondition(subjectID, sessionID, signalName)
    condActivity = struct(); errorMessage = ''; sfi = struct();
    % FIX: initialize all outputs up front so a caught error (or a
    % missing field) never leaves an output argument unassigned, which
    % would otherwise throw an uncatchable "output not assigned" error
    % on function exit.
    SMI = []; SMR = []; movementVisualGain = []; flaggedLaps = [];
    try
        sessionInfoPath = findSessionFileInfoFilePath(subjectID, sessionID);
        D = load(sessionInfoPath, 'sessionFileInfo');
        sfi = D.sessionFileInfo; 
        stimNames = string({sfi.stimFiles.name});
        
%         targetIdx = find(contains(stimNames, "LandManipCorridor") & contains(stimNames, "CombinedRuns"), 1);
%         
%         if isempty(targetIdx)
%             targetIdx = find(contains(stimNames, "LandManipCorridor"), 1);
%         end
%         
%         if isempty(targetIdx)
%             targetIdx = find(contains(stimNames, "BaselineCorridor") & contains(stimNames, "CombinedRuns"), 1);
%         end
%         
%         if isempty(targetIdx)
%             targetIdx = find(contains(stimNames, "BaselineCorridor"), 1);
%         end
% 
%         if isempty(targetIdx)
%             error('No valid LandManip or Baseline corridor found.');
%         end

        % select combined runs if present 
        targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);

        if isempty(targetIdx)
            allCorridorIdx = find(contains(stimNames, "Corridor"));

            if isscalar(allCorridorIdx)
           
                targetIdx = allCorridorIdx;
            elseif length(allCorridorIdx) > 1
                
                targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "00002"), 1);
            end


        end
        
        responseFilePath = sfi.stimFiles(targetIdx).Response;
    
        % FIX: 'SMR_metric' -> 'SMR_Metrics' (must match the variable name
        % exactly as saved by computeSpatialModulationRatio, since load()
        % silently skips names that don't match -- it does not error).
        response = load(responseFilePath, 'lapPositionActivity', 'movementVisualGain', 'trialIndicesByCondition', 'SMI_Metrics', 'flaggedLaps', 'stimName', 'SMR_Metrics');
        
        
        % response = load(sfi.stimFiles(targetIdx).Response);        
        if ~isfield(response.lapPositionActivity, signalName)
            error('Signal %s missing from data.', signalName);
         end
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

        % FIX: SMI_Metrics.(signalName) is a struct with fields
        % {SMI, Rp, Rn, RpBin, RnBin, ...}. Index into .SMI to get the
        % plain per-ROI index vector, matching how SMI is used downstream
        % (currentSessionData.SMI = SMI;) as a plain numeric array.
        if isfield(response, 'SMI_Metrics') && ~isempty(response.SMI_Metrics) && isfield(response.SMI_Metrics, signalName)
           SMI = response.SMI_Metrics.(signalName);       
        end

        % FIX: same issue as SMI above -- index into .SMR for the plain
        % ratio vector rather than assigning the whole metrics struct.
        if isfield(response, 'SMR_Metrics') && ~isempty(response.SMR_Metrics) && isfield(response.SMR_Metrics, signalName)
           SMR = response.SMR_Metrics.(signalName);       
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

% gemini functions to load fields 
% s is the sturct; sub is the sub-struct and and then the field name 
function val = safeGet(s, field)
    if isstruct(s) && isfield(s, field)
        val = s.(field); 
    else
        val = []; 
    end
end

function val = safeGetField(s, sub, field, default)
    if isstruct(s) && isfield(s, sub) && isfield(s.(sub), field)
        val = s.(sub).(field); 
    else
        val = default; 
    end
end