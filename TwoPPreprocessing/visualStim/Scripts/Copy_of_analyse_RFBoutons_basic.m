% This script loads the rf mapping stimulus for all recordings with unique
% fovs - across 3 mice; these are from day 1 of experience
 % Gaussian-fit quality cutoff for "well-fit" RF (used by gaussianR2 method)
r2Thresh        = 0.1;           
selectionMethod = 'isResponsive'; 
% 'isResponsive' or 'gaussianR2' -- which criterion feeds
% Panel A example selection and Panels B/C/D populations.
% 'isResponsive' = analyseRFMapping.m's ANOVA(p<0.05) AND
%  prefVal>blankMean+2*blankStd criterion (currently 53 boutons).

crossValVarName = 'crossValExpVar';   

% identify the correct name.
signalName = 'dFFNeuropilCorrected';

% this will always be neu for the boutons but just to be safe 

%% sessions to load 

pairs=struct;
% fov sweep sessions: 16 grids (including blanks) 
pairs.M25132 = {...
%                 '20260214A','20260214B','20260214C','20260214D', ...
%                 '20260216A','20260216B','20260216C', ...
                '20260219','20260223','20260226','20260228','20260303','20260313','20260306'};
pairs.M25133 = {...
%                 '20260216A','20260216B','20260216C','20260216D', ...
                 '20260219','20260223','20260221'};
pairs.M26003 = {...
%                 '20260307A','20260307B', ...
%                 '20260313A','20260313B','20260313C', ...
                '20260316','20260322','20260324','20260325'};


filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, 'HasStimulus', {'RFMapping', 'BaselinCorridor', 'LandManipCorridor'});
allMice    = filteredTable.MouseID;
uniqueMice = unique(allMice, 'stable');
%% initialize pooled outputs
allRFMapping        = [];   % concatenated RFMapping struct array across ALL boutons/sessions
RFMappingMetadata    = [];  % metadata from the first session that has it — used as the shared reference
sessionLabels        = {};  % track which mouse/session each bouton chunk came from (for debugging/QC)
allMedianCrossValR2  = [];  % cross-validated explained variance from the vr per bouton
allMeanCrossValR2    = [];  % 
allShuffPVal         = [];  % 
haveWarnedMissingVar  = false; % only print the "available variables" note once
%% loop over mice
for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));
    fprintf('MOUSE: %s | %d sessions\n', thisMouse, length(mouseSessIdx));
    %% loop over sessions for this mouse
    for iSess = 1:length(mouseSessIdx)
        tableRow        = filteredTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        thisDay         = tableRow.DayOfExperience;
        thisArea        = tableRow.TargetArea;
        fprintf('\n--- Session: %s | Day %d | Area: %s ---\n', ...
            thisSessionName, thisDay, char(thisArea));
        %% load sessionFileInfo
        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath)
            warning('sfi missing for %s — skipping.', thisSessionName);
            continue;
        end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        % find all RFMapping runs
        RFMapIdx = find(contains(stimNames, 'RFMapping'));
        if isempty(RFMapIdx)

            
            warning('No RFMapping files for %s — skipping.', thisSessionName);

            continue;
        end
        %% load RFMapping + metadata for this session (from sessionROIData)
        try
            sessionRFData = load(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
                'RFMapping', 'RFMappingMetadata');
        catch ME
            warning('  Could not load RFMapping for %s: %s', thisSessionName, ME.message);
            continue;
        end
        if ~isfield(sessionRFData, 'RFMapping') || isempty(sessionRFData.RFMapping)
            warning('  RFMapping empty/missing for %s — skipping.', thisSessionName);
            continue;
        end
        thisRFMapping = sessionRFData.RFMapping;

        try 
            responseData = load(sessionFileInfo.stimFiles(RFMapIdx).Response, 'response');
        catch ME
            warning('  Could not load Response for %s: %s', thisSessionName, ME.message);
            continue;
        end
        thisRFMapping = thisRFMapping(:); 

        % Build the RAW (pre-baseline-subtraction) aligned trial matrix directly from the
        % Response file's psthData, rather than relying on sessionROIData's baselineSubtracted
        % field (which analyseRFMapping.m computes by subtracting a SECOND, per-trial
        % pre-stimulus baseline on top of whatever global dF/F normalization was already
        % applied). This attaches rawAlignedTrials/rawAlignedTrialsBlank to thisRFMapping
        % BEFORE concatenation into allRFMapping, so it flows through the rest of the pooling
        % pipeline for free -- no need to re-run or re-save analyseRFMapping.m for any session.
        response = responseData.response;
        psthData = response.psthData;
        stimVs   = vertcat(psthData.stimValue);
        nROI_thisSess = size(psthData(1).alignedResponses, 1);

        % Separate grid vs blank positions -- blank = stimValue [200 0], matching analyseRFMapping.m
        blankIdxPos = find(stimVs(:,1) == 200 & stimVs(:,2) == 0, 1);
        gridMaskPos = stimVs(:,1) ~= 200;
        gridPSTH_thisSess  = psthData(gridMaskPos);
        gridStim_thisSess  = stimVs(gridMaskPos, :);
        blankPSTH_thisSess = psthData(blankIdxPos);

        % Local grid arrangement for THIS session: use sessionRFData.RFMappingMetadata (already
        % loaded above), NOT a freshly-recomputed uAz/uEl from the raw response file. Different
        % mice/sessions can have the SAME grid SHAPE (e.g. 3x5) at DIFFERENT absolute
        % azimuth/elevation offsets (depending on cranial window position) -- recomputing
        % uAz/uEl independently here risked matching sizes while silently misordering which
        % column/row corresponds to which real position, scrambling azimuth/elevation labels
        % across sessions. Using this session's own already-verified metadata keeps it
        % self-consistent with baselineSubtracted/meanGridResponse for the SAME session; the
        % grid-consistency check further below still catches any session that doesn't match
        % the shared reference grid before pooling.
        uAz_thisSess      = sessionRFData.RFMappingMetadata.uAz;
        uEl_thisSess_plot = sessionRFData.RFMappingMetadata.uEl;
        nAz_thisSess = numel(uAz_thisSess);
        nEl_thisSess = numel(uEl_thisSess_plot);

        if numel(thisRFMapping) ~= nROI_thisSess
            warning(['  ROI count mismatch between sessionROIData RFMapping (%d) and Response ' ...
                     'psthData (%d) for %s -- skipping raw-trial attachment for this session.'], ...
                     numel(thisRFMapping), nROI_thisSess, thisSessionName);
        else
            for iROI = 1:nROI_thisSess
                rawTrialMatrix = cell(nEl_thisSess, nAz_thisSess);
                for thisPos = 1:numel(gridPSTH_thisSess)
                    allTrialsAtPos = gridPSTH_thisSess(thisPos).alignedResponses(iROI, :, :); % [1 x Time x Trials]
                    rowIdx = find(uEl_thisSess_plot == gridStim_thisSess(thisPos, 2), 1);
                    colIdx = find(uAz_thisSess == gridStim_thisSess(thisPos, 1), 1);
                    rawTrialMatrix{rowIdx, colIdx} = squeeze(allTrialsAtPos)'; % [Trials x Time], RAW, no baseline subtraction
                end
                thisRFMapping(iROI).rawAlignedTrials = rawTrialMatrix;
                thisRFMapping(iROI).rawAlignedTrialsBlank = squeeze(blankPSTH_thisSess.alignedResponses(iROI, :, :))'; % [Trials x Time], RAW
            end
            fprintf('  Attached raw (pre-baseline-subtraction) aligned trials for %d ROIs from Response file.\n', nROI_thisSess);
        end

        %% load cross-validated explained variance (from the spatial tuning curve analyses)
        thisMedianCrossValR2 = nan(numel(thisRFMapping), 1);
        thisMeanCrossValR2 = nan(numel(thisRFMapping), 1);
        thisPVals      = nan(numel(thisRFMapping), 1);
        try
            vrSessionData = load(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
                crossValVarName);
       catch ME
            warning('  Could not load RFMapping for %s: %s', thisSessionName, ME.message);
            continue;
        end
        if isfield(vrSessionData, crossValVarName) && ~isempty(vrSessionData.(crossValVarName)) ...
                && isfield(vrSessionData.(crossValVarName), signalName)
            sigData = vrSessionData.(crossValVarName).(signalName);
            if isfield(sigData, 'medianExpVar')
                candidateMedianExpVar = sigData.medianExpVar;
                if length(candidateMedianExpVar) == numel(thisRFMapping)
                    thisMedianCrossValR2 = candidateMedianExpVar;
                else
                    warning('  %s.%s.cvExpVar length (%d) does not match number of boutons (%d) for %s -- filling with NaN.', ...
                        crossValVarName, signalName, numel(candidateMedianExpVar), numel(thisRFMapping), thisSessionName);
                end
            end
            if isfield(sigData, 'pValues')
                candidatePVal = sigData.pValues(:);
                if numel(candidatePVal) == numel(thisRFMapping)
                    thisPVals = candidatePVal;
                else
                    warning('  %s.%s.pValues length (%d) does not match number of boutons (%d) for %s -- filling with NaN.', ...
                        crossValVarName, signalName, numel(candidatePVal), numel(thisRFMapping), thisSessionName);
                end
            end
            if isfield(sigData, 'meanExpVar')
                candidateMeanExpVar = sigData.meanExpVar(:);
                if numel(candidateMeanExpVar) == numel(thisRFMapping)
                    thisMeanCrossValR2 = candidateMeanExpVar;
                else
         
                    warning('  %s.%s.meanExpVar length (%d) does not match number of boutons (%d) for %s -- filling with NaN.', ...
                        crossValVarName, signalName, numel(candidateMeanExpVar), numel(thisRFMapping), thisSessionName);
                end
            end
        else
            fprintf('  No %s.%s for %s (other stimulus/signal not present this session) -- filling with NaN.\n', ...
                crossValVarName, signalName, thisSessionName);
        end

        %% load smi data from response
%          RFMapIdx = find(contains(stimNames, ''));
%         if isempty(RFMapIdx)
%             warning('No RFMapping files for %s — skipping.', thisSessionName);
%             continue;
%         end
       
        %% grid/time-vector consistency check against the reference metadata
        if isempty(RFMappingMetadata)
            % first session with data becomes the shared reference grid
            RFMappingMetadata = sessionRFData.RFMappingMetadata;
        else
            refMeta  = RFMappingMetadata;
            thisMeta = sessionRFData.RFMappingMetadata;
            gridMatches = isequal(refMeta.uAz, thisMeta.uAz) && ...
                          isequal(refMeta.uEl, thisMeta.uEl) && ...
                          isequal(refMeta.timeVector, thisMeta.timeVector);
            if ~gridMatches
                warning(['  Stimulus grid/time vector for %s does not match reference — ' ...
                         'skipping this session to keep pooling consistent.'], thisSessionName);
                continue;
            end
        end
        %% concatenate this session's boutons into the pooled struct array
        if isempty(allRFMapping)
            allRFMapping = thisRFMapping;
        else
            allRFMapping = [allRFMapping; thisRFMapping];
        end
        
        allMedianCrossValR2 = [allMedianCrossValR2; thisMedianCrossValR2]; 
        allMeanCrossValR2   = [allMeanCrossValR2; thisMeanCrossValR2]; 
        allShuffPVal        = [allShuffPVal; thisPVals]; 
        
        sessionLabels = [sessionLabels; repmat({sprintf('%s_%s', thisMouse, thisSessionName)}, ...
            numel(thisRFMapping), 1)]; 
        fprintf('  Added %d boutons (running total: %d).\n', ...
            numel(thisRFMapping), numel(allRFMapping));
    end
end

fprintf('Pooling complete: %d total boutons across %d sessions.\n', ...
    numel(allRFMapping), numel(unique(sessionLabels)));

%% basic responsiveness bookkeeping
numBoutons = numel(allRFMapping);
uAz        = RFMappingMetadata.uAz;
uEl_plot   = RFMappingMetadata.uEl;
timeVector = RFMappingMetadata.timeVector;
respWin    = [0.1 3];   


nAz = length(uAz);
nEl = length(uEl_plot);


%% Recompute responsiveness and RF centers manually for all pooled boutons
% Re-extract the matching indices from the reference metadata
respIdx = timeVector >= respWin(1) & timeVector <= respWin(2);

fprintf('Recomputing responsiveness for %d boutons using spatial grids...\n', numBoutons);

for iROI = 1:numBoutons
    trialMatrix      = allRFMapping(iROI).rawAlignedTrials;       % Cell array [nEl x nAz] 3 x 5 cell array
    bTrials = allRFMapping(iROI).rawAlignedTrialsBlank;           % Matrix [Trials x Time] 
    meanGridResponse = nan(size(allRFMapping(iROI).meanGridResponse));   % will be recomputed below, not reused

    % cheks
    if isempty(trialMatrix) || isempty(bTrials) || isempty(meanGridResponse)
        allRFMapping(iROI).pValANOVA    = NaN;
        allRFMapping(iROI).isResponsive = false;
        allRFMapping(iROI).centerAz     = NaN;
        allRFMapping(iROI).centerEl     = NaN;
        continue;
    end

    % Force types to double numeric matrices
    bTrials = double(bTrials);
    meanGridResponse = double(meanGridResponse);

    % In analyseRFMapping, blankTrialMeans = squeeze(mean(bTrialsCorrected(1, respIdx, :), 2))
    % Depending on how it squeezed, check format: [Trials x Time]
    blankTrialMeans = mean(bTrials(:, respIdx), 2, 'omitnan');

    allTrialMeans = blankTrialMeans(:);
    groupLabels   = repmat(17, numel(blankTrialMeans), 1); % Group 17 is Blank % is this just a label that was assigned


    % loop through grid locations to run ANOVA
    [nEl_curr, nAz_curr] = size(trialMatrix);
    posCounter = 1;

    for r = 1:nEl_curr
        for c = 1:nAz_curr
            correctedTrials = trialMatrix{r, c};

            if ~isempty(correctedTrials)
                correctedTrials = double(correctedTrials);
                posTrialMeans = mean(correctedTrials(:, respIdx), 2, 'omitnan');

                allTrialMeans = [allTrialMeans; posTrialMeans(:)];
                groupLabels   = [groupLabels; repmat(posCounter, numel(posTrialMeans), 1)];

                meanGridResponse(r, c) = mean(posTrialMeans, 'omitnan');   % <-- ADD THIS LINE
            end
            posCounter = posCounter + 1;
        end
    end

    % Clean out NaN trials to preserve ANOVA
    validIdx = ~isnan(allTrialMeans) & ~isnan(groupLabels);
    allTrialMeans = allTrialMeans(validIdx);
    groupLabels   = groupLabels(validIdx);

    if isempty(allTrialMeans) || length(unique(groupLabels)) < 2
        allRFMapping(iROI).pValANOVA    = NaN;
        allRFMapping(iROI).isResponsive = false;
        continue;
    end

    % run inclusion critera
    pValANOVA = anova1(allTrialMeans, groupLabels, 'off');
%     pValANOVA = kruskalwallis(allTrialMeans, groupLabels, 'off');

    blankMean = mean(blankTrialMeans, 'omitnan');
    blankStd  = std(blankTrialMeans, 'omitnan');
    
    % mean responses (avg across trials) at the preferred position
    [prefVal, mI] = max(meanGridResponse(:));
    
    % critera 
    isResponsive = (pValANOVA < 0.05)   && (prefVal > (blankMean + 1 * blankStd));

    % update
    allRFMapping(iROI).pValANOVA    = pValANOVA;
    allRFMapping(iROI).isResponsive = isResponsive;
    
    % recompute spatial coordaintes 
    if isResponsive
        [rPeak, cPeak] = ind2sub(size(meanGridResponse), mI);
        allRFMapping(iROI).centerAz = uAz(cPeak);
        allRFMapping(iROI).centerEl = uEl_plot(rPeak);
    else
        allRFMapping(iROI).centerAz = NaN;
        allRFMapping(iROI).centerEl = NaN;
    end
end


%% Basic responsiveness bookkeeping 
isResp        = [allRFMapping.isResponsive];
respIdxList   = find(isResp);
unrespIdxList = find(~isResp);

fprintf('\n--- After Recomputation ---\n');
fprintf('Total boutons: %d | Responsive (isResponsive): %d (%.1f%%)\n', ...
    numBoutons, numel(respIdxList), 100 * numel(respIdxList) / numBoutons);