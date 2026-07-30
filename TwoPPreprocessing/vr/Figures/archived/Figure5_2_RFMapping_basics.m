% This script loads the rf mapping stimulus for all recordings with unique
% fovs - across 3 mice; these are from day 1 of experience
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

%% sessions to load (include all learning and post learning sessions)

pairs=struct;
pairs.M25132 = {'20260219','20260223','20260226','20260228','20260303','20260313','20260306'};
pairs.M25133 = {'20260219','20260223','20260221'};
pairs.M26003 = {'20260316','20260322','20260324','20260325'};


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
        thisRFMapping = thisRFMapping(:);  % force column for consistent concatenation
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
        
        % FIX: Replaced non-existent variable 'allCrossR2' with the correct target arrays,
        % and fixed the concatenation of 'thisPVals' into 'allShuffPVal'.
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


% change response window here 0.1 to 3s 
uAz        = RFMappingMetadata.uAz;
uEl_plot   = RFMappingMetadata.uEl;
timeVector = RFMappingMetadata.timeVector;
% respWin    = RFMappingMetadata.respWin;
respWin = [0.1 3]

nAz = length(uAz);
nEl = length(uEl_plot);

% isResp        = [allRFMapping.isResponsive];
% respIdxList   = find(isResp);
% unrespIdxList = find(~isResp);
% 
% fprintf('Total boutons: %d | Responsive (isResponsive): %d (%.1f%%)\n', ...
%     numel(allRFMapping), numel(respIdxList), 100*numel(respIdxList)/numel(allRFMapping));

%% Recompute responsiveness and RF centers manually for all pooled boutons
numBoutons = numel(allRFMapping);

% Re-extract the matching indices from the reference metadata
respIdx = timeVector >= respWin(1) & timeVector <= respWin(2);

fprintf('Recomputing responsiveness for %d boutons using spatial grids...\n', numBoutons);

for iROI = 1:numBoutons
    trialMatrix      = allRFMapping(iROI).baselineSubtracted;       % Cell array [nEl x nAz]
    bTrialsCorrected = allRFMapping(iROI).baselineSubtractedBlank;  % Matrix [Trials x Time]
    meanGridResponse = allRFMapping(iROI).meanGridResponse;         % Matrix [nEl x nAz]
    
    % cheks 
    if isempty(trialMatrix) || isempty(bTrialsCorrected) || isempty(meanGridResponse)
        allRFMapping(iROI).pValANOVA    = NaN;
        allRFMapping(iROI).isResponsive = false;
        allRFMapping(iROI).centerAz     = NaN;
        allRFMapping(iROI).centerEl     = NaN;
        continue;
    end
    
    % Force types to double numeric matrices
    bTrialsCorrected = double(bTrialsCorrected);
    meanGridResponse = double(meanGridResponse);
    
    % In analyseRFMapping, blankTrialMeans = squeeze(mean(bTrialsCorrected(1, respIdx, :), 2))
    % Depending on how it squeezed, check format: [Trials x Time]
    blankTrialMeans = mean(bTrialsCorrected(:, respIdx), 2, 'omitnan');
    
    allTrialMeans = blankTrialMeans(:);
    groupLabels   = repmat(17, numel(blankTrialMeans), 1); % Group 17 is Blank
    
    % loop through grid locations to run ANOVA
    [nEl_curr, nAz_curr] = size(trialMatrix);
    posCounter = 1;
    
    for r = 1:nEl_curr
        for c = 1:nAz_curr
            % Extract the [Trials x Time] matrix for this specific grid location
            correctedTrials = trialMatrix{r, c};
            
            if ~isempty(correctedTrials)
                correctedTrials = double(correctedTrials);
                
                % Compute the mean response inside the window for each trial
                posTrialMeans = mean(correctedTrials(:, respIdx), 2, 'omitnan');
                
                % Append to vectors
                allTrialMeans = [allTrialMeans; posTrialMeans(:)];
                groupLabels   = [groupLabels; repmat(posCounter, numel(posTrialMeans), 1)];
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
    
    blankMean = mean(blankTrialMeans, 'omitnan');
    blankStd  = std(blankTrialMeans, 'omitnan');
    
    % peak preferred value across the grid
    [prefVal, mI] = max(meanGridResponse(:));
    
    % critera 
    isResponsive = (pValANOVA < 0.05)   && (prefVal > (blankMean + 1 * blankStd));
    
    % update update
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
%%  how much is the amplitude gate costing us vs. ANOVA alone?
anovaOnlyMask   = [allRFMapping.pValANOVA] < 0.05;
lostToAmplitude = anovaOnlyMask & ~isResp;

fprintf('Pass ANOVA alone (p<0.05): %d (%.1f%%)\n', ...
    sum(anovaOnlyMask), 100*sum(anovaOnlyMask)/numel(allRFMapping));
fprintf('Pass ANOVA but excluded by amplitude gate: %d (%.1f%%)\n', ...
    sum(lostToAmplitude), 100*sum(lostToAmplitude)/numel(allRFMapping));

%% are responsive boutons concentrated in a few sessions? yes.... m25132:/
respSessionLabels = sessionLabels(respIdxList);
[uniqueSessResp, ~, sessGroupIdx] = unique(respSessionLabels);
sessCounts = accumarray(sessGroupIdx, 1);
[sessCountsSorted, sortOrderSess] = sort(sessCounts, 'descend');
uniqueSessRespSorted = uniqueSessResp(sortOrderSess);

fprintf('\nSession breakdown of %d responsive (isResponsive) boutons:\n', numel(respIdxList));
for s = 1:numel(uniqueSessRespSorted)
    fprintf('  %s: %d (%.1f%% of responsive boutons)\n', uniqueSessRespSorted{s}, sessCountsSorted(s), ...
        100*sessCountsSorted(s)/numel(respIdxList));
end
if ~isempty(sessCountsSorted) && sessCountsSorted(1) / numel(respIdxList) > 0.5
    fprintf(' Woops: more than half of all responsive boutons come from a single session (%s).\n', ...
        uniqueSessRespSorted{1});
end

%% select which criterion feeds Panel A / B / C / D 
switch selectionMethod
    case 'isResponsive'
        candidateRespIdxList = respIdxList;
    case 'gaussianR2'
        candidateRespIdxList = gaussianRespIdxList;
    otherwise
        error('selectionMethod must be ''isResponsive'' or ''gaussianR2'', got ''%s''.', selectionMethod);
end
candidateUnrespIdxList = setdiff(1:numel(allRFMapping), candidateRespIdxList);

fprintf('\nUsing selectionMethod = ''%s'' --> %d boutons feeding Panel A/B/C/D as "responsive"\n', ...
    selectionMethod, numel(candidateRespIdxList));


%%
% careful when running this - it plots all selected 
plotAllRFResponsiveBoutons
% check
respIdxList   = find([allRFMapping.isResponsive]);
unrespIdxList = find(~[allRFMapping.isResponsive]);
%% manually chosen (picked by eye from plotAllResponsiveBoutons)
responsiveBoutons = [2330,3850, 3189];   %  3189, 
unresponsiveBoutons  = [2637, 69];           % 

exampleSet    = [responsiveBoutons, unresponsiveBoutons];
exampleLabels = {'Responsive', 'Responsive', 'Responsive', ...
                  'Unresponsive', 'Unresponsive'};

nEx  = numel(exampleSet);
colW = 1 / nEx;

figA = figure('Color', 'w', 'Position', [50 50 nEx*550 750], 'Name', 'Panel A: hand-picked examples');

for i = 1:nEx
    b = allRFMapping(exampleSet(i));
    xBase = (i-1)*colW;
    axPanel = axes('Position', [xBase+0.03 0.15 colW-0.05 0.70]);
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    title(axPanel, sprintf('%s\n(Bouton %d)', exampleLabels{i}, exampleSet(i)), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
end


for i = 1:nEx
    b = allRFMapping(exampleSet(i));
    xBase = (i-1)*colW;
    
    axPanel = axes('Position', [xBase+0.03 0.15 colW-0.05 0.70]);
    
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    
    title(axPanel, sprintf('%s\n(Bouton %d)', exampleLabels{i}, exampleSet(i)), ...
        'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
end

% save
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\eg_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figA, 'Visible', 'off');
saveFigureFormats(figA, fullfile(outputDir, 'responsive_unresponsive_examplegrid_linesuperimposed'));

%%  PANEL B 
% Bar chart(s): % responsive vs unresponsive boutons (using selectionMethod)

pctResponsive   = 100 * numel(candidateRespIdxList)   / numel(allRFMapping);
pctUnresponsive = 100 * numel(candidateUnrespIdxList) / numel(allRFMapping);

figB = figure('Color', 'w', 'Position', [50 50 400 400], 'Name', 'Panel B: responsive vs unresponsive');
bar([pctResponsive, pctUnresponsive], 'FaceColor', [0.4 0.4 0.4]);
set(gca, 'XTickLabel', {'Responsive', 'Unresponsive'}, 'Box', 'off', 'TickDir', 'out', ...
    'FontName', 'Arial', 'FontSize', 9);
ylabel('% of boutons', 'FontName', 'Arial', 'FontSize', 9);
ylim([0 100]);
title(sprintf('n = %d boutons (%s)', numel(allRFMapping), selectionMethod), ...
    'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
defaultAxesProperties(gca,true);


% save
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\overall_rep_unresp_hist';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figB, 'Visible', 'off');
saveFigureFormats(figB, fullfile(outputDir, 'responsive_unresponsive_overall_percent'));

%% Panel B2: % responsive per animal (bar = mean across sessions, dots = individual sessions)

boutonMouseID = cellfun(@(s) extractBefore(s, '_'), sessionLabels, 'UniformOutput', false);
nTotalResp = numel(candidateRespIdxList);

respMouseID       = boutonMouseID(candidateRespIdxList);
respSessionLabels = sessionLabels(candidateRespIdxList);

uniqueMiceInData = unique(boutonMouseID, 'stable');
nMiceInData = numel(uniqueMiceInData);

pctRespByMouse = nan(nMiceInData, 1);

figComp = figure('Color', 'w', 'Position', [50 50 500 400], 'Name', 'Composition of responsive pool');
hold on;

fprintf('\nComposition of %d responsive boutons:\n', nTotalResp);
for m = 1:nMiceInData
    pctRespByMouse(m) = 100 * sum(strcmp(respMouseID, uniqueMiceInData{m})) / nTotalResp;
    bar(m, pctRespByMouse(m), 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.7, 'EdgeColor', 'none');

    sessionsForMouse = unique(respSessionLabels(strcmp(respMouseID, uniqueMiceInData{m})), 'stable');
    nSessThisMouse = numel(sessionsForMouse);
    pctPerSession = nan(nSessThisMouse, 1);
    for s = 1:nSessThisMouse
        pctPerSession(s) = 100 * sum(strcmp(respSessionLabels, sessionsForMouse{s})) / nTotalResp;
    end

    % fixed, deterministic horizontal spread (no randomness) -- if only
    % one session, it sits exactly on the bar's x position
    if nSessThisMouse > 1
        xOffsets = linspace(-0.25, 0.25, nSessThisMouse);
    else
        xOffsets = 0;
    end
    xFixed = m + xOffsets;

    scatter(xFixed, pctPerSession, 30, [0.8 0.1 0.1], 'filled', 'MarkerFaceAlpha', 0.7);

    fprintf('  %s: %.1f%% of pool total, from %d session(s): %s\n', ...
        uniqueMiceInData{m}, pctRespByMouse(m), nSessThisMouse, ...
        strjoin(arrayfun(@(x) sprintf('%.1f%%',x), pctPerSession, 'UniformOutput', false), ', '));
end

set(gca, 'XTick', 1:nMiceInData, 'XTickLabel', uniqueMiceInData, 'Box', 'off', 'TickDir', 'out', ...
    'FontName', 'Arial', 'FontSize', 9);
ylabel('% of responsive pool', 'FontName', 'Arial', 'FontSize', 9);
ytickformat(gca, '%g%%');
xlim([0.5, nMiceInData+0.5]);
title(sprintf('n = %d responsive boutons -- bars = per mouse, dots = per session', nTotalResp), ...
    'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
defaultAxesProperties(gca, true);

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\responsive_percent_per_animal';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figComp, 'Visible', 'off');
saveFigureFormats(figComp, fullfile(outputDir, 'responsive_percent_per_animal_with_sessions'));

%%  PANEL C 

prefAz = [allRFMapping(candidateRespIdxList).centerAz];
prefEl = [allRFMapping(candidateRespIdxList).centerEl];

figC = figure('Position', [100, 100, 850, 280]);

% az histogram
axC1 = subplot(1,3,1);
if length(uAz) > 1
    azStep = mean(diff(uAz));
    azEdges = (min(uAz) - azStep/2) : azStep : (max(uAz) + azStep/2);
else
    azEdges = uAz;
end
azCounts  = histcounts(prefAz, azEdges);
azPct     = 100 * azCounts / numel(prefAz);
azCenters = azEdges(1:end-1) + diff(azEdges)/2;
bar(azCenters, azPct, 1, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Preferred Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel('% of boutons', 'FontName', 'Arial', 'FontSize', 9);
set(axC1, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', sort(uAz));
drawnow;
defaultAxesProperties(axC1, true);

% el histogram
axC2 = subplot(1,3,2);
if length(uEl_plot) > 1
    elStep = abs(mean(diff(uEl_plot)));
    elEdges = (min(uEl_plot) - elStep/2) : elStep : (max(uEl_plot) + elStep/2);
else
    elEdges = uEl_plot;
end
elCounts  = histcounts(prefEl, elEdges);
elPct     = 100 * elCounts / numel(prefEl);
elCenters = elEdges(1:end-1) + diff(elEdges)/2;
bar(elCenters, elPct, 1, 'FaceColor', [0.8 0.4 0.3]);
xlabel('Preferred Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel('% of boutons', 'FontName', 'Arial', 'FontSize', 9);
set(axC2, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', sort(uEl_plot));
drawnow;
defaultAxesProperties(axC2, true);

% bubble scatter
axC3 = subplot(1,3,3);
[uniqueCoords, ~, idx] = unique([prefAz(:), prefEl(:)], 'rows');
counts = accumarray(idx, 1);
bubbleSizes = counts * 80;
scatter(uniqueCoords(:,1), uniqueCoords(:,2), bubbleSizes, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Preferred Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel('Preferred Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
set(axC3, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
xlim(axC3, [min(uAz)-10, max(uAz)+10]);
ylim(axC3, [min(uEl_plot)-10, max(uEl_plot)+10]);
set(axC3, 'XTick', sort(uAz), 'YTick', sort(uEl_plot));
hold(axC3, 'on');
for i = 1:size(uniqueCoords,1)
    text(uniqueCoords(i,1), uniqueCoords(i,2), sprintf('%d', counts(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
end
axis(axC3, 'square');
drawnow;
defaultAxesProperties(axC3, true);


outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\pref_az_el';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figC, 'Visible', 'off');
saveFigureFormats(figC, fullfile(outputDir, 'responsive_pref_az_el'));


%% Test: is there any within session differences of preferred Az/El positions
uniqueSessResp = unique(respSessionLabels, 'stable');  % from your composition code
nSessResp = numel(uniqueSessResp);
nPossiblePositions = numel(uAz) * numel(uEl_plot);

sessionDiversity = table('Size', [nSessResp, 4], ...
    'VariableTypes', {'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Session', 'nBoutons', 'nUniquePositions', 'pctPositionsCovered'});

for s = 1:nSessResp
    sessMask = strcmp(respSessionLabels, uniqueSessResp{s});
    sessAz = prefAz(sessMask);
    sessEl = prefEl(sessMask);
    uniquePos = unique([sessAz(:), sessEl(:)], 'rows');

    sessionDiversity.Session(s) = uniqueSessResp{s};
    sessionDiversity.nBoutons(s) = sum(sessMask);
    sessionDiversity.nUniquePositions(s) = size(uniquePos, 1);
    sessionDiversity.pctPositionsCovered(s) = 100 * size(uniquePos,1) / nPossiblePositions;
end

disp(sessionDiversity);

% one bubble plot per session (small multiples)
nCols = 4;
nRows = ceil(nSessResp / nCols);
figDiv = figure('Color', 'w', 'Position', [50 50 nCols*250 nRows*250], ...
    'Name', 'Preferred position diversity by session');

for s = 1:nSessResp
    sessMask = strcmp(respSessionLabels, uniqueSessResp{s});
    sessAz = prefAz(sessMask);
    sessEl = prefEl(sessMask);
    [uniquePos, ~, idx] = unique([sessAz(:), sessEl(:)], 'rows');
    counts = accumarray(idx, 1);

    axS = subplot(nRows, nCols, s);
    scatter(uniquePos(:,1), uniquePos(:,2), counts*60, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.6);
    xlim([min(uAz)-10, max(uAz)+10]);
    ylim([min(uEl_plot)-10, max(uEl_plot)+10]);
    set(axS, 'XTick', sort(uAz), 'YTick', sort(uEl_plot), 'FontSize', 6, ...
        'Box', 'off', 'TickDir', 'out');
    axis(axS, 'square');
    title(sprintf('%s\nn=%d, %d/%d pos', uniqueSessResp{s}, sum(sessMask), ...
        size(uniquePos,1), nPossiblePositions), 'FontSize', 7, 'FontWeight', 'normal');
end

% inspectClusteredBoutons
% inspectSessionDiversity
%%  Panel E-ish: cross-val R2 relationship
% is RF-mapping responsiveness related to how well-explained a
% bouton is by the OTHER (e.g. moving-dot) stimulus's cross-validated model?
if any(~isnan(allMedianCrossValR2))
    figE = figure('Color', 'w', 'Position', [50 50 500 400], 'Name', 'Cross-val R2 vs RF responsiveness');
    axC3 = axes(figE);
    hold(axC3, 'on');
    edgesShared = linspace(-0.8, 1, 40);
    histogram(axC3, allMedianCrossValR2(candidateRespIdxList), edgesShared, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', [0.1 0.5 0.1], 'LineWidth', 2, ...
        'DisplayName', 'RF-responsive');
    histogram(axC3, allMedianCrossValR2(candidateUnrespIdxList), edgesShared, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', [0.6 0.2 0.2], 'LineWidth', 2, ...
        'DisplayName', 'RF-non-responsive');
    xline(axC3, 0.1, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
        'Label', '0.1', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);

    % --- combined criterion: EV > 0.1 AND significant vs shuffle ---
    % computed for BOTH groups, so you can see how many RF-non-responsive
    % boutons still pass the other stimulus's criteria (i.e. "responding
    % there but not here")
    sigThresh = 0.01;  % matches the 99th-percentile shuffle threshold used elsewhere

    passing_resp = (allMedianCrossValR2(candidateRespIdxList) > 0.1) & (allShuffPVal(candidateRespIdxList) < sigThresh);
    pct_passing_resp = 100 * sum(passing_resp) / numel(passing_resp);

    passing_unresp = (allMedianCrossValR2(candidateUnrespIdxList) > 0.1) & (allShuffPVal(candidateUnrespIdxList) < sigThresh);
    pct_passing_unresp = 100 * sum(passing_unresp) / numel(passing_unresp);

    fprintf('\nBoutons meeting BOTH criteria (EV>0.1 AND p<%.2f vs shuffle):\n', sigThresh);
    fprintf('  RF-responsive:     %.1f%% (%d/%d)\n', pct_passing_resp, sum(passing_resp), numel(passing_resp));
    fprintf('  RF-non-responsive: %.1f%% (%d/%d)\n', pct_passing_unresp, sum(passing_unresp), numel(passing_unresp));

    % --- annotations on the figure, color-matched to each group ---
    yl = ylim(axC3);
    text(axC3, 0.12, yl(2)*0.90, sprintf('Passing both: %.1f%%', pct_passing_resp), ...
        'Color', [0.1 0.5 0.1], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(axC3, 0.12, yl(2)*0.82, sprintf('Passing both: %.1f%%', pct_passing_unresp), ...
        'Color', [0.6 0.2 0.2], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    xlabel(axC3, 'Cross-validated R^2 (other stimulus)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(axC3, 'Proportion of boutons', 'FontName', 'Arial', 'FontSize', 9);
    legend(axC3, 'Location', 'best');
    set(axC3, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
    
    hasCrossVal = ~isnan(allMedianCrossValR2);
    nSessionsWithCrossVal = numel(unique(sessionLabels(hasCrossVal)));
    nSessionsTotal = numel(unique(sessionLabels));

    fprintf('\nCross-val data present for %d/%d sessions.\n', nSessionsWithCrossVal, nSessionsTotal);

    drawnow;
    defaultAxesProperties(axC3, true);
    title(axC3, sprintf('Do RF-responsive boutons also have better cross-val tuning elsewhere? (%d/%d sessions)', ...
        nSessionsWithCrossVal, nSessionsTotal), ...
        'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'normal');
else
    fprintf('\nallCrossValR2 is all-NaN -- check crossValVarName at the top of this script.\n');
end


outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\spatialTuning_EVComparions';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
set(figE, 'Visible', 'off');
saveFigureFormats(figE, fullfile(outputDir, 'spatialEV_vs_resp_unresp'));

%% gaussian fit (allen sdk)
% % Fit on all boutons that pass ANOVA alone (923 in your last run), independent
% % of the amplitude gate -- this lets 'gaussianR2' serve as a genuinely
% % alternative/independent responsiveness criterion to compare against isResponsive.
% 
% gaussianFitCandidates = find(anovaOnlyMask);
% nCand = numel(gaussianFitCandidates);
% 
% fitR2_all            = nan(nCand, 1);
% fitFWHM_az_all       = nan(nCand, 1);
% fitFWHM_el_all       = nan(nCand, 1);
% fitFWHM_combined_all = nan(nCand, 1);
% 
% for k = 1:nCand
%     b = allRFMapping(gaussianFitCandidates(k));
%     Z = b.meanGridResponse;
% 
%     try
%         [fitFun, fwhmX, fwhmY] = fitGaussian2D(uAz, uEl_plot, Z);
% 
%         [Xg, Yg] = meshgrid(uAz, uEl_plot);
%         Zhat = fitFun(Xg, Yg);
%         validMask = ~isnan(Z);
% 
%         SSres = sum((Z(validMask) - Zhat(validMask)).^2);
%         SStot = sum((Z(validMask) - mean(Z(validMask))).^2);
%         r2 = 1 - SSres/SStot;
% 
%         fitR2_all(k)            = r2;
%         fitFWHM_az_all(k)       = fwhmX;
%         fitFWHM_el_all(k)       = fwhmY;
%         fitFWHM_combined_all(k) = sqrt(fwhmX * fwhmY);
%     catch ME
%         warning('Gaussian fit failed for bouton %d: %s', gaussianFitCandidates(k), ME.message);
%     end
% end
% 
% figure('Name', 'Gaussian fit R2 distribution (ANOVA-passing boutons)');
% histogram(fitR2_all, 30);
% xlabel('R^2 of 2D Gaussian fit'); ylabel('# boutons');
% title(sprintf('n = %d ANOVA-passing boutons -- check for a bimodal split before setting r2Thresh', nCand));
% 
% gaussianRespMaskLocal = fitR2_all > r2Thresh;
% gaussianRespIdxList   = gaussianFitCandidates(gaussianRespMaskLocal);
% 
% fprintf('\nGaussian-fit criterion (R2 > %.2f, fit on %d ANOVA-passing boutons): %d well-fit (%.1f%% of candidates)\n', ...
%     r2Thresh, nCand, numel(gaussianRespIdxList), 100*numel(gaussianRespIdxList)/nCand);
% 
% overlapCount = numel(intersect(respIdxList, gaussianRespIdxList));
% fprintf('Overlap between isResponsive (%d) and gaussianR2 (%d) sets: %d boutons in both\n', ...
%     numel(respIdxList), numel(gaussianRespIdxList), overlapCount);

%% 
%% Panel A: per-bouton Gaussian fit  (allen SDK: fitgaussian2D) https://github.com/AllenInstitute/AllenSDK
% % Fit Gaussian specifically to the chosen candidate pool for width-based
% % example selection (re-fit here since gaussianFitCandidates above was the
% % broader ANOVA-passing pool, not necessarily == candidateRespIdxList).
% 
% % Set up arrays for the optimization loop
% nResp = numel(candidateRespIdxList);
% fitR2            = nan(nResp, 1);
% fitFWHM_combined = nan(nResp, 1);
% 
% for k = 1:nResp
%     %  data for this responsive bouton
%     b = allRFMapping(candidateRespIdxList(k));
%     Z = b.meanGridResponse;
%     
%     try
%         % Fit 2D Gaussian to get the spatial tuning curves and width parameters
%         [fitFun, fwhmX, fwhmY] = fitGaussian2D(uAz, uEl_plot, Z);
%         
%         % Generate a grid to evaluate our model fit
%         [Xg, Yg] = meshgrid(uAz, uEl_plot);
%         Zhat = fitFun(Xg, Yg);
%         
%         % Handle any missing data blocks cleanly
%         validMask = ~isnan(Z);
%         
%         % Compute R-squared to evaluate fit quality
%         SSres = sum((Z(validMask) - Zhat(validMask)).^2);
%         SStot = sum((Z(validMask) - mean(Z(validMask))).^2);
%         fitR2(k) = 1 - SSres/SStot;
%         
%         %  geometric mean of X/Y widths as a single size metric
%         fitFWHM_combined(k) = sqrt(fwhmX * fwhmY);
%         
%     catch ME
%         % Don't crash the script if a noisy bouton fails optimization
%         warning('Gaussian fit failed for bouton %d: %s', candidateRespIdxList(k), ME.message);
%     end
% end
% 
% % Filter for boutons that actually match a Gaussian profile well
% wellFitIdxInResp = find(fitR2 > r2Thresh);
% wellFitBoutonIdx = candidateRespIdxList(wellFitIdxInResp);
% widths           = fitFWHM_combined(wellFitIdxInResp);
% 
% % Sort remaining population from narrowest to widest receptive fields
% [sortedWidths, sortOrder] = sort(widths);
% nWellFit = numel(sortedWidths);
% 
% % Select specific percentiles to grab diverse RF sizes for Panel A
% pctTargets = [10 50 90];  % narrow, medium, broad representatives
% exampleRespIdxList = nan(1,3);
% for i = 1:3
%     targetRank = round(pctTargets(i)/100 * nWellFit);
%     targetRank = max(1, min(nWellFit, targetRank)); % boundary guard
%     exampleRespIdxList(i) = wellFitBoutonIdx(sortOrder(targetRank));
% end
% 
% % Print selection summary to command window
% fprintf('\nSelected responsive examples (bouton IDs): narrow=%d, medium=%d, broad=%d\n', ...
%     exampleRespIdxList);
% fprintf('Corresponding combined FWHM (deg): %.1f, %.1f, %.1f\n', ...
%     sortedWidths(max(1,min(nWellFit,round(pctTargets/100*nWellFit)))));
% 
% % Pick two non-responsive boutons at different noise/amplitude levels as controls
% unrespAmps       = [allRFMapping(candidateUnrespIdxList).peakAmplitude];
% [~, medRankIdx]  = sort(unrespAmps);
% exampleUnrespIdx1 = candidateUnrespIdxList(medRankIdx(round(numel(medRankIdx) * 0.25)));
% exampleUnrespIdx2 = candidateUnrespIdxList(medRankIdx(round(numel(medRankIdx) * 0.75)));
% 
% % Combine examples and labels for the final multi-panel plot
% exampleSet    = [exampleRespIdxList, exampleUnrespIdx1, exampleUnrespIdx2];
% exampleLabels = {'Narrow RF', 'Medium RF', 'Broad RF', 'Non-responsive', 'Non-responsive'};
% nEx  = numel(exampleSet);
% colW = 1 / nEx;
% 
% figA = figure('Color', 'w', 'Position', [50 50 nEx*550 750], 'Name', 'Panel A: RF examples');
% 
% for i = 1:nEx
%     b = allRFMapping(exampleSet(i));
%     xBase = (i-1)*colW;
%     
%     % Format axes positions dynamically across the canvas width
%     axPanel = axes('Position', [xBase+0.03 0.12 colW-0.06 0.75]);
%     
%     % Generate the final raw trace map visualization
%     plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
%     
%     title(axPanel, sprintf('%s (bouton %d)', exampleLabels{i}, exampleSet(i)), ...
%         'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');
% end

%%%%  PANEL D 
% Population receptive field widths (re-centered on each bouton's own peak).
% nResp = numel(candidateRespIdxList);
% azProfiles = nan(nResp, nAz);
% elProfiles = nan(nResp, nEl);
% grid2D_stack = nan(nEl, nAz, nResp);
% 
% for k = 1:nResp
%     b = allRFMapping(candidateRespIdxList(k));
%     g = b.meanGridResponse;
% 
%     peakVal = max(g(:), [], 'omitnan');
%     if peakVal <= 0 || isnan(peakVal), continue; end
%     gNorm = g / peakVal;
% 
%     [~, mI] = max(gNorm(:), [], 'omitnan');
%     [rPeak, cPeak] = ind2sub(size(gNorm), mI);
% 
%     azProfiles(k, :) = gNorm(rPeak, :);
%     elProfiles(k, :) = gNorm(:, cPeak)';
% 
%     shiftRows = round(nEl/2) - rPeak;
%     shiftCols = round(nAz/2) - cPeak;
%     grid2D_stack(:, :, k) = circshift(gNorm, [shiftRows, shiftCols]);
% end
% 
% azOffset = uAz - uAz(round(nAz/2));
% elOffset = uEl_plot - uEl_plot(round(nEl/2));
% 
% meanAzProfile = mean(azProfiles, 1, 'omitnan');
% meanElProfile = mean(elProfiles, 1, 'omitnan');
% meanGrid2D    = mean(grid2D_stack, 3, 'omitnan');
% 
% [azFit, azFWHM] = fitGaussian1D(azOffset, meanAzProfile);
% [elFit, elFWHM] = fitGaussian1D(elOffset, meanElProfile);
% [fit2D, fwhmX, fwhmY] = fitGaussian2D(azOffset, elOffset, meanGrid2D);
% 
% fprintf('\nPopulation RF widths (n = %d %s boutons):\n', nResp, selectionMethod);
% fprintf('  Azimuth   FWHM (1D fit): %.1f deg\n', azFWHM);
% fprintf('  Elevation FWHM (1D fit): %.1f deg\n', elFWHM);
% fprintf('  Azimuth   FWHM (2D fit): %.1f deg\n', fwhmX);
% fprintf('  Elevation FWHM (2D fit): %.1f deg\n', fwhmY);
% 
% figD = figure('Color', 'w', 'Position', [50 50 1200 400], 'Name', 'Panel D: population RF widths');
% 
% subplot(1,3,1);
% plot(azOffset, meanAzProfile, 'o', 'Color', [0.3 0.5 0.8], 'MarkerFaceColor', [0.3 0.5 0.8]);
% hold on;
% xFine = linspace(min(azOffset), max(azOffset), 200);
% plot(xFine, azFit(xFine), 'k-', 'LineWidth', 1.5);
% xlabel('Azimuth offset from peak (\circ)', 'FontName', 'Arial', 'FontSize', 9);
% ylabel('Normalized response', 'FontName', 'Arial', 'FontSize', 9);
% title(sprintf('FWHM = %.1f\\circ', azFWHM), 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
% set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
% 
% subplot(1,3,2);
% plot(elOffset, meanElProfile, 'o', 'Color', [0.8 0.4 0.3], 'MarkerFaceColor', [0.8 0.4 0.3]);
% hold on;
% xFine = linspace(min(elOffset), max(elOffset), 200);
% plot(xFine, elFit(xFine), 'k-', 'LineWidth', 1.5);
% xlabel('Elevation offset from peak (\circ)', 'FontName', 'Arial', 'FontSize', 9);
% ylabel('Normalized response', 'FontName', 'Arial', 'FontSize', 9);
% title(sprintf('FWHM = %.1f\\circ', elFWHM), 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
% set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
% 
% subplot(1,3,3);
% imagesc(azOffset, elOffset, meanGrid2D);
% axis image; set(gca, 'YDir', 'normal');
% colormap(gca, 'parula'); colorbar;
% xlabel('Azimuth offset (\circ)', 'FontName', 'Arial', 'FontSize', 9);
% ylabel('Elevation offset (\circ)', 'FontName', 'Arial', 'FontSize', 9);
% title(sprintf('2D fit: FWHM_x=%.1f\\circ, FWHM_y=%.1f\\circ', fwhmX, fwhmY), ...
%     'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal');
% set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
% 


