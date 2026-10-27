RSPSessions = filterMasterTable('MouseID', {'M25132', 'M25133', 'M26003'}, 'Suite2PPreprocessing', 1, ...
    'DayOfExperience', [5 200], 'HasStimulus', {'LandManipCorridor', 'BaselineCorridor'});
uniqueMice = unique(RSPSessions.MouseID);
numMice = length(uniqueMice);
% Find the maximum number of sessions any one mouse has to set the grid width
sessionCounts = groupsummary(RSPSessions, 'MouseID');
maxSessions = max(sessionCounts.GroupCount);
fileName = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\FOVs\FOV_ComparisonRSP_BaselineAndManip.png';
similarityFileName = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\FOVs\FOV_SimilarityRSP_BaselineAndManip.png';

% Collect reference images + labels alongside the grid plot so we can run
% the pairwise similarity comparison (compareFOVSimilarity, defined as a
% local function at the bottom of this file) afterwards.
collectedRefImages     = {};
collectedMouseIDs      = {};
collectedSessionLabels = {};

% Create figure
figure('Name', 'FOV Sequential Sessions', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout(numMice, maxSessions, 'TileSpacing', 'compact', 'Padding', 'tight');
for m = 1:numMice
    currMouse = uniqueMice(m);
    % Get all rows for this specific mouse, sorted by Session
    mouseSessions = RSPSessions(RSPSessions.MouseID == currMouse, :);
    mouseSessions = sortrows(mouseSessions, 'Session');
    numMouseSessions = height(mouseSessions);
    for s = 1:maxSessions
        nexttile;
        if s <= numMouseSessions
            % Extract data
            mouseID = char(mouseSessions.MouseID(s));
            sessionStr = char(mouseSessions.Session(s)); % e.g., '20260219'
            dayVal  = mouseSessions.DayOfExperience(s);
            % Load and find data (same as before)
            sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID, sessionStr);
            load(sessionFileInfoFilePath);
            try
                stimNames = {sessionFileInfo.stimFiles.name};
                stimIdx = find(contains(stimNames, 'RFMapping', 'IgnoreCase', true), 1);
                dataPath = sessionFileInfo.stimFiles(stimIdx).mergedBonsai2PSuite2pData;
                data = load(dataPath);
                if isfield(data.twoPData.ops, 'refImg')
                    imagesc(data.twoPData.ops.refImg);
                    colormap gray;
                    axis image off;
                    % UPDATED TITLE: Includes Session ID
                    title({sprintf('%s | %s', mouseID, sessionStr); ...
                           sprintf('S%d (Day %d)', s, dayVal)}, 'FontSize', 8);
                    
                    % Stash this reference image for the similarity comparison below
                    collectedRefImages{end+1}     = data.twoPData.ops.refImg; %#ok<SAGROW>
                    collectedMouseIDs{end+1}      = mouseID;                  %#ok<SAGROW>
                    collectedSessionLabels{end+1} = sprintf('%s_%s', mouseID, sessionStr); %#ok<SAGROW>
                else
                    text(0.5, 0.5, 'No refImg', 'HAlign', 'center');
                    axis off;
                end
            catch ME
                axis off;
                title(sprintf('Error %s', sessionStr));
                fprintf('Could not load %s for %s: %s\n', sessionStr, mouseID, ME.message);
            end
        else
            axis off;
        end
    end
end
% Save high-res
exportgraphics(gcf, fileName, 'Resolution', 600);

%% Pairwise FOV similarity across all collected sessions
% Confirms sessions sampled distinct FOVs (not accidental re-imaging of the
% same location) using translation-tolerant normalized cross-correlation.
if length(collectedRefImages) >= 2
    [simMatrix, figSim, fovSimilaritySummary] = compareFOVSimilarity( ...
        collectedRefImages, collectedMouseIDs, collectedSessionLabels);
    exportgraphics(figSim, similarityFileName, 'Resolution', 600);
    
    %% Flag same-mouse session pairs above a similarity threshold for manual review
    % Set this based on what you see in the similarity matrix / any known-duplicate
    % self-match check -- there's no universal "correct" value, it depends on your
    % data's noise floor. Cross-mouse pairs are not checked since they can't be
    % duplicates of each other.
    duplicateThreshold = 0.7;
    duplicateCandidates = flagPotentialDuplicateFOVs(simMatrix, collectedMouseIDs, ...
        collectedSessionLabels, duplicateThreshold, collectedRefImages);
else
    warning('Fewer than 2 reference images were collected; skipping FOV similarity comparison.');
end


%% ===================== LOCAL FUNCTIONS =====================
% Local functions must come after all script code in the file. MATLAB
% supports this (script with local functions) since R2016b -- keeping
% compareFOVSimilarity here means it has no separate path dependency.

function [simMatrix, fig, summary] = compareFOVSimilarity(refImages, mouseIDs, sessionLabels, templateFraction)
    % COMPAREFOVSIMILARITY Computes pairwise similarity between imaging FOVs
    % to demonstrate that sessions sampled distinct fields of view (i.e. are
    % not accidental re-imaging of the same location).
    %
    % INPUTS:
    %   refImages       - cell array of 2D reference images, one per session
    %                      (e.g. suite2p's ops.refImg / meanImg / meanImgE)
    %   mouseIDs        - cell array of strings, one per session, e.g.
    %                      {'M1','M1','M2','M2',...}. Used to separate
    %                      same-mouse vs different-mouse pairs in the summary.
    %   sessionLabels   - cell array of strings for axis labels, e.g.
    %                      {'M1_s1','M1_s2','M2_s1',...}. Defaults to
    %                      'Session1','Session2',... if not provided.
    %   templateFraction- fraction of the image (centered) used as the
    %                      matching template, default 0.7. A value <1 lets
    %                      normxcorr2 search over small translational
    %                      offsets rather than only comparing at zero shift,
    %                      which matters if a revisited FOV was remounted
    %                      with a slight X/Y offset.
    %
    % OUTPUTS:
    %   simMatrix - n x n matrix of peak normalized cross-correlation values
    %               (1 on the diagonal, symmetric off-diagonal)
    %   fig       - handle to the similarity matrix heatmap
    %   summary   - struct with mean/SD of same-mouse vs different-mouse
    %               off-diagonal similarity values (if mouseIDs provided)
    %
    % Requires the Image Processing Toolbox (normxcorr2).
    
    n = length(refImages);
    if nargin < 3 || isempty(sessionLabels)
        sessionLabels = arrayfun(@(x) sprintf('Session%d', x), 1:n, 'UniformOutput', false);
    end
    if nargin < 4 || isempty(templateFraction), templateFraction = 0.7; end
    
    simMatrix = eye(n);
    
    for i = 1:n
        for j = (i+1):n
            imgA = double(refImages{i});
            imgB = double(refImages{j});
            
            % Crop a central template out of imgB to search within imgA,
            % allowing for small translational misalignment between sessions.
            [hB, wB] = size(imgB);
            th = round(hB * templateFraction);
            tw = round(wB * templateFraction);
            rowStart = floor((hB - th) / 2) + 1;
            colStart = floor((wB - tw) / 2) + 1;
            template = imgB(rowStart:rowStart+th-1, colStart:colStart+tw-1);
            
            c = normxcorr2(template, imgA);
            peakVal = max(c(:));
            
            simMatrix(i, j) = peakVal;
            simMatrix(j, i) = peakVal;
        end
    end
    
    %% --- Plot ---
    fig = figure('Color', 'w', 'Position', [100 100 650 600]);
    imagesc(simMatrix, [0 1]);
    axis image;
    colormap(flipud(gray));
    cb = colorbar;
    cb.Label.String = 'Peak normalized cross-correlation';
    cb.TickDirection = 'out';
    cb.Box = 'off';
    
    set(gca, 'XTick', 1:n, 'XTickLabel', sessionLabels, 'XTickLabelRotation', 90, ...
             'YTick', 1:n, 'YTickLabel', sessionLabels, 'TickDir', 'out', 'FontSize', 8);
    title('Pairwise FOV similarity across sessions', 'FontWeight', 'normal', 'FontName', 'Arial');
    box off;
    
    %% --- Summary: same-mouse vs different-mouse pairs ---
    summary = struct();
    if nargin >= 2 && ~isempty(mouseIDs)
        offDiagMask = ~eye(n);
        sameMouseMask = false(n);
        for i = 1:n
            for j = 1:n
                if i ~= j
                    sameMouseMask(i, j) = strcmp(mouseIDs{i}, mouseIDs{j});
                end
            end
        end
        
        sameMouseVals = simMatrix(sameMouseMask & offDiagMask);
        diffMouseVals = simMatrix(~sameMouseMask & offDiagMask);
        
        summary.sameMouse_mean = mean(sameMouseVals);
        summary.sameMouse_sd   = std(sameMouseVals);
        summary.sameMouse_n    = length(sameMouseVals) / 2; % pairs counted twice (symmetric)
        
        summary.diffMouse_mean = mean(diffMouseVals);
        summary.diffMouse_sd   = std(diffMouseVals);
        summary.diffMouse_n    = length(diffMouseVals) / 2;
        
        fprintf('\n=== FOV Similarity Summary ===\n');
        fprintf('Same-mouse, different-session pairs (n=%d):      %.3f +/- %.3f (mean +/- SD)\n', ...
            summary.sameMouse_n, summary.sameMouse_mean, summary.sameMouse_sd);
        fprintf('Different-mouse pairs (n=%d):                    %.3f +/- %.3f (mean +/- SD)\n', ...
            summary.diffMouse_n, summary.diffMouse_mean, summary.diffMouse_sd);
    else
        offDiagVals = simMatrix(~eye(n));
        summary.allPairs_mean = mean(offDiagVals);
        summary.allPairs_sd   = std(offDiagVals);
        fprintf('\n=== FOV Similarity Summary ===\n');
        fprintf('All session pairs: %.3f +/- %.3f (mean +/- SD)\n', summary.allPairs_mean, summary.allPairs_sd);
    end
end

function duplicateCandidates = flagPotentialDuplicateFOVs(simMatrix, mouseIDs, sessionLabels, threshold, refImages)
    % FLAGPOTENTIALDUPLICATEFOVS Flags same-mouse session pairs whose FOV
    % similarity exceeds a threshold, for manual visual review before
    % deciding whether to exclude one of the sessions from analysis.
    %
    % Only same-mouse pairs are checked (cross-mouse pairs can't be
    % duplicates of each other by definition). Flagging is NOT automatic
    % exclusion -- this returns candidates and shows them side-by-side so
    % you can confirm by eye before excluding anything.
    %
    % INPUTS:
    %   simMatrix     - n x n similarity matrix from compareFOVSimilarity
    %   mouseIDs      - cell array of mouse IDs, one per session (same order
    %                   as simMatrix rows/columns)
    %   sessionLabels - cell array of session labels, same order
    %   threshold     - similarity value at/above which a same-mouse pair
    %                   is flagged as a potential duplicate FOV
    %   refImages     - (optional) cell array of the same reference images
    %                   used to build simMatrix. If provided, each flagged
    %                   pair is plotted side-by-side for visual comparison.
    %
    % OUTPUT:
    %   duplicateCandidates - table with MouseID, Session1, Session2,
    %                         Correlation, sorted descending by Correlation
    
    n = length(mouseIDs);
    mouseIDCol  = strings(0, 1);
    session1Col = strings(0, 1);
    session2Col = strings(0, 1);
    corrCol     = [];
    idx1Col     = [];
    idx2Col     = [];
    
    for i = 1:n
        for j = (i+1):n
            if strcmp(mouseIDs{i}, mouseIDs{j}) && simMatrix(i, j) >= threshold
                mouseIDCol(end+1, 1)  = string(mouseIDs{i});    %#ok<AGROW>
                session1Col(end+1, 1) = string(sessionLabels{i}); %#ok<AGROW>
                session2Col(end+1, 1) = string(sessionLabels{j}); %#ok<AGROW>
                corrCol(end+1, 1)     = simMatrix(i, j);        %#ok<AGROW>
                idx1Col(end+1, 1)     = i;                      %#ok<AGROW>
                idx2Col(end+1, 1)     = j;                      %#ok<AGROW>
            end
        end
    end
    
    duplicateCandidates = table(mouseIDCol, session1Col, session2Col, corrCol, ...
        'VariableNames', {'MouseID', 'Session1', 'Session2', 'Correlation'});
    if height(duplicateCandidates) > 0
        [duplicateCandidates, sortIdx] = sortrows(duplicateCandidates, 'Correlation', 'descend');
        idx1Col = idx1Col(sortIdx);
        idx2Col = idx2Col(sortIdx);
    end
    
    fprintf('\n=== Potential duplicate FOVs (same mouse, correlation >= %.2f) ===\n', threshold);
    if height(duplicateCandidates) == 0
        fprintf('None found.\n');
    else
        disp(duplicateCandidates);
    end
    
    % Sessions that never appear in a flagged pair -- i.e. confirmed unique FOVs
    flaggedIdx = unique([idx1Col; idx2Col]);
    uniqueIdx  = setdiff(1:n, flaggedIdx);
    
    fprintf('\n=== Unique (non-flagged) sessions ===\n');
    if isempty(uniqueIdx)
        fprintf('None -- every session was flagged in at least one pair above threshold.\n');
    else
        uniqueMiceHere = unique(mouseIDs(uniqueIdx));
        for mIdx = 1:length(uniqueMiceHere)
            thisMouse = uniqueMiceHere{mIdx};
            theseSessions = sessionLabels(uniqueIdx(strcmp(mouseIDs(uniqueIdx), thisMouse)));
            fprintf('%s: %s\n', thisMouse, strjoin(theseSessions, ', '));
        end
    end
    
    % Side-by-side visual review of each flagged pair, if reference images were given
    if nargin >= 5 && ~isempty(refImages) && height(duplicateCandidates) > 0
        for k = 1:height(duplicateCandidates)
            figure('Name', sprintf('Potential duplicate FOV: %s', duplicateCandidates.MouseID(k)), 'Color', 'w');
            subplot(1, 2, 1);
            imagesc(refImages{idx1Col(k)}); colormap gray; axis image off;
            title(duplicateCandidates.Session1(k), 'Interpreter', 'none', 'FontSize', 9);
            subplot(1, 2, 2);
            imagesc(refImages{idx2Col(k)}); colormap gray; axis image off;
            title(duplicateCandidates.Session2(k), 'Interpreter', 'none', 'FontSize', 9);
            sgtitle(sprintf('%s: r = %.3f', duplicateCandidates.MouseID(k), duplicateCandidates.Correlation(k)), ...
                'FontName', 'Arial', 'FontWeight', 'normal');
        end
    end
end