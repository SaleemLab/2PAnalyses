function [roisToKeep, roisToDiscard, groups, corrMatrixCombined] = findHightlyCorrelatedROIs(sessionFileInfo, signalToUse, VRThreshold_ForBoutonMatch, IncludeVRRuns, NonVRThreshold_ForBoutonMatch, plotFlag)
% FINDHIGHTLYCORRELATEDROIS identifies groups of highly correlated ROIs 
% using thresholds on VR and/or NonVR data. 
% 
% NOTE: This function calculates matrices and APPENDS them to the .mat file 
% specified in sessionFileInfo.otherSessFilePaths.sessionROIData.

%% Set defaults and initialize
if nargin < 2 || isempty(signalToUse); signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 3 || isempty(VRThreshold_ForBoutonMatch); VRThreshold_ForBoutonMatch = 0.4; end
if nargin < 4 || isempty(IncludeVRRuns); IncludeVRRuns = false; end % 
if nargin < 5 || isempty(NonVRThreshold_ForBoutonMatch); NonVRThreshold_ForBoutonMatch = 0.7; end
if nargin < 6 || isempty(plotFlag); plotFlag = true; end

% Initialize trace containers
traces_VR = [];
traces_NonVR = [];

if IncludeVRRuns
    dataScope = 'CombinedVR_NonVR';
else
    dataScope = 'Only_NonVR';
    fprintf('IncludeVRRuns is FALSE. Skipping VR data processing entirely.\n');
end

%% Load VR Data (Conditional)
if IncludeVRRuns
    VRStimIndices = find(contains({sessionFileInfo.stimFiles.name}, 'Corridor') & ...
        ~contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );
    numVRFiles = length(VRStimIndices);
    
    if numVRFiles >= 1
        fprintf('Loading and concatenating %d VR files...\n', numVRFiles);
        allDFFData_VR = cell(1, numVRFiles);
        for thisVRFile = 1:numVRFiles
            thisStimIdx = VRStimIndices(thisVRFile);
            processedTwoPDataPath = sessionFileInfo.stimFiles(thisStimIdx).processedMergedBonsaiSuite2pData;
            if exist(processedTwoPDataPath, "file")
                processedTwoPData = load(processedTwoPDataPath, "zScoredProcessedSignals");
                allDFFData_VR{thisVRFile} = processedTwoPData.zScoredProcessedSignals.(signalToUse);
            else
                warning('Missing processedTwoPData file for VR: %s \n', sessionFileInfo.stimFiles(thisStimIdx).name);
            end
        end
        traces_VR = [allDFFData_VR{:}];
    else
        warning('No VRCorr files found. VR correlation matrix will be empty.');
    end
end

%% Load NonVR Data
isNotVRStim = find(~contains({sessionFileInfo.stimFiles.name}, 'Corridor') & ...
    ~contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );
numNonVRFiles = length(isNotVRStim);
if numNonVRFiles >= 1
    fprintf('Loading and concatenating %d NonVR files...\n', numNonVRFiles);
    allDFFData_NonVR = cell(1, numNonVRFiles);
    for thisNonVRFile = 1:numNonVRFiles
        try
            thisOtherStimIdx = isNotVRStim(thisNonVRFile);
            processedTwoPDataPath = sessionFileInfo.stimFiles(thisOtherStimIdx).processedMergedBonsaiSuite2pData;
            if exist(processedTwoPDataPath, "file")
                processedTwoPData = load(processedTwoPDataPath);
                allDFFData_NonVR{thisNonVRFile} = processedTwoPData.zScoredProcessedSignals.(signalToUse);
            else
                sprintf('Missing processedTwoPData file for NonVR: %s \n', sessionFileInfo.stimFiles(thisOtherStimIdx).name);
            end
        catch
        end
    end
    traces_NonVR = [allDFFData_NonVR{:}];
else
    warning('No NonVR files found. NonVR correlation matrix will be empty.');
end

%% Prepare Combined Traces and Check Consistency
traces_combined = [traces_VR, traces_NonVR];
if isempty(traces_combined)
    error('Traces variable is empty after loading data.');
end
[nROIs, ~] = size(traces_combined);
fprintf('Loaded %d ROIs. Starting correlation analysis...\n', nROIs);

%% Calculate Correlation Matrices and Adjacency Matrices
% VR Correlation and Thresholding 
if IncludeVRRuns && ~isempty(traces_VR)
    VRCorrMatrix = corr(traces_VR');
    adjMatrix_VR = VRCorrMatrix > VRThreshold_ForBoutonMatch;
    fprintf('VR correlation matrix calculated. Applying VR Threshold: %.2f\n', VRThreshold_ForBoutonMatch);
else
    VRCorrMatrix = zeros(nROIs);
    adjMatrix_VR = false(nROIs);
end

% NonVR Correlation and Thresholding
if ~isempty(traces_NonVR)
    NonVRCorrMatrix = corr(traces_NonVR');
    adjMatrix_NonVR = NonVRCorrMatrix > NonVRThreshold_ForBoutonMatch;
    fprintf('NonVR correlation matrix calculated. Applying NonVR Threshold: %.2f\n', NonVRThreshold_ForBoutonMatch);
else
    NonVRCorrMatrix = zeros(nROIs);
    adjMatrix_NonVR = false(nROIs);
end

% Max correlation across active conditions
if IncludeVRRuns
    corrMatrixCombined = max(VRCorrMatrix, NonVRCorrMatrix);
    adjMatrix_final = adjMatrix_VR | adjMatrix_NonVR; 
    fprintf('Combining adjacency matrices using OR logic.\n');
else
    corrMatrixCombined = NonVRCorrMatrix;
    adjMatrix_final = adjMatrix_NonVR;
    fprintf('Using NonVR adjacency matrix only.\n');
end

adjMatrix_final(logical(eye(nROIs))) = false; % Remove self-correlation

%% Find Groups and Pruning
fprintf('Calculating SNR scores (max z-score) for %d ROIs using dynamic trace scope...\n', nROIs);
snrScores = max(traces_combined, [], 2);

% Use graph theory to find connected components
fprintf('Finding robustly correlated groups...\n');
G = graph(adjMatrix_final, 'upper'); 
groups = conncomp(G);
uniqueGroups = unique(groups);
fprintf('Found %d unique groups (including singletons).\n', length(uniqueGroups));

% Keep only highest SNR ROI from each group (Pruning logic)
fprintf('Pruning groups to keep highest SNR ROI...\n');
roisToKeep_mask = false(nROIs, 1);
roisToDiscard_mask = false(nROIs, 1);
for thisGroup = 1:length(uniqueGroups)
    thisGroupID = uniqueGroups(thisGroup);
    roisInGroup_idx = find(groups == thisGroupID);
    
    if isscalar(roisInGroup_idx)
        roisToKeep_mask(roisInGroup_idx) = true;
    else
        groupSnrScores = snrScores(roisInGroup_idx);
        [~, localWinnerIdx] = max(groupSnrScores); 
        globalWinnerIdx = roisInGroup_idx(localWinnerIdx);
        
        roisToKeep_mask(globalWinnerIdx) = true;
        roisInGroup_idx(localWinnerIdx) = []; 
        roisToDiscard_mask(roisInGroup_idx) = true;
    end
end
roisToKeep = find(roisToKeep_mask);
roisToDiscard = find(roisToDiscard_mask);

%% Copy to save within a structure 
outputFilePath = sessionFileInfo.otherSessFilePaths.sessionROIData;
highlyCorrBoutons.includeVRRuns = IncludeVRRuns;
highlyCorrBoutons.vrThreshold = VRThreshold_ForBoutonMatch;
highlyCorrBoutons.nonvrThreshold = NonVRThreshold_ForBoutonMatch;
highlyCorrBoutons.vrAdjMatrix = adjMatrix_VR;
highlyCorrBoutons.nonvrAdjMatrix = adjMatrix_NonVR;
highlyCorrBoutons.dataScope = dataScope; 
highlyCorrBoutons.roisToKeep = roisToKeep;
highlyCorrBoutons.roisToDiscard = roisToDiscard; 
highlyCorrBoutons.adjMatrixFinal = adjMatrix_final; 
highlyCorrBoutons.roisWithAssignedGroups = groups';
highlyCorrBoutons.corrMatrixCombined = corrMatrixCombined;  

fprintf('Analysis complete. Kept %d ROIs, discarded %d ROIs.\n', ...
        length(highlyCorrBoutons.roisToKeep), length(highlyCorrBoutons.roisToDiscard));

discardedGroupIDs = groups(highlyCorrBoutons.roisToDiscard);
highlyCorrBoutons.discardedGroupIDs = discardedGroupIDs';
uniqueDiscardedGroups = unique(discardedGroupIDs);
highlyCorrBoutons.uniqueDiscardedGroups = uniqueDiscardedGroups; 

%% Save to sessionROIData 
if exist(outputFilePath, 'file') == 2
    fprintf('Appending new correlation variables to existing file: %s\n', outputFilePath);
    save(outputFilePath, "highlyCorrBoutons", '-append');
    fprintf('Successfully appended correlation data.\n');
else
    error('sessionROIData file not found at: %s. Cannot append correlation data.', outputFilePath);
end

%% Plotting 
if plotFlag
    fig1 = figure('Name', 'Correlation Matrix and Adjacency Matrices', 'NumberTitle', 'off');
    
    subplot(1, 2, 1);
    imagesc(corrMatrixCombined);
    axis image;
    colorbar;
    title(sprintf('Correlation Matrix (%d ROIs) [%s]', nROIs, dataScope));
    xlabel('ROI Index');
    ylabel('ROI Index');
    
    subplot(1, 2, 2);
    imagesc(adjMatrix_final);
    colormap gray;
    axis image;
    if IncludeVRRuns
        title(sprintf('Final Adjacency Matrix (VR > %.2f OR NonVR > %.2f)', VRThreshold_ForBoutonMatch, NonVRThreshold_ForBoutonMatch));
    else
        title(sprintf('Final Adjacency Matrix (NonVR > %.2f Only)', NonVRThreshold_ForBoutonMatch));
    end
    xlabel('ROI Index');
    ylabel('ROI Index'); 
    
    % Plotting traces for two unique discarded groups
    if numel(uniqueDiscardedGroups) >= 2
        fig2 = figure('Name', 'Correlated ROI Groups', 'NumberTitle', 'off');
        subplot(211)
        gM1 = find(groups == uniqueDiscardedGroups(1)); 
        plot(traces_combined(gM1, :)');
        title(['Discarded Group ID: ' num2str(uniqueDiscardedGroups(1))]);
        ylabel('\Delta F/F');
    
        subplot(212)
        gM2 = find(groups == uniqueDiscardedGroups(2)); 
        plot(traces_combined(gM2, :)');
        title(['Discarded Group ID: ' num2str(uniqueDiscardedGroups(2))]);
        ylabel('\Delta F/F');
    else
        warning('Could not find 2 unique groups in roisToDiscard to plot.');
    end
end
end
    
    

%Plot spatial tuning curves
% load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_CombinedRuns.mat")
% lapActivity = response.lapPositionActivity.dFFNeuropilCorrected; 
% 
% % 
% maxValidROI = size(lapActivity, 1); 
% 
% roisToKeep_filtered = roisToKeep(roisToKeep <= maxValidROI);
% roisToDiscard_filtered = roisToDiscard(roisToDiscard <= maxValidROI);
% 
% % --- Continue plotting using the filtered indices ---
% meanAll = squeeze(mean(lapActivity, 2, 'omitnan'));
% normAll = normalize(meanAll, 2, 'range');
% [~, peakIdx] = max(normAll, [], 2);
% [~, sortIdx] = sort(peakIdx);
% 
% meanROIsToKeep = squeeze(mean(lapActivity(roisToKeep_filtered,:,:), 2, 'omitnan'));
% normROIsToKeep = normalize(meanROIsToKeep, 2, 'range');
% [~, peakIdxToKeep] = max(normROIsToKeep, [], 2);
% [~, sortIdxToKeep] = sort(peakIdxToKeep);
% meanROIsToDiscard = squeeze(mean(lapActivity(roisToDiscard_filtered,:,:), 2, 'omitnan'));
% normROIsToDiscard = normalize(meanROIsToDiscard, 2, 'range');
% [~, peakIdxToDiscard] = max(normROIsToDiscard, [], 2);
% [~, sortIdxToDiscard] = sort(peakIdxToDiscard);
% 
% fig3 = figure;
% subplot(131)
% imagesc(normAll(sortIdx, :));
% xlabel('Position (cm)');
% ylabel('ROIs');
% title('All')
% ylabel(colorbar, 'Activity (normalised)');
% 
% subplot(132)
% imagesc(normROIsToKeep(sortIdxToKeep, :));
% xlabel('Position (cm)');
% ylabel('ROIs');
% title('Single pringle')
% ylabel(colorbar, 'Activity (normalised)');
% 
% subplot(133)
% imagesc(normROIsToDiscard(sortIdxToDiscard, :));
% xlabel('Position (cm)');
% ylabel('ROIs');
% title('Matched and discarded')
% ylabel(colorbar, 'Activity (normalised)');
% 
% saveas(fig3, fullfile(['Z:\ibn-vision\USERS\Sonali\Figures\BoutonCorr\' sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' 'SpatialTuningPlots_SingleMatchedROIs_' dataScope '.png']))