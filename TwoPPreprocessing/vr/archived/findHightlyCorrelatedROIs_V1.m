function [roisToKeep, roisToDiscard, groups, corrMatrix] = findHightlyCorrelatedROIs_V1(sessionFileInfo, signalToUse, useNonVRRecOnly, plotFlag)
% Finds highly correlated ROIs (e.g., from the same axon) and keeps only
% the one with the highest SNR (max z-score) from each group.
%
% OUTPUTS:
%   roisToKeep    - Column vector of indices for ROIs to KEEP.
%   roisToDiscard - Column vector of indices for ROIs to DISCARD.
%   groups        - Vector (1 x nROIs) mapping each original ROI to a
%                   component group ID.

% VRthreshold = 0.5 
% Set defaults
if nargin < 2; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 3; useNonVRRecOnly = true; end
if nargin < 3; plotFlag = true; end
% if nargin < 3; plotFlag = true; end % Uncomment for debugging
%%
% Find NonVR StimIdx and load response struture load processedTwoPData 
% Initialise the output variable
traces = [];
dataScope = '';

if useNonVRRecOnly
    % Grab indices 
    isNotVRStim = find(~contains({sessionFileInfo.stimFiles.name}, 'VRCorr'));
    if any(isNotVRStim)
        % Get the number of files found
        numNonVRFiles = length(isNotVRStim);
        allDFFData = cell(1, numNonVRFiles);
        for thisNonVRFile = 1:numNonVRFiles
            thisOtherStimIdx = isNotVRStim(thisNonVRFile);
            fprintf('Loading file %d/%d (Index: %d)...\n', thisNonVRFile, numNonVRFiles, thisOtherStimIdx);
            processedTwoPDataPath = sessionFileInfo.stimFiles(thisOtherStimIdx).processedMergedBonsaiSuite2pData;
            if exist(processedTwoPDataPath,"file")
                processedTwoPData = load(processedTwoPDataPath);
                allDFFData{thisNonVRFile} = processedTwoPData.processedTwoPData.zScoredProcessedSignals.(signalToUse);
            else
                fprintf('Missing processedTwoPData file: %s \n', sessionFileInfo.stimFiles(thisOtherStimIdx).name);
            end
        end
        % Concatenate all data
        traces = [allDFFData{:}];
        dataScope = 'NonVRStimOnly';
    end

% Find the right VRStimIdx and load processedTwoPData
else
     % Find all VR stim files that are NOT CombinedRuns
    VRStimIndices = find(contains({sessionFileInfo.stimFiles.name}, 'VRCorr') & ...
        ~contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );

    % Get the number of files found
    numFiles = length(VRStimIndices);
    if numFiles == 1
        % Exactly one file found
        fprintf('Loading single processedTwoPData VRCorr file (Index: %d)...\n', VRStimIndices);

        processedTwoPDataPath = sessionFileInfo.stimFiles(VRStimIndices).processedMergedBonsaiSuite2pData;
        processedTwoPData = load(processedTwoPDataPath);
        % Use dynamic field name to load the correct signal
        traces = processedTwoPData.processedTwoPData.zScoredProcessedSignals.(signalToUse);
        dataScope = 'VRStim';

    elseif numFiles > 1
        % Multiple files found (Concatenate)
        fprintf('Multiple processedTwoPData VRCorr files found. Loading and concatenating %d ...\n', numFiles);

        allDFFData = cell(1, numFiles);
        for thisVRFile = 1:numFiles
            thisStimIdx = VRStimIndices(thisVRFile);
            fprintf('Loading file %d/%d (Index: %d)...\n', thisVRFile, numFiles, thisStimIdx);
            processedTwoPDataPath = sessionFileInfo.stimFiles(thisStimIdx).processedMergedBonsaiSuite2pData;
            processedTwoPData = load(processedTwoPDataPath);
            allDFFData{thisVRFile} = processedTwoPData.processedTwoPData.zScoredProcessedSignals.(signalToUse);
        end
        % Concatenate all data
        traces = [allDFFData{:}];
        dataScope = 'VRStimRuns';
    else
        error('No VRCorr files found.\n');
    end
end
% Check if traces were loaded
if isempty(traces)
    error('Traces variable is empty after loading. Check file paths and signalToUse name.');
end

[nROIs, ~] = size(traces);
fprintf('Loaded %d ROIs. Starting correlation analysis...\n', nROIs);

%% Calculate SNR Score for all ROIs
% Using max z-score as the "signal-to-noise" metric 
fprintf('Calculating SNR scores (max z-score) for %d ROIs...\n', nROIs);
snrScores = max(traces, [], 2); % Max value of each row (dim 2)

%% Prepare Traces for Correlation
% Data is assumed to be pre-smoothed. Using traces directly.
fprintf('Using loaded traces directly for correlation...\n');
tracesToCorrelate = traces;

%% Calculate Pairwise Correlation
fprintf('Calculating %dx%d correlation matrix...\n', nROIs, nROIs);
% corr() expects columns as variables, so we transpose our (ROI x Time) matrix
% to (Time x ROI)
corrMatrix = corr(tracesToCorrelate');

%% Find Groups (Connected Components)
% Create an adjacency matrix based on the 0.5 threshold
correlationThreshold = 0.3;
% Filtering highly correlated ROIs based on threshold. 
adjMatrix = corrMatrix > correlationThreshold;
% An ROI should not be 'grouped' with itself; 
% Matrices across the diagonal are 0. 
% Nodes cannot have connections to themselves i.e., looped.  
adjMatrix(logical(eye(nROIs))) = false; % Set diagonal to false; eye - creates an identity matrix 

% Use graph theory to find all connected components
fprintf('Finding correlated groups (rho > %.2f)...\n', correlationThreshold);
% % Would be the same if selecting 'lower' because adjacency matrices are symmetrical across the diagonal 
G = graph(adjMatrix, 'upper'); 
  
groups = conncomp(G); % Assigns a group ID to each ROI
uniqueGroups = unique(groups);
fprintf('Found %d unique groups (including singletons).\n', length(uniqueGroups));

%% Keep only highest SNR ROI from each group
fprintf('Pruning groups to keep highest SNR ROI...\n');

% Logical masks for ROIs to keep vs. discard
roisToKeep_mask = false(nROIs, 1);
roisToDiscard_mask = false(nROIs, 1);

for thisGroup = 1:length(uniqueGroups)
    thisGroupID = uniqueGroups(thisGroup);
    
    % Get the global indices of all ROIs in this group
    roisInGroup_idx = find(groups == thisGroupID);
    
    if isscalar(roisInGroup_idx)
        % This is a single (uncorrelated ROI) 
        % Keep it
        roisToKeep_mask(roisInGroup_idx) = true;
    else
        % This is a correlated group (2 or more ROIs) 
        % Find the winner (highest SNR) within this group
        groupSnrScores = snrScores(roisInGroup_idx);
        [~, localWinnerIdx] = max(groupSnrScores); % Max zscoreddFFNeuroCorr @Aman? 
        % Get the global index of the winner
        globalWinnerIdx = roisInGroup_idx(localWinnerIdx);
        % Keep only the winner
        roisToKeep_mask(globalWinnerIdx) = true;
        % Mark all the other ROIs in this group as the discarded group. 
        roisInGroup_idx(localWinnerIdx) = []; % remove winner from list!! 
        roisToDiscard_mask(roisInGroup_idx) = true;
    end
end

% Convert logical masks to roi index numbers 
roisToKeep = find(roisToKeep_mask);
roisToDiscard = find(roisToDiscard_mask);

fprintf('Analysis complete. Kept %d ROIs, discarded %d ROIs.\n', ...
        length(roisToKeep), length(roisToDiscard));
        
%% Plotting (Optional)
if plotFlag
    fig1 = figure('Name', 'Correlation Matrix and Adjacency Matrices', 'NumberTitle', 'off');
    
    % Full correlation matrix
    subplot(1, 2, 1);
    imagesc(corrMatrix);
    axis image;
    % colormap gray;
    colorbar;
    title(sprintf('Full Correlation Matrix (%d ROIs)', nROIs));
    xlabel('ROI Index');
    ylabel('ROI Index');
    
    % Thresholded adjacency matrix
    % This would mean that dots in black (0) do not have connections. 
    % Correlation between rois are below threshold 

    subplot(1, 2, 2);
    imagesc(adjMatrix);
    colormap gray;
    axis image;
    title(sprintf('Adjacency/Connectivity Matrix (r > %d)', correlationThreshold));
    xlabel('ROI Index');
    ylabel('ROI Index'); 
    saveas(fig1, fullfile(['Z:\ibn-vision\USERS\Sonali\Figures\BoutonCorr\' sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' 'CorrAndAdjacencyMatrix_' dataScope '.png']))

    
    % Why are the dff values so high? CHECK 
    discardedGroupIDs = groups(roisToDiscard);
    uniqueDiscardedGroups = unique(discardedGroupIDs);

    if numel(uniqueDiscardedGroups) >= 2
    
        fig2 = figure('Name', 'Correlated ROI Groups', 'NumberTitle', 'off');
        subplot(211)
        gM1 = find(groups == uniqueDiscardedGroups(1)); 
        plot(traces(gM1, :)');
        title(['Discarded Group ID: ' num2str(uniqueDiscardedGroups(1))]);
        ylabel('dF/F');
    
    
        subplot(212)
        gM2 = find(groups == uniqueDiscardedGroups(2)); 
        plot(traces(gM2, :)');
        title(['Discarded Group ID: ' num2str(uniqueDiscardedGroups(2))]);
        ylabel('dF/F');
        saveas(fig2, fullfile(['Z:\ibn-vision\USERS\Sonali\Figures\BoutonCorr\' sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' 'HighlyCorreltedROIGroups_' dataScope '.png',]))

    
    else
        warning('Could not find 2 unique groups in roisToDiscard to plot.');
    end
    % Plot spatial tuning curves
    % Change in the morning: use idex to load either the second response
    % file or combined if present. 
    load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25040\Analysis\20250511B\M25040_20250511B_Response_M25040_VRCorr_20250511_00002.mat")
    lapActivity = response.lapPositionActivity.dFFNeuropilCorrected; 
    meanAll = squeeze(mean(lapActivity, 2, 'omitnan'));
    normAll = normalize(meanAll, 2, 'range');
    [~, peakIdx] = max(normAll, [], 2);
    [~, sortIdx] = sort(peakIdx);
    
    meanROIsToKeep = squeeze(mean(lapActivity(roisToKeep,:,:), 2, 'omitnan'));
    normROIsToKeep = normalize(meanROIsToKeep, 2, 'range');
    [~, peakIdxToKeep] = max(normROIsToKeep, [], 2);
    [~, sortIdxToKeep] = sort(peakIdxToKeep);

    meanROIsToDiscard = squeeze(mean(lapActivity(roisToDiscard,:,:), 2, 'omitnan'));
    normROIsToDiscard = normalize(meanROIsToDiscard, 2, 'range');
    [~, peakIdxToDiscard] = max(normROIsToDiscard, [], 2);
    [~, sortIdxToDiscard] = sort(peakIdxToDiscard);

    fig3 = figure;
    subplot(131)
    imagesc(normAll(sortIdx, :));
    xlabel('Position (cm)');
    ylabel('ROIs');
    title('All')
    % caxis([0 1]); colormap(flipud(gray));
    ylabel(colorbar, 'Activity (normalised)');

    subplot(132)
    imagesc(normROIsToKeep(sortIdxToKeep, :));
    xlabel('Position (cm)');
    ylabel('ROIs');
    title('Single pringle')
    % caxis([0 1]); colormap(flipud(gray));
    ylabel(colorbar, 'Activity (normalised)');

    subplot(133)
    imagesc(normROIsToDiscard(sortIdxToDiscard, :));
    xlabel('Position (cm)');
    ylabel('ROIs');
    title('Matched and discarded')
    % caxis([0 1]); colormap(flipud(gray));
    ylabel(colorbar, 'Activity (normalised)');
    saveas(fig3, fullfile(['Z:\ibn-vision\USERS\Sonali\Figures\BoutonCorr\' sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_' 'SpatialTuningPlots_SingleMatchedROIs_' dataScope '.png']))
end

end


% [roisToKeepNonVR, roisToDiscardNonVR, groupsNonVR, corrMatrixNonVR] = findHightlyCorrelatedROIs(sessionFileInfo, 'dFFNeuropilCorrected', true, true);
% [roisToKeepVR, roisToDiscardVR, groupsVR, corrMatrixVR] = findHightlyCorrelatedROIs(sessionFileInfo, 'dFFNeuropilCorrected', false, true);
% idx = tril(true(size(corrMatrixNonVR)), -1);
% corrVR = corrMatrixVR(idx);
% corrNonVR = corrMatrixNonVR(idx);
% figure; scatter(corrVR, corrNonVR, 50, 'filled', 'MarkerFaceAlpha', 0.6);
% xlabel('CorrVR');
% ylabel('CorrNonVR');
% axis equal
% figure; scatter(corrVR, corrNonVR, 50, 'filled', 'MarkerFaceAlpha', 0.6);
% xlabel('CorrVR');
% ylabel('CorrNonVR');
% axis equal