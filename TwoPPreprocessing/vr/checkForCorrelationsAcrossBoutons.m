function checkForCorrelationsAcrossBoutons(mouseID, session, checkVRRecOnly, onlyIncludeROIs, targetArea)

% How similar is the activity trace of ROI #1 to the activity trace of 
% ROI #2 over the entire time-point of recording or VRSession?
rootDir = 'Z:\ibn-vision\DATA\SUBJECTS';
fileName = sprintf('%s_%s_sessionFileInfo.mat', mouseID, session);
sessionFilePath = fullfile(rootDir, mouseID, 'Analysis', session, fileName);

% Load SessionFileInfo 
load(sessionFilePath);

if checkVRRecOnly
    fprintf('checkVRRecOnly=true. Searching for VRCorr stim files...\n');
    isVRStim = find( contains({sessionFileInfo.stimFiles.name}, 'VRCorr') & ~contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns') );
    if length(isVRStim) > 1
        % Select the last one 
        isVRStim = isVRStim(end);
    else 
        isVRStim = isVRStim; 
    end 

    fprintf('Multiple files found. Selecting last one (Index: %d).\n', isVRStim);
    % Load for VR stimulus 
    twoPDataPath = sessionFileInfo.stimFiles(isVRStim).mergedBonsai2PSuite2pData;
    % Loading without assignment 
    fprintf('Loading 2P data file: %s\n', twoPDataPath);
    load(twoPDataPath)
    isROI = twoPData.iscell;
    F = twoPData.F; 
    
else
    fprintf('checkVRRecOnly=false. Finding and loading Fall.mat...\n');
    FallPath = findFile(sessionFileInfo.suite2pFiles.planes, 'Fall.mat');
    fprintf('Loading Fall.mat file: %s\n', FallPath);
    fAll = load(FallPath);
    isROI = fAll.iscell;
    F = fAll.F; 
end

% Select ROIs
if onlyIncludeROIs
    % Find indices where the first column of iscell is 1
    ROIs = find(isROI(:, 1) == 1);
    fprintf('Including %d ROIs (iscell == 1)\n', length(ROIs));
else
    % Include all potential ROIs
    ROIs = 1:size(isROI, 1);
    fprintf('Including all %d potential ROIs\n', length(ROIs));
end

numCells = length(ROIs);
if numCells < 2
    fprintf('Fewer than 2 ROIs found. Cannot calculate correlations.\n');
    return;
end

% Extract F
% F contains fluorescence traces (ROIs x Time)
% We select only the rows corresponding to our chosen ROIs
fSelected = F(ROIs, :);

% Calculate Pairwise Correlations 
% corrcoef calculates the correlation matrix. 
% Each row of F_selected is a variable (an ROI's trace), 
% and each column is time point.
% Transpose fSelected so that ROIs are columns (variables).
% R will be a numCells x numCells matrix 
fprintf('Calculating %d x %d correlation matrix...\n', numCells, numCells);
R = corrcoef(fSelected');

% Find distribution of unique pairwise correlations.
% Create a mask for the upper triangle (k=1 means exclude the main diagonal)
upperTriangleMask = triu(true(numCells), 1);

% Extract the correlation values using the mask
correlationValues = R(upperTriangleMask);

fprintf('Extracted %d unique pairwise correlation values.\n', length(correlationValues));


% Plot
hFig = figure; 
histogram(correlationValues, 100); 
if checkVRRecOnly
    title(sprintf('Distribution of Pairwise ROI Correlations on VRStimulusOnly\n %s Mouse: %s, Session: %s (n=%d ROIs)', ...
              targetArea, mouseID, session, numCells));
else 
    title(sprintf('Distribution of Pairwise ROI Correlations on concatenated stimuli\n %s Mouse: %s, Session: %s (n=%d ROIs)', ...
              targetArea, mouseID, session, numCells));
end
xlabel('Pearson Correlation Coefficient (r)');
ylabel('Count');
grid on;
box off;


fprintf('Saving figure...\n');

savePath = fullfile(sessionFileInfo.Directories.save_folder, 'Figures');
if ~exist(savePath, 'dir')
    mkdir(savePath);
    fprintf('Created save directory: %s\n', savePath);
end

if checkVRRecOnly
    fileName = [mouseID '_' session '_BoutonCorrelations_OnVRStimulusOnly.png'];
else 
    fileName = [mouseID '_' session '_BoutonCorrelations_OnConcatenatedStimuli.png'];
end
fullSavePath = fullfile(savePath, fileName);
saveas(hFig, fullSavePath);
fprintf('Figure saved to: %s\n', fullSavePath);

end