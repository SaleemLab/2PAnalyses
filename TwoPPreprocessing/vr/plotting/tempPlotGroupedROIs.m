% NOTE: This code assumes the variables 'roisToKeep', 'roisToDiscard', and 'groups' 
% containing the final global ROI indices/groups from your analysis are currently 
% available in your MATLAB workspace.
% It also assumes the external function 'plot_discarded_groups' is available.
%
% CRITICAL ASSUMPTION: ALL ROIs in sessionROIData belong to the selectedPlane, 
% AND a full volume cycle is only 8 directories long (no separate channels).
% =========================================================================
% 1. DEFINE FILE PATHS AND PARAMETERS 
% =========================================================================
roiDataPath = "\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25040\Analysis\20250511B\M25040_20250511B_sessionROIData.mat"; % Path to your saved sessionROIData file
tiffFilePath ="\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25040\OPhys\20250511B\M25040_VRCorr_20250511_00002\M25040_VRCorr_20250511_00002_00004 - Copy.tif"; % Path to your raw TIFF file
startFrame = 10;   % Frame index to start averaging (1-based volume index)
numFramesToAverage = 1000; % Number of volumes to average for background image
numPlanes = 8; % Total number of planes in the TIF acquisition (FULL_CYCLE_SIZE is now 8)

% =========================================================================
% 2. LOAD ROI DATA
% =========================================================================
fprintf('Loading aggregated ROI data from: %s\n', roiDataPath);
try
    sessionROIData = load(roiDataPath);
catch ME
    error('Failed to load data. Check file path/name. Error: %s', ME.message);
end

% **CRITICAL CHECK: Ensure selectedPlane exists before use.**
if ~isfield(sessionROIData, 'selectedPlane')
    error('Required field sessionROIData.selectedPlane is missing from the loaded file.');
end

% Check for necessary fields
if ~all(isfield(sessionROIData, {'roiPlaneIdentity', 'ypix', 'xpix'}))
    error('sessionROIData structure is missing required fields (roiPlaneIdentity, ypix, or xpix).');
end

% Define the target plane index (used only for TIF reading in Section 4).
targetPlaneIdx = sessionROIData.selectedPlane; 

% =========================================================================
% 3. FILTER ROIS TO THE TARGET PLANE
% =========================================================================
planeROIsToKeep = intersect(roisToKeep, sessionROIData.roiGlobalIdx);
planeROIsToDiscard = intersect(roisToDiscard, sessionROIData.roiGlobalIdx);
fprintf('Found %d ROIs to Keep and %d ROIs to Discard on Plane %d.\n', ...
         length(planeROIsToKeep), length(planeROIsToDiscard), targetPlaneIdx);
         
% =========================================================================
% 4. LOAD AND PROCESS TIFF IMAGE FOR GREEN CHANNEL BACKGROUND (Type Casting Fix)
% =========================================================================
backgroundImg = [];
try
    t = Tiff(tiffFilePath, 'r');
    
    % Get dimensions and file info
    fileInfo = imfinfo(tiffFilePath);
    nTotalFrames = numel(fileInfo);
    imgH = fileInfo(1).Height;
    imgW = fileInfo(1).Width;
    
    % *** FIX 1: The full volume cycle is 16 directories (8 planes * 2 channels) ***
    FULL_CYCLE_SIZE = double(numPlanes * 2); 

    % 1. Calculate the offset (initial directory) for the Green channel
    % FIX 2: Offset is based on 2 directories per plane (2 * plane_index + 1 for Green)
    targetOffset = (double(targetPlaneIdx) * 2) + 1; 

    % 2. Determine the volume range (0-indexed volumes)
    startVolume = double(startFrame) - 1; 
    volumesToRead = startVolume : (startVolume + numFramesToAverage - 1); 
    
    % 3. Calculate the final 1-based directory indices 
    % The addition/multiplication is now safe due to the double() casting.
    frameIndicesToRead = targetOffset + volumesToRead * FULL_CYCLE_SIZE;
    frameIndicesToRead = frameIndicesToRead(frameIndicesToRead <= nTotalFrames);
    
    % Preallocate background image
    backgroundImg = zeros(imgH, imgW, 'single');
    frameCounter = 0;
    
    % Load and sum the frames using Tiff object
    for k = 1:length(frameIndicesToRead)
        currentDirectory = frameIndicesToRead(k);
        
        try
            t.setDirectory(currentDirectory);
            frameData = t.read(); 

            if ~isempty(frameData) && all(size(frameData) == [imgH, imgW])
                backgroundImg = backgroundImg + single(frameData); 
                frameCounter = frameCounter + 1;
            else
                warning('Skipping directory %d (Frame %d). Data size mismatch.', currentDirectory, k);
            end
        catch readErr
            warning('Skipping directory %d (Frame %d). Read failed: %s', currentDirectory, k, readErr.message);
            continue; 
        end
    end
    t.close();
    
    if frameCounter > 0
        backgroundImg = backgroundImg / frameCounter; 
        fprintf('Averaged %d Green Channel frames for Plane %d background image.\n', frameCounter, targetPlaneIdx);
    else
        error('Could not load any valid frames for the Green channel background image.');
    end
    
catch ME
    warning('Error loading or processing TIFF file: %s. Using default black background.', ME.message);
    backgroundImg = zeros(256, 256); 
    imgH = 256; 
    imgW = 256;
end
% =========================================================================
% 5. VISUAL CONFIRMATION: PLOT DISCARDED GROUPS ONLY
% =========================================================================
% Define image dimensions (Updated if default was used)
imgH = size(backgroundImg, 1);
imgW = size(backgroundImg, 2);

% SIMPLIFICATION: Assuming all loaded ROIs belong to this plane, 
% we use the whole ypix/xpix arrays directly without filtering.
planeYpixAll = sessionROIData.ypix;
planeXpixAll = sessionROIData.xpix;

figure('Name', sprintf('Discarded Correlated ROIs on Plane %d', targetPlaneIdx));
imagesc(backgroundImg);
colormap('gray');
axis image;
hold on;
% Call the external function: Plots ONLY the removed ROIs, color-coded by group ID.
plot_discarded_groups(roisToDiscard, groups, ...
                      sessionROIData.roiGlobalIdx, planeYpixAll, planeXpixAll, imgH, imgW); 
title(sprintf('Plane %d: Discarded Correlated ROIs (Color-coded by Group ID)', targetPlaneIdx));
hold off;