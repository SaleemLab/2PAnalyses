function [sessionFileInfo] = createSessionROIData(sessionFileInfo)
%% Create new file
% Constructs the file name and path for the session ROI data file
fileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_sessionROIData.mat'];
sessionFileInfo.otherSessFilePaths.sessionROIData = fullfile(sessionFileInfo.Directories.save_folder, fileName);

%% Save imagedType from the masterdata sheet for each access
% Assumes filterMasterTable is an external, accessible function
filteredTable = filterMasterTable('MouseID', sessionFileInfo.animal_name, 'Session', sessionFileInfo.session_name);
typeImaged = filteredTable.TypeImaged;
targetArea = filteredTable.TargetArea;

%% Load Data
% Load 2PMetaData from one session; Using the VR session here
isVRStim = contains({sessionFileInfo.stimFiles.name}, 'VRCorr');
iStim = find(isVRStim==1);
if length(iStim) > 1
    iStim = iStim(1); % Take the first match if multiple VRStim files exist
end
disp('Loading TwoPMetaData and TwoPDara structs...');
load(sessionFileInfo.stimFiles(iStim).TwoPMetaData, "twopMetadata");
load(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, "twoPData");

%% TwoPMetaData
% Load information needed from twoPMetadata
scanZoomFactor = twopMetadata.scanZoomFactor;
if scanZoomFactor == 6
    micronsPerPix = 0.1318287;
    selectedPlane = twoPData.ops.selected_plane;
elseif scanZoomFactor == 4
    micronsPerPix = 0.0988692;
    selectedPlane = [];
else
    fprintf('Microns to pixel conversion missing for zoom: %d \n', scanZoomFactor);
    micronsPerPix = [];
    selectedPlane = [];
end
pixelsPerLine = twopMetadata.pixelsPerLine;
linesPerFrame = twopMetadata.linesPerFrame;

%% TwoPData (Info about planes, iscell, npixels location in FOV and radius)
% Vectors (Concatenated as column vectors [;])
roiGlobalIdx     = [];   % Added variable to store ROI index within its plane
roiPlaneIdentity = [];
iscell           = [];
redcell          = [];
radiuspix        = [];
npix             = [];

% cells (Concatenated as cell arrays)
% column cell arrays {0x1}
ypix = cell(0, 1);
xpix = cell(0, 1);

numPlanes = numel(twoPData);
for thisPlaneIdx = 1:numPlanes
    currentPlane = twoPData(thisPlaneIdx);
    numROIs = size(currentPlane.F, 1);
    
    roiIdx_currentPlane = (1:numROIs)';         % Create column vector [1; 2; ... N]
    roiGlobalIdx = [roiGlobalIdx; roiIdx_currentPlane]; % Concatenate to global index
    

    iscell            = [iscell; currentPlane.iscell(:,1)];
    redcell           = [redcell; currentPlane.redcell(:,1)];
    
    % --- Plane Identity and Stats ---
    roiPlaneIdentity  = [roiPlaneIdentity; repmat(thisPlaneIdx - 1, numROIs, 1)];
    allStatStructuresCell = currentPlane.stat;
    
    radiiVector = cellfun(@(c) c.radius, allStatStructuresCell, 'UniformOutput', true);
    npixVector  = cellfun(@(c) c.npix, allStatStructuresCell, 'UniformOutput', true);
    radiuspix         = [radiuspix; radiiVector'];
    npix              = [npix; npixVector'];
    
    ypixCellRow = cellfun(@(c) c.ypix(:), allStatStructuresCell, 'UniformOutput', false);
    xpixCellRow = cellfun(@(c) c.xpix(:), allStatStructuresCell, 'UniformOutput', false);
    ypix              = [ypix; ypixCellRow'];
    xpix              = [xpix; xpixCellRow'];
end

%% Assemble and Save sessionROIData Structure
sessionROIData.typeImaged        = typeImaged;
sessionROIData.targetArea        = targetArea;
sessionROIData.scanZoomFactor    = scanZoomFactor;
sessionROIData.micronsPerPix     = micronsPerPix;
sessionROIData.pixelsPerLine     = pixelsPerLine;
sessionROIData.linesPerFrame     = linesPerFrame;
sessionROIData.roiPlaneIdentity  = roiPlaneIdentity;
sessionROIData.roiGlobalIdx      = roiGlobalIdx;   % New field
sessionROIData.iscell            = iscell;
sessionROIData.redcell           = redcell;
sessionROIData.radiuspix         = radiuspix;
sessionROIData.npix              = npix;
sessionROIData.ypix              = ypix;
sessionROIData.xpix              = xpix;
sessionROIData.selectedPlane     = selectedPlane;

% Save individual variables to the .mat file for selective loading
save(sessionFileInfo.otherSessFilePaths.sessionROIData, ...
     'scanZoomFactor', 'micronsPerPix', 'pixelsPerLine', 'linesPerFrame', ...
     'typeImaged', 'targetArea', ...
     'roiPlaneIdentity', 'roiGlobalIdx', 'iscell', 'redcell', 'radiuspix', 'npix', ... % Added roiGlobalIdx
     'ypix', 'xpix', 'selectedPlane');
     
% Save sessionFileInfo
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end