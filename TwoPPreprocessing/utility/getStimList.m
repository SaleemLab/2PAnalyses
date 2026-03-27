function stimList = getStimList(mouseID, session)
%% Directories 
%% Uses fall.mat to find the order or stimuli that were concatenated before processing through suite2p. 
rootDir = ['Z:' filesep fullfile('ibn-vision','DATA','SUBJECTS',mouseID)];
processedFolder = fullfile(rootDir, 'Processed', session);
suite2pFolder = fullfile(processedFolder, 'suite2p');

%% Load fall.mat to access the list of all TIFF files processed
planeItems = dir(fullfile(suite2pFolder, 'plane*'));
planeFolders = planeItems([planeItems.isdir]);
if isempty(planeFolders)
    error('No plane folder found in %s', suite2pFolder);
end
firstPlaneFolder = planeFolders(1).name; 
firstPlaneFall_FilePath = fullfile(suite2pFolder, firstPlaneFolder, 'Fall.mat');

% Load the .mat file to access the ops struct
fAll = load(firstPlaneFall_FilePath);
ops = fAll.ops;
if isfield(ops, 'filelist')
    tifFileList = ops.filelist; 
elseif isfield(ops, 'file_list')
    tifFileList = ops.file_list; 
end
%% Parse the file list to get the final, correctly formatted stimList
stimList = parseStimNames(tifFileList);
end


%% LOCAL HELPER FUNCTION
% This function now correctly identifies stimulus names by extracting the name
% of the parent folder for each TIFF file. This is robust to complex
% filenames with sequence numbers or '_trimmed' suffixes.
function uniqueNames = parseStimNames(fileList)
    if ischar(fileList), fileList = cellstr(fileList); end

    % Create a cell array to hold the parent folder names for each file
    parentFolderNames = cell(size(fileList));

    for i = 1:length(fileList)
        % Get the full path of the immediate parent directory of the TIFF file
        parentPath = fileparts(fileList{i});
        
        % From that directory path, extract just the folder name itself.
        % This name is the true stimulus run name.
        [~, folderName, ~] = fileparts(parentPath);
        
        parentFolderNames{i} = folderName;
    end

    % The unique list of these folder names is our desired stimList.
    % 'stable' preserves the order of appearance.
    uniqueNames = unique(parentFolderNames, 'stable')';
end