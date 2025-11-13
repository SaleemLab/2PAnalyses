function [ResponseVRRuns, processed2PDataVRRuns] = checkForAndGetMultiVRRuns(sessionFileInfo)
% Finds and loads multiple VR 'response' and 'processed2PData' files.
%
% 1. Searches for stimulus files containing 'VRCorr' to identify VR runs.
% 2. Excludes any previously combined runs.
% 3. Checks that BOTH the response and processed2PData files exist for a run.
% 4. If more than one valid run is found, it loads the corresponding files.

%% Initialise and Find Runs
% Keywords to identify stimulus types and runs to exclude
vrKeywords = {'VRCorr'};
excludeKeywords = {'CombinedRuns'}; % Prevent merging already combined VR runs

% Get all stimulus names from the sessionFileInfo struct
allStimNames = {sessionFileInfo.stimFiles.name};

% Find the indices of the stimulus files to process
allVrIndices = find(contains(allStimNames, vrKeywords, 'IgnoreCase', true));
indicesToExclude = find(contains(allStimNames, excludeKeywords, 'IgnoreCase', true));
vrRunIndices = setdiff(allVrIndices, indicesToExclude);

% Initialize empty cell arrays for the outputs
ResponseVRRuns = {};
processed2PDataVRRuns = {};

%% Check number of VR runs before proceeding
if length(vrRunIndices) <= 1
    disp('Found 1 or 0 VRCorr runs to process. No combination necessary.');
    return; % Exit the function, returning the empty cell arrays
end

fprintf('Found %d VRCorr runs. Preparing to load files...\n', length(vrRunIndices));

%% Load data for each valid VR run
% Pre-allocate cell arrays for speed
ResponseVRRuns = cell(1, length(vrRunIndices));
processed2PDataVRRuns = cell(1, length(vrRunIndices));

for i = 1:length(vrRunIndices)
    % Get the actual index in the main sessionFileInfo.stimFiles array
    runIndex = vrRunIndices(i);
    stimName = allStimNames{runIndex};
    
    % Get the file paths for the current run
    responseFile = sessionFileInfo.stimFiles(runIndex).Response;
    processedFile = sessionFileInfo.stimFiles(runIndex).processedMergedBonsaiSuite2pData;
    
    % --- Robustness Check: Ensure both files exist before loading ---
    if isempty(responseFile) || ~exist(responseFile, 'file')
        warning('Response file not found or path is empty for stimulus "%s". Skipping this run.', stimName);
        continue; % continue skips to the next iteration of the loop
    end
    if isempty(processedFile) || ~exist(processedFile, 'file')
        warning('processedMergedBonsaiSuite2pData file not found or path is empty for stimulus "%s". Skipping this run.', stimName);
        continue;
    end
    
    % --- Load Response Data ---
    fprintf('Loading Response for: %s\n', stimName);
    loadedData = load(responseFile);
    if isfield(loadedData, 'response')
        ResponseVRRuns{i} = loadedData.response;
    else
        warning('Variable "response" not found in file %s. Skipping.', responseFile);
        continue;
    end

    % --- Load processedTwoPData ---
    fprintf('Loading processed2PData for: %s\n', stimName);
    loadedData = load(processedFile);
    if isfield(loadedData, 'processedTwoPData')
        processed2PDataVRRuns{i} = loadedData.processedTwoPData;
    else
        warning('Variable "processedTwoPData" not found in file %s. Skipping.', processedFile);
        % Invalidate the corresponding response data to keep the pairs consistent
        ResponseVRRuns{i} = [];
        continue;
    end
end

%% Final Cleanup
% Remove any empty cells from runs that were skipped due to errors
emptyCells = cellfun(@isempty, ResponseVRRuns);
ResponseVRRuns(emptyCells) = [];
processed2PDataVRRuns(emptyCells) = [];

if isempty(ResponseVRRuns)
    disp('No valid pairs of VR run files were loaded.');
else
    fprintf('Successfully loaded %d pairs of VR run files.\n', length(ResponseVRRuns));
end

end