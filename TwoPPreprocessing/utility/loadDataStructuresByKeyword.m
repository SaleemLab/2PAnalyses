function loadedData = loadDataByKeyword(sessionFileInfo, stimKeyword, varargin)
% Flexibly loads specific data for stimuli matching a keyword.
%
% This function finds all stimuli containing 'stimKeyword', excludes combined runs,
% and then loads only the data specified by the field names in varargin
% (e.g., 'Response', 'processedMergedBonsaiSuite2pData').
%
% NOTE: If only one field is requested, it returns the cell array directly.
% If multiple fields are requested, it returns a struct containing the cell arrays.
%
% SYNTAX:
%   loadedData = loadDataByKeyword(sessionFileInfo, stimKeyword, 'field1', 'field2', ...)
%
% EXAMPLE (Single Variable):
%   % Returns a cell array directly
%   responseCellArray = loadDataByKeyword(sessionFileInfo, 'VRCorr', 'Response');
%
% EXAMPLE (Multiple Variables):
%   % Returns a struct
%   dataStruct = loadDataByKeyword(sessionFileInfo, 'RFMapping', ...
%       'Response', 'processedMergedBonsaiSuite2pData');

%% Input Validation
if nargin < 3
    error('You must specify at least one sessionFileInfo field to load (e.g., ''Response'').');
end
fieldsToLoad = varargin;

%% Find All Valid Stimulus Runs
excludeKeywords = {'CombinedRuns'};
allStimNames = {sessionFileInfo.stimFiles.name};

indicesToProcess = find(contains(allStimNames, stimKeyword, 'IgnoreCase', true));
indicesToExclude = find(contains(allStimNames, excludeKeywords, 'IgnoreCase', true));
stimRunIndices = setdiff(indicesToProcess, indicesToExclude);

%% Initialize Output Struct
outputStruct = struct();
for i = 1:length(fieldsToLoad)
    outputStruct.(fieldsToLoad{i}) = {}; % e.g., outputStruct.Response = {};
end

if isempty(stimRunIndices)
    fprintf('No stimuli found matching keyword "%s". Returning empty.\n', stimKeyword);
    % Handle the output format for the no-data case
    if length(fieldsToLoad) == 1
        loadedData = {}; % Return empty cell if one thing was requested
    else
        loadedData = outputStruct; % Return empty struct if multiple things were requested
    end
    return;
end
fprintf('Found %d runs matching keyword "%s". Preparing to load data...\n', length(stimRunIndices), stimKeyword);

%% Loop Through Runs and Load Requested Data
for i = 1:length(stimRunIndices)
    runIndex = stimRunIndices(i);
    stimName = allStimNames{runIndex};
    
    tempRunData = struct(); % Temp storage for all data for this run
    isRunComplete = true;   % Flag to check if all requested files are found

    for j = 1:length(fieldsToLoad)
        fieldName = fieldsToLoad{j};
        
        if ~isfield(sessionFileInfo.stimFiles, fieldName)
            warning('Field "%s" does not exist in sessionFileInfo.stimFiles. Skipping.', fieldName);
            isRunComplete = false;
            break; 
        end
        
        filePath = sessionFileInfo.stimFiles(runIndex).(fieldName);

        if isempty(filePath) || ~exist(filePath, 'file')
            warning('File path for "%s" in stimulus "%s" not found or is empty. Skipping this entire run.', fieldName, stimName);
            isRunComplete = false;
            break;
        end
        
        switch fieldName
            case 'Response'
                varNameToLoad = 'response';
            case 'processedMergedBonsaiSuite2pData'
                varNameToLoad = 'processedTwoPData';
            otherwise
                varNameToLoad = fieldName;
        end
        
        try
            matObj = matfile(filePath);
            if isprop(matObj, varNameToLoad)
                tempRunData.(fieldName) = matObj.(varNameToLoad);
            else
                warning('Variable "%s" not found in file: %s. Skipping this entire run.', varNameToLoad, filePath);
                isRunComplete = false;
                break;
            end
        catch ME
            warning('Could not load from file %s. Error: %s. Skipping this entire run.', filePath, ME.message);
            isRunComplete = false;
            break;
        end
    end
    
    if isRunComplete
        fprintf('Successfully loaded all requested data for stimulus: %s\n', stimName);
        for j = 1:length(fieldsToLoad)
            fieldName = fieldsToLoad{j};
            outputStruct.(fieldName){end+1} = tempRunData.(fieldName);
        end
    end
end

fprintf('Finished. Loaded complete data for %d of %d runs.\n', length(outputStruct.(fieldsToLoad{1})), length(stimRunIndices));

%% Final Output Formatting (Modified Section)
% If the user only asked for one variable, return the cell array directly.
% Otherwise, return the full struct.
if length(fieldsToLoad) == 1
    loadedData = outputStruct.(fieldsToLoad{1});
else
    loadedData = outputStruct;
end

end