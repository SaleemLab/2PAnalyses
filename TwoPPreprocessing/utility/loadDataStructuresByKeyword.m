function loadedData = loadDataStructuresByKeyword(sessionFileInfo, stimKeyword, varargin)
% Flexibly loads specific data for stimuli matching a keyword.
% Updated Jan 2026 for Flattened .mat files.

%% Input Validation
if nargin < 3
    error('You must specify at least one sessionFileInfo field to load (e.g., ''Response'').');
end
fieldsToLoad = varargin;

%% Find All Valid Stimulus Runs
excludeKeywords = {'Combined'};
allStimNames = {sessionFileInfo.stimFiles.name};
indicesToProcess = find(contains(allStimNames, stimKeyword, 'IgnoreCase', true));
indicesToExclude = find(contains(allStimNames, excludeKeywords, 'IgnoreCase', true));
stimRunIndices = setdiff(indicesToProcess, indicesToExclude);

%% Initialize Output Struct
outputStruct = struct();
for i = 1:length(fieldsToLoad)
    outputStruct.(fieldsToLoad{i}) = {}; 
end

if isempty(stimRunIndices)
    fprintf('No stimuli found matching keyword "%s". Returning empty.\n', stimKeyword);
    if isscalar(fieldsToLoad), loadedData = {}; else, loadedData = outputStruct; end
    return;
end

fprintf('Found %d runs matching keyword "%s". Loading flattened variables...\n', length(stimRunIndices), stimKeyword);

%% Loop Through Runs and Load Requested Data
for i = 1:length(stimRunIndices)
    runIndex = stimRunIndices(i);
    stimName = allStimNames{runIndex};
    
    tempRunData = struct(); 
    isRunComplete = true;   
    
    for j = 1:length(fieldsToLoad)
        fieldName = fieldsToLoad{j};
        
        if ~isfield(sessionFileInfo.stimFiles, fieldName)
            warning('Field "%s" does not exist in sessionFileInfo.stimFiles. Skipping.', fieldName);
            isRunComplete = false;
            break; 
        end
        
        filePath = sessionFileInfo.stimFiles(runIndex).(fieldName);
        if isempty(filePath) || ~exist(filePath, 'file')
            warning('File path for "%s" in stimulus "%s" not found. Skipping run.', fieldName, stimName);
            isRunComplete = false;
            break;
        end
        
        try
            % Since files are now FLAT, load(filePath) pulls all root-level 
            % variables and packs them into a struct. This reconstructs 
            % 'response' or 'processedTwoPData' locally.
            tempRunData.(fieldName) = load(filePath);
            
        catch ME
            warning('Could not load from file %s. Error: %s. Skipping run.', filePath, ME.message);
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

%% Final Output Formatting
if isscalar(fieldsToLoad)
    loadedData = outputStruct.(fieldsToLoad{1});
else
    loadedData = outputStruct;
end

end