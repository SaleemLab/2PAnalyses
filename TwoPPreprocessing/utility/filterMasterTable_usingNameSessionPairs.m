function filteredTable = filterMasterTable_usingNameSessionPairs(varargin)
% Filters the master data table based on specified criteria.
% It loads and parses the table using a nested function and applies
% flexible, name-value pair filters. Supports paired Mouse/Session tracking.
%  % Define pairs using a struct
%   pairs.M25032 = [1, 2, 3];
%   pairs.M26001 = [1, 3, 5];
%ExpRSPSessions = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);
% Define the path and load the parsed master table using the nested function
filepath = "Z:\ibn-vision\USERS\Sonali\datatable\MasterExpDatatable.csv";
masterTable = loadAndParseMasterTable(filepath);

p = inputParser;
p.KeepUnmatched = true;

% Define default values
defaultMousePairs = []; % Expects a struct or cell array mapping mice to sessions
defaultMouseID = string.empty;
defaultSession = string.empty;
defaultDay = [];
defaultHasStimulus = [];
defaultExclude = [];

% Add the parameters
addParameter(p, 'MousePairs', defaultMousePairs);
addParameter(p, 'MouseID', defaultMouseID);
addParameter(p, 'Session', defaultSession);
addParameter(p, 'DayOfExperience', defaultDay, @isnumeric);
addParameter(p, 'HasStimulus', defaultHasStimulus);
addParameter(p, 'Exclude', defaultExclude, @(x) isscalar(x) && (x==0 || x==1));

% Parse the inputs
parse(p, varargin{:});
criteria = p.Results;
unmatchedCriteria = p.Unmatched;

%% Start with an index that includes all rows ---
rowIndex = true(height(masterTable), 1);

%% --- Apply Paired Mouse-Session Filter if Provided ---
if ~isempty(criteria.MousePairs)
    % Create an initial mask of all zeros for the pairing filter
    pairIndex = false(height(masterTable), 1);

    if isstruct(criteria.MousePairs)
        mice = fieldnames(criteria.MousePairs);
        for m = 1:length(mice)
            currentMouse = mice{m};
            rawSessions = criteria.MousePairs.(currentMouse);

            % --- PARSING FIX FOR CONCATENATED 8-DIGIT DATE STRINGS ---
            % If it's a single string/char array longer than 8 characters, and looks like joined dates
            if (ischar(rawSessions) || isstring(rawSessions)) && strlength(rawSessions) > 8
                sessionStr = string(rawSessions);
                % Use regexp to slice the string into chunks of exactly 8 digits
                allowedSessions = regexp(sessionStr, '\d{8}', 'match');
            else
                allowedSessions = string(rawSessions);
            end

            % Find rows that match BOTH this specific mouse AND its parsed sessions
            mouseMask = (string(masterTable.MouseID) == string(currentMouse)) & ...
                ismember(string(masterTable.Session), allowedSessions);

            pairIndex = pairIndex | mouseMask;
        end
    elseif iscell(criteria.MousePairs)
        for r = 1:size(criteria.MousePairs, 1)
            currentMouse = criteria.MousePairs{r, 1};
            rawSessions = criteria.MousePairs{r, 2};

            % Apply the same parsing fix for cell array formats
            if (ischar(rawSessions) || isstring(rawSessions)) && strlength(rawSessions) > 8
                sessionStr = string(rawSessions);
                allowedSessions = regexp(sessionStr, '\d{8}', 'match');
            else
                allowedSessions = string(rawSessions);
            end

            mouseMask = (string(masterTable.MouseID) == string(currentMouse)) & ...
                ismember(string(masterTable.Session), allowedSessions);

            pairIndex = pairIndex | mouseMask;
        end
    end

    rowIndex = rowIndex & pairIndex;
else
    %% --- Fallback to Independent Filters if MousePairs is not used ---
    if ~isempty(criteria.MouseID)
        rowIndex = rowIndex & ismember(string(masterTable.MouseID), string(criteria.MouseID));
    end

    if ~isempty(criteria.Session)
        rowIndex = rowIndex & ismember(string(masterTable.Session), string(criteria.Session));
    end
end

%% --- Apply Remaining Filters ---
% Filter by DayOfExperience
if ~isempty(criteria.DayOfExperience)
    rowIndex = rowIndex & ismember(masterTable.DayOfExperience, criteria.DayOfExperience);
end

% Filter by presence of data in a stimulus column
if ~isempty(criteria.HasStimulus)
    stimNamesToCheck = cellstr(criteria.HasStimulus);
    validStimCols = stimNamesToCheck(ismember(stimNamesToCheck, masterTable.Properties.VariableNames));

    if ~isempty(validStimCols)
        hasStimIndex = false(height(masterTable), 1);
        for i = 1:height(masterTable)
            for j = 1:length(validStimCols)
                colName = validStimCols{j};
                if ~isempty(masterTable.(colName){i})
                    hasStimIndex(i) = true;
                    break;
                end
            end
        end
        rowIndex = rowIndex & hasStimIndex;
    end
end

% Filter by Exclude flag
if ~isempty(criteria.Exclude)
    rowIndex = rowIndex & (masterTable.Exclude == criteria.Exclude);
end

% Apply generic filters for any other table columns provided
otherFieldNames = fieldnames(unmatchedCriteria);
for i = 1:length(otherFieldNames)
    fieldName = otherFieldNames{i};
    if ismember(fieldName, masterTable.Properties.VariableNames)
        filterValue = unmatchedCriteria.(fieldName);
        rowIndex = rowIndex & ismember(string(masterTable.(fieldName)), string(filterValue));
    else
        warning('Ignoring unrecognized parameter: ''%s'' is not a valid column name.', fieldName);
    end
end

%% Create the final filtered table ---
filteredTable = masterTable(rowIndex, :);
disp(['Filtering complete. ' num2str(height(filteredTable)) ' rows selected.']);
end

%% Nested function to load and parse mastertable
function masterTable = loadAndParseMasterTable(filepath)
opts = detectImportOptions(filepath, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'string');
masterTable = readtable(filepath, opts);

columnsToParse = {
    'VRCorr', 'BaselineCorridor', 'LandManipCorridor' ,'RFMapping', 'GrayScreen', 'DirTuning', ...
    'CheckerBoard', 'DotMotion_SpeedTuning', 'DriftingBar', ...
    'FullFiledFlash', 'SparseNoiseTexture', 'zStack', 'LapsRecorded', 'LapsCompleted'
    };
columnsToParse = columnsToParse(ismember(columnsToParse, masterTable.Properties.VariableNames));

for i = 1:length(columnsToParse)
    colName = columnsToParse{i};
    newColumn = cell(height(masterTable), 1);
    originalColumnData = masterTable.(colName);

    for j = 1:height(masterTable)
        cellContent = originalColumnData(j);
        if ~ismissing(cellContent) && strlength(cellContent) > 0
            newColumn{j} = strtrim(strsplit(cellContent, ','));
        else
            newColumn{j} = {};
        end
    end
    masterTable.(colName) = newColumn;
end

if ismember('DayOfExperience', masterTable.Properties.VariableNames)
    masterTable.DayOfExperience = str2double(masterTable.DayOfExperience);
end
if ismember('Exclude', masterTable.Properties.VariableNames)
    masterTable.Exclude = str2double(masterTable.Exclude);
end
end