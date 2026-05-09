function filteredTable = filterMasterTable(varargin)
% Filters the master data table based on specified criteria.
% It loads and parses the table using a nested function and applies
% flexible, name-value pair filters.

    % Define the path and load the parsed master table using the nested function
    filepath = "Z:\ibn-vision\USERS\Sonali\datatable\MasterExpDatatable.csv";
    masterTable = loadAndParseMasterTable(filepath);
    

    p = inputParser;
    p.KeepUnmatched = true;
    
    % Define default values
    defaultMouseID = string.empty;
    defaultSession = string.empty;
    defaultDay = [];
    defaultHasStimulus =[];
    defaultExclude = []; %
    
    % Add the parameters
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
    
    %% --- Apply filters based on criteria ---
    % Filter by MouseID
    if ~isempty(criteria.MouseID)
        % --- THIS IS THE FIX ---
        % Ensure both sides of the comparison are string arrays for robustness.
        rowIndex = rowIndex & ismember(string(masterTable.MouseID), string(criteria.MouseID));
    end
    
    % Filter by Session
    if ~isempty(criteria.Session)
        rowIndex = rowIndex & ismember(string(masterTable.Session), string(criteria.Session));
    end
    
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
    % Loads the MasterExpDatatable.csv and parses relevant columns.
    
    % Import all as string type to start
    opts = detectImportOptions(filepath, 'VariableNamingRule', 'preserve');
    opts = setvartype(opts, 'string'); 
    masterTable = readtable(filepath, opts);
    
    % Identify columns that need to be converted to cell arrays of strings
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
    
    % Convert specific known columns to numeric types for easier filtering
    if ismember('DayOfExperience', masterTable.Properties.VariableNames)
        masterTable.DayOfExperience = str2double(masterTable.DayOfExperience);
    end
    if ismember('Exclude', masterTable.Properties.VariableNames)
        masterTable.Exclude = str2double(masterTable.Exclude);
    end
end
