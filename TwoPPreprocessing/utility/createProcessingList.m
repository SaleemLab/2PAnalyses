function [mouseInfo, sessionsToProcess] = createProcessingList(varargin)
% Creates a processing list for the pipeline by filtering a master table.
%
% This function now loads the master table from a hardcoded path, applies
% filters, and generates the 'mouseInfo' cell array for the pipeline.
%
% SYNTAX:
%   [mouseInfo, sessionsToProcess] = createProcessingList('FilterName1', value1, ...)
%
% EXAMPLE:
%   [mouseInfo, filteredTable] = createProcessingList('TypeImaged', 'Boutons', ...
%       'HasStimulus', 'VRCorr');

%% Path to mastertable and data loading
masterTablePath = fullfile("Z:\ibn-vision\USERS\Sonali\datatable\MasterExpDatatable.csv");
disp('Loading master table from file...');
masterTable = readtable(masterTablePath);

%% Filter the Master Table
disp("Filtering table to create processing list...");
sessionsToProcess = filterMasterTable(masterTable, varargin{:});

%% Handle the Case of No Matching Sessions
if isempty(sessionsToProcess)
    warning('No sessions matched the specified filters. Returning an empty list.');
    mouseInfo = {}; 
    return;
end

%% Dynamically Generate the mouseInfo List for the Pipeline
disp('Dynamically creating mouseInfo list from filtered sessions...');
uniqueMiceInSelection = unique(sessionsToProcess.MouseID);
mouseInfo = cell(length(uniqueMiceInSelection), 2);

for i = 1:length(uniqueMiceInSelection)
    currentMouseID = uniqueMiceInSelection(i);
    
    % --- THIS IS THE CORRECTED LINE ---
    % BEFORE (Incorrect for cell arrays): isCurrentMouse = (sessionsToProcess.MouseID == currentMouseID);
    % AFTER (Correct): Use strcmp for comparing strings within cell arrays.
    isCurrentMouse = strcmp(sessionsToProcess.MouseID, currentMouseID);
    % --- END OF FIX ---
    
    sessionsForThisMouse = sessionsToProcess.Session(isCurrentMouse);
    
    mouseInfo{i, 1} = char(currentMouseID);
    mouseInfo{i, 2} = cellstr(sessionsForThisMouse);
end

fprintf('Processing list created: %d mice with a total of %d sessions.\n', size(mouseInfo, 1), height(sessionsToProcess));

end