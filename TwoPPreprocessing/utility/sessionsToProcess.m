function mouseInfo = sessionsToProcess(filteredTable)
% Get the unique mice from your filtered list
uniqueMice = unique(filteredTable.MouseID);

% Pre-allocate the mouseInfo cell array for efficiency
mouseInfo = cell(length(uniqueMice), 2);

% Loop through each unique mouse to gather their corresponding sessions
for i = 1:length(uniqueMice)
    currentMouseID = uniqueMice(i);

    % Find all rows in the 'filteredTable' that match the current mouse
    isCurrentMouse = (filteredTable.MouseID == currentMouseID);

    % Get all corresponding sessions for this mouse from the filtered list
    sessionsForThisMouse = filteredTable.Session(isCurrentMouse);

    % Populate the mouseInfo array in the format your pipeline needs
    mouseInfo{i, 1} = char(currentMouseID);           % First column is the Mouse ID
    mouseInfo{i, 2} = cellstr(sessionsForThisMouse); % Second column is a cell array of sessions
end

end