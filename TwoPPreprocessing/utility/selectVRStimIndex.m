function selectedIndex = selectVRStimIndex(sessionFileInfo)
% Selects the correct VR stimulus file index from sessionFileInfo.
%
% 1. Finds all stimFiles with 'VRCorr' in their name.
% 2. If 0 found -> Prints error, returns [].
% 3. If 1 found -> Prints message, returns that index.
% 4. If >1 found:
%    a. Looks for a file with 'CombinedRuns'. If found, returns the
%       *first* 'CombinedRuns' index.
%    b. If 'CombinedRuns' is not found, returns the index of the
%       *last* 'VRCorr' file found.

% Start by initializing the output as empty.
% This will be the return value if an error occurs.
selectedIndex = []; 

fprintf('Processing Mouse %s session %s...\n', sessionFileInfo.animal_name, sessionFileInfo.session_name);

% Find all VRStim indices
VRStimIndices = find(contains({sessionFileInfo.stimFiles.name}, 'VRCorr'));

if isempty(VRStimIndices)
    fprintf('Error: No VRCorr stim file found.\n');
    return; % Exit the function, returning selectedIndex = []
    
elseif isscalar(VRStimIndices)
    fprintf('VRCorr file found. Loading (Index: %d).\n', VRStimIndices);
    selectedIndex = VRStimIndices;
    
else % Multiple VRCorr files found, need to decide which one
    
    % Check if a 'CombinedRuns' file exists
    combinedRunIndices = find(contains({sessionFileInfo.stimFiles.name}, 'CombinedRuns'));
    
    if ~isempty(combinedRunIndices)
        % A 'CombinedRuns' file exists. Use the first one.
        selectedIndex = combinedRunIndices(1); 
        fprintf('Multiple VRCorr files found. Selecting CombinedResponseRun (Index: %d).\n', selectedIndex);
        
        if ~isscalar(combinedRunIndices)
             fprintf('Warning: Multiple CombinedRuns files found. Using the first one.\n');
        end
          
    else
        % No 'CombinedRuns' file. Use the last VRCorr file.
        selectedIndex = VRStimIndices(end); 
        fprintf('Multiple VRCorr files found and Response Structs not combined. Choosing the last Response Struct (Index: %d).\n', selectedIndex);
    end
end

% The function automatically ends here, returning the value of 'selectedIndex'
end