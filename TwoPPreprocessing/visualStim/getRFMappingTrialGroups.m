function [bonsaiData, sessionFileInfo] = getRFMappingTrialGroups(sessionFileInfo, stimName)
% Updated version: Matches stimulus parameters to individual triggers
% Works with the new Response structure where 1 Trigger = 1 Trial.

isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name});
iStim = find(isStim, 1);
if isempty(iStim), error('Stimulus "%s" not found.', stimName); end

load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
load(sessionFileInfo.stimFiles(iStim).Response, 'response');

% We expect 160 for RFMapping (10 repeats); This gets trial counts 
nStimParams = length(bonsaiData.positionX);
nTrials = length(response.responseFrameIdx);

% Sync counts to handle early session stops
actualCount = min(nStimParams, nTrials);

% Extract trial-by-trial coordinates
%Use (:) to force column vectors so [X, Y] is N-by-2
trialX = bonsaiData.positionX(1:actualCount);
trialY = bonsaiData.positionY(1:actualCount);
trialPos = [trialX(:), trialY(:)]; 

% Group by unique (X, Y) positions
[uniquePositions, ~, ic] = unique(trialPos, 'rows');

% Build trialGroups struct
bonsaiData.trialGroups = []; 
trialGroups = struct('value', {}, 'stimTypeName', {}, 'trials', {}, 'stimIndices', {});

for thisUniquePosition = 1:size(uniquePositions, 1)
    % Find every trigger index associated with this X,Y coordinate
    tIdx = find(ic == thisUniquePosition);  
    
    trialGroups(thisUniquePosition).value = uniquePositions(thisUniquePosition, :);  % [X Y]
    trialGroups(thisUniquePosition).stimTypeName = sprintf('X%d_Y%d', uniquePositions(thisUniquePosition, 1), uniquePositions(thisUniquePosition, 2));
    trialGroups(thisUniquePosition).trials = tIdx; 
    
    % trial index and stim index are identical
    trialGroups(thisUniquePosition).stimIndices = tIdx;
end

if nStimParams ~= nTrials
    warning('Count mismatch! Bonsai has %d trials, but Photodiode found %d triggers. Check for dropped triggers!', ...
        nStimParams, nTrials);
end


% Save
bonsaiData.trialGroups = trialGroups;
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');

fprintf('Success: Grouped %d trials into %d unique (X,Y) positions for %s.\n', ...
    actualCount, size(uniquePositions, 1), stimName);
end