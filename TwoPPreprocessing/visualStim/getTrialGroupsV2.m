function [bonsaiData, sessionFileInfo] = getTrialGroupsV2(sessionFileInfo, stimName)
% Automatically groups RFMapping trials into 4-stim blocks,
% and assigns them to trial groups based on unique (X, Y) combinations.
%
% Output:
%   - trialGroups(i).value = [X Y]
%   - trialGroups(i).trials = trial indices (1–120)
%   - trialGroups(i).stimIndices = stimulus index sets (e.g. [1 2 3 4])

% Locate stimulus
isStim = strcmp(stimName, {sessionFileInfo.stimFiles.name});
iStim = find(isStim, 1);
if isempty(iStim)
    error('Stimulus name "%s" not found in sessionFileInfo.', stimName);
end

% Load Bonsai data
load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
load(sessionFileInfo.stimFiles(iStim).Response, 'response');

nStim = length(bonsaiData.positionX);
nStimPerTrial = 4;
nTrials = floor(nStim / nStimPerTrial);  % usually 480/4 = 120

% Build trial stim index mapping
stimGroups = reshape(1:nTrials*nStimPerTrial, nStimPerTrial, [])';  % [120 x 4]
firstStim = stimGroups(:, 1);  % use first stim per trial

% Extract position for each trial
trialX = bonsaiData.positionX(firstStim);
trialY = bonsaiData.positionY(firstStim);
trialPos = [trialX(:), trialY(:)];

% Group by unique (X, Y)
[uniquePositions, ~, ic] = unique(trialPos, 'rows');

% Build trialGroups struct
trialGroups = struct();
for i = 1:size(uniquePositions, 1)
    trialIdx = find(ic == i);  % trials with this (X, Y)
    trialGroups(i).value = uniquePositions(i, :);  % [X Y]
    trialGroups(i).stimTypeName = sprintf('X%d_Y%d', ...
        uniquePositions(i, 1), uniquePositions(i, 2));
    trialGroups(i).trials = trialIdx;  % trial indices 1–120
    trialGroups(i).stimIndices = num2cell(stimGroups(trialIdx, :), 2);  % e.g. [25 26 27 28]
end

% Save into bonsaiData and file
bonsaiData.trialGroups = trialGroups;
save(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');

end
