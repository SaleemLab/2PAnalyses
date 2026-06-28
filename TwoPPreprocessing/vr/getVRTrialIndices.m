function [response, sessionFileInfo] = getVRTrialIndices(sessionFileInfo, VRStimName)
% getVRTrialIndices
% Matches landmark sequences to matrix rows using a simple 1-to-1 mapping.
% Previous version used absolute index 
% Sonali - Feb 2026


stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('VRStimName not found.'); end
response = load(sessionFileInfo.stimFiles(stimIdx).Response);

% Conditions
trialIndices.Baseline = [];
trialIndices.Swap_2_3 = [];
trialIndices.Swap_3_4 = [];
trialIndices.Omit_2   = [];
trialIndices.Omit_3   = [];
trialIndices.Omit_4   = [];

% nLaps is the count of completed sequences available
nLaps = length(response.completedStartTimes);

if contains(VRStimName, 'LandManipCorridor') 
    for thisLap = 1:nLaps
        % The sequence for the current lap
        currentSeq = response.completedLandmarkNames{thisLap};
        if isempty(currentSeq) 
            continue;
        end
        
        % The Index (thisLap) matches the row in response.lapPositionActivity
        % idx = thisLap; 
        
        % Sequence matching 
        if length(currentSeq) == 4 && ...
           strcmp(currentSeq{1}, 'grating_vertical') && ...
           strcmp(currentSeq{2}, 'plaid') && ...
           strcmp(currentSeq{3}, 'grating_vertical') && ...
           strcmp(currentSeq{4}, 'plaid')
            trialIndices.Baseline = [trialIndices.Baseline; thisLap];
            
        elseif length(currentSeq) == 4 && ... % Swap 2-3
               strcmp(currentSeq{1}, 'grating_vertical') && ...
               strcmp(currentSeq{2}, 'grating_vertical') && ...
               strcmp(currentSeq{3}, 'plaid') && strcmp(currentSeq{4}, 'plaid')
            trialIndices.Swap_2_3 = [trialIndices.Swap_2_3; thisLap];
            
        elseif length(currentSeq) == 4 && ... % Swap 3-4
               strcmp(currentSeq{1}, 'grating_vertical') && ...
               strcmp(currentSeq{2}, 'plaid') && ...
               strcmp(currentSeq{3}, 'plaid') && ...
               strcmp(currentSeq{4}, 'grating_vertical')
            trialIndices.Swap_3_4 = [trialIndices.Swap_3_4; thisLap];
            
        elseif length(currentSeq) == 3 && ... % Omit 2
               strcmp(currentSeq{1}, 'grating_vertical') && ...
               strcmp(currentSeq{2}, 'grating_vertical') && ...
               strcmp(currentSeq{3}, 'plaid')
            trialIndices.Omit_2 = [trialIndices.Omit_2; thisLap];
            
        elseif length(currentSeq) == 3 && ... % Omit 3
               strcmp(currentSeq{1}, 'grating_vertical') && ...
               strcmp(currentSeq{2}, 'plaid') && ...
               strcmp(currentSeq{3}, 'plaid')
            trialIndices.Omit_3 = [trialIndices.Omit_3; thisLap];
            
        elseif length(currentSeq) == 3 && ... % Omit 4
               strcmp(currentSeq{1}, 'grating_vertical') && ...
               strcmp(currentSeq{2}, 'plaid') && ...
               strcmp(currentSeq{3}, 'grating_vertical')
            trialIndices.Omit_4 = [trialIndices.Omit_4; thisLap];
        end
    end
else
    % 2025 RSP bouton corridor ('VRCorr' and new BaselineCorridor) are 100% Baseline
    trialIndices.Baseline = (1:nLaps)'; %completed laps
end


response.trialIndicesByCondition = trialIndices;

% Save and update
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

fprintf('Classified %d laps for: %s\n', nLaps, VRStimName);
end