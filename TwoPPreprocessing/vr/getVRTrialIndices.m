% function [response, sessionFileInfo] = getVRTrialIndices(sessionFileInfo, VRStimName)
% % getVRTrialIndices: Classifies trials strictly by the sequence of landmark names.
% % Baseline: G-P-G-P  Swap 2-3: G-G-P-P  Swap 3-4: G-P-P-G
% % Omit 2: G-G-P  Omit 3: G-P-P  Omit 4: G-P-G
% 
% stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
% if isempty(stimIdx), error('Specified VRStimName not found.'); end
% disp(['Characterizing trials for: ', VRStimName]);
% 
% response = load(sessionFileInfo.stimFiles(stimIdx).Response);
% 
% trialIndices.Baseline = [];
% trialIndices.Swap_2_3 = [];
% trialIndices.Swap_3_4 = [];
% trialIndices.Omit_2   = [];
% trialIndices.Omit_3   = [];
% trialIndices.Omit_4   = [];
% 
% 
% nLaps = length(response.completedLandmarkNames);
% 
% if contains(VRStimName, 'LandManipCorridor')
%     for iLap = 1:nLaps
%         % Get the sequence of landmark names for this lap
%         currentSeq = response.completedLandmarkNames{iLap};
% 
%         if isempty(currentSeq)
%             continue; 
%         end
% 
%         % BASELINE: Grating, Plaid, Grating, Plaid
%         if length(currentSeq) == 4 && ...
%            strcmp(currentSeq{1}, 'grating_vertical') && strcmp(currentSeq{2}, 'plaid') && ...
%            strcmp(currentSeq{3}, 'grating_vertical') && strcmp(currentSeq{4}, 'plaid')
%             trialIndices.Baseline = [trialIndices.Baseline; iLap];
% 
%         % SWAP 2-3: Grating, Grating, Plaid, Plaid
%         elseif length(currentSeq) == 4 && ...
%                strcmp(currentSeq{1}, 'grating_vertical') && strcmp(currentSeq{2}, 'grating_vertical') && ...
%                strcmp(currentSeq{3}, 'plaid') && strcmp(currentSeq{4}, 'plaid')
%             trialIndices.Swap_2_3 = [trialIndices.Swap_2_3; iLap];
% 
%         % SWAP 3-4: Grating, Plaid, Plaid, Grating
%         elseif length(currentSeq) == 4 && ...
%                strcmp(currentSeq{1}, 'grating_vertical') && strcmp(currentSeq{2}, 'plaid') && ...
%                strcmp(currentSeq{3}, 'plaid') && strcmp(currentSeq{4}, 'grating_vertical')
%             trialIndices.Swap_3_4 = [trialIndices.Swap_3_4; iLap];
% 
%         % OMIT 2: Grating, Grating, Plaid 
%         elseif length(currentSeq) == 3 && ...
%                strcmp(currentSeq{1}, 'grating_vertical') && ...
%                strcmp(currentSeq{2}, 'grating_vertical') && ...
%                strcmp(currentSeq{3}, 'plaid')
%             trialIndices.Omit_2 = [trialIndices.Omit_2; iLap];
% 
%         % OMIT 3: Grating, Plaid, Plaid 
%         elseif length(currentSeq) == 3 && ...
%                strcmp(currentSeq{1}, 'grating_vertical') && ...
%                strcmp(currentSeq{2}, 'plaid') && ...
%                strcmp(currentSeq{3}, 'plaid')
%             trialIndices.Omit_3 = [trialIndices.Omit_3; iLap];
% 
%         % OMIT 4: Grating, Plaid, Grating
%         elseif length(currentSeq) == 3 && ...
%                strcmp(currentSeq{1}, 'grating_vertical') && ...
%                strcmp(currentSeq{2}, 'plaid') && ...
%                strcmp(currentSeq{3}, 'grating_vertical')
%             trialIndices.Omit_4 = [trialIndices.Omit_4; iLap];
%         end
%     end
%     disp('LandManipCorridor: Trials categorised by raw name sequence.');
% else
%     % Standard corridors (VRCorr/Baseline) are 100% Baseline
%     trialIndices.Baseline = (1:nLaps)';
% end
% 
% response.trialIndicesByCondition = trialIndices;
% save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
% save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
% 
% % Log results
% fprintf('Done. Total Laps: %d | Classified: %d\n', nLaps, ...
%     length(trialIndices.Baseline) + length(trialIndices.Swap_2_3) + ...
%     length(trialIndices.Swap_3_4) + length(trialIndices.Omit_2) + ...
%     length(trialIndices.Omit_3) + length(trialIndices.Omit_4));
% end

function [response, sessionFileInfo] = getVRTrialIndices(sessionFileInfo, VRStimName)
% getVRTrialIndices
% Matches landmark sequences to matrix rows using a simple 1-to-1 mapping.
% Previous version used absolute index 

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
               strcmp(currentSeq{1}, 'grating_vertical') && strcmp(currentSeq{2}, 'plaid') && ...
               strcmp(currentSeq{3}, 'plaid') && strcmp(currentSeq{4}, 'grating_vertical')
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
    % Standard corridors (VRCorr/Baseline) are 100% Baseline
    trialIndices.Baseline = (1:nLaps)'; %completed laps
end


response.trialIndicesByCondition = trialIndices;

% Save and update
save(sessionFileInfo.stimFiles(stimIdx).Response, '-struct', 'response', '-append');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

fprintf('Classified %d laps for: %s\n', nLaps, VRStimName);
end