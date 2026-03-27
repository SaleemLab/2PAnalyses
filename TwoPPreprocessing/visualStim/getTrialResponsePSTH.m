function [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, stimName, signalToUse)

if nargin < 3, signalToUse = 'dFF'; end 

% Locate stimulus
iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
assert(~isempty(iStim), 'Stimulus "%s" not found.', stimName);

% load
mergedPath = sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData;
load(mergedPath, 'processedSignals');
load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
load(sessionFileInfo.stimFiles(iStim).Response, 'response');

signalMatrix = processedSignals.(signalToUse); 
[nNeurons, ~] = size(signalMatrix);
nGroups = numel(bonsaiData.trialGroups);

% define trial window 
% Since responseFrameIdx are masks, we sum them to find the number of "true" frames per trial
allFrameCounts = cellfun(@sum, response(1).responseFrameIdx); 
maxFrames = mode(allFrameCounts); 

fprintf('Processing %s: %d neurons, ~%d frames per trial.\n', stimName, nNeurons, maxFrames);

pd = struct('stimValue', [], 'alignedResponses', [], 'meanResponse', [], ...
            'stdResponse', [], 'semResponse', [], 'timeVector', [], 'responseType', []);
response(1).psthData = repmat(pd, nGroups, 1);

% main processing loop; as long as trial groups are saved for all visual
% stimuli this should work.. 
for g = 1:nGroups
    grp     = bonsaiData.trialGroups(g);
    trIdxs  = grp.trials;
    nTrials = numel(trIdxs);
    
    %  [Neurons x TrialFrames x NumTrials]
    aligned = nan(nNeurons, maxFrames, nTrials);
    
    for ti = 1:nTrials
        trialID = trIdxs(ti);
        fMask = response(1).responseFrameIdx{trialID}; % This is the logical mask 
        
        if isempty(fMask) || ~any(fMask), continue; end
        
        % Use the mask to pull ONLY the trial frames from signalMatrix
        % This results in a matrix of size [nNeurons x nFramesInTrial]
        trialData = signalMatrix(:, fMask); 
        
        % Determine how much of trialData fits into our 'aligned' window
        nFramesToCopy = min(size(trialData, 2), maxFrames);
        
        %  Assign to the 3D matrix
        aligned(:, 1:nFramesToCopy, ti) = trialData(:, 1:nFramesToCopy);
    end
    
    % Compute 
    mResp = nanmean(aligned, 3);
    sResp = nanstd(aligned, 0, 3);
    semR  = sResp ./ sqrt(sum(~isnan(aligned), 3));
    
    % Time vector
    sampleTrial = trIdxs(find(allFrameCounts(trIdxs) >= maxFrames, 1));
    if ~isempty(sampleTrial)
        tVec = response(1).responseFrameRelTimes{sampleTrial};
        % Ensure tVec matches maxFrames length
        if length(tVec) > maxFrames, tVec = tVec(1:maxFrames); end
    else
        tVec = 1:maxFrames;
    end
    
    % store in response 
    response(1).psthData(g).stimValue        = grp.value;
    response(1).psthData(g).alignedResponses = aligned;
    response(1).psthData(g).meanResponse     = mResp;
    response(1).psthData(g).stdResponse      = sResp;
    response(1).psthData(g).semResponse      = semR;
    response(1).psthData(g).timeVector       = tVec;
    response(1).psthData(g).responseType     = 'Raw Frames (Logical Mask Indexing)';
end

% 8. Save
save(sessionFileInfo.stimFiles(iStim).Response, 'response');
fprintf('PSTH processing complete for %s.\n', stimName);
end