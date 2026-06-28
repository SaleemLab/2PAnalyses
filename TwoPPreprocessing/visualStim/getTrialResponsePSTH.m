% function [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, stimName, signalToUse)
% if nargin < 3, signalToUse = 'dFFNeuropilCorrected'; end 
% % Locate stimulus
% iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
% assert(~isempty(iStim), 'Stimulus "%s" not found.', stimName);
% % load
% mergedPath = sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData;
% load(mergedPath, 'processedSignals');
% load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
% load(sessionFileInfo.stimFiles(iStim).Response, 'response');
% 
% signalMatrix = processedSignals.(signalToUse); 
% [nNeurons, ~] = size(signalMatrix);
% nGroups = numel(bonsaiData.trialGroups);
% % define trial window 
% % Since responseFrameIdx are masks, we sum them to find the number of "true" frames per trial
% allFrameCounts = cellfun(@sum, response(1).responseFrameIdx); 
% maxFrames = mode(allFrameCounts); 
% fprintf('Processing %s: %d neurons, ~%d frames per trial.\n', stimName, nNeurons, maxFrames);
% pd = struct('stimValue', [], 'alignedResponses', [], 'meanResponse', [], ...
%             'stdResponse', [], 'semResponse', [], 'timeVector', [], 'responseType', []);
% response(1).psthData = repmat(pd, nGroups, 1);
% % main processing loop; as long as trial groups are saved for all visual
% % stimuli this should work.. 
% for g = 1:nGroups
%     grp     = bonsaiData.trialGroups(g);
%     trIdxs  = grp.trials;
%     nTrials = numel(trIdxs);
% 
%     %  [Neurons x TrialFrames x NumTrials]
%     aligned = nan(nNeurons, maxFrames, nTrials);
% 
%     for ti = 1:nTrials
%         trialID = trIdxs(ti);
%         fMask = response(1).responseFrameIdx{trialID}; % This is the logical mask 
% 
%         if isempty(fMask) || ~any(fMask), continue; end
% 
%         % Use the mask to pull ONLY the trial frames from signalMatrix
%         % This results in a matrix of size [nNeurons x nFramesInTrial]
%         trialData = signalMatrix(:, fMask); 
% 
%         % Determine how much of trialData fits into our 'aligned' window
%         nFramesToCopy = min(size(trialData, 2), maxFrames);
% 
%         %  Assign to the 3D matrix
%         aligned(:, 1:nFramesToCopy, ti) = trialData(:, 1:nFramesToCopy);
%     end
% 
%     % Compute 
%     mResp = nanmean(aligned, 3);
%     sResp = nanstd(aligned, 0, 3);
%     semR  = sResp ./ sqrt(sum(~isnan(aligned), 3));
% 
%     % Time vector
%     sampleTrial = trIdxs(find(allFrameCounts(trIdxs) >= maxFrames, 1));
%     if ~isempty(sampleTrial)
%         tVec = response(1).responseFrameRelTimes{sampleTrial};
%         % Ensure tVec matches maxFrames length
%         if length(tVec) > maxFrames, tVec = tVec(1:maxFrames); end
%     else
%         tVec = 1:maxFrames;
%     end
% 
%     % store in response 
%     response(1).psthData(g).stimValue        = grp.value;
%     response(1).psthData(g).alignedResponses = aligned;
%     response(1).psthData(g).meanResponse     = mResp;
%     response(1).psthData(g).stdResponse      = sResp;
%     response(1).psthData(g).semResponse      = semR;
%     response(1).psthData(g).timeVector       = tVec;
%     response(1).psthData(g).responseType     = 'Raw Frames (Logical Mask Indexing)';
% end
% % Save
% save(sessionFileInfo.stimFiles(iStim).Response, 'response');
% fprintf('PSTH processing complete for %s.\n', stimName);
% end


function [response, sessionFileInfo] = getTrialResponsePSTH(sessionFileInfo, stimName, signalToUse)
if nargin < 3, signalToUse = 'dFFNeuropilCorrected'; end 


% Load data

iStim = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
if ~isempty(iStim)
    filePath = sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData;
    fileInfo = whos('-file', filePath);
    if any(strcmp({fileInfo.name}, 'spks'))
        processedTwoPData = load(filePath, 'processedSignals', 'spks');
    else
        processedTwoPData = load(filePath, 'processedSignals');
    end 
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).Response, 'response'); 
else
    warning('Stimulus name "%s" not found in sessionFileInfo.', stimName);
end

if contains(signalToUse, 'spks')
    signalMatrix = processedTwoPData.spks;
else
    signalMatrix = processedTwoPData.processedSignals.(signalToUse); 
end 

[nNeurons, ~] = size(signalMatrix);
nGroups = numel(bonsaiData.trialGroups);

% define grid
% Use the window defined in the slicing function
pre = response(1).preStimTime;
post = response(1).postStimTime;
Time = -pre : (1/60) : post; 
nFrames = length(Time);

pd = struct('stimValue', [], 'alignedResponses', [], 'meanResponse', [], ...
            'stdResponse', [], 'semResponse', [], 'timeVector', [], 'responseType', []);
response(1).psthData = repmat(pd, nGroups, 1);

for g = 1:nGroups
    grp = bonsaiData.trialGroups(g);
    trIdxs = grp.trials;
    nTrialsInGrp = numel(trIdxs);

    % Initialize 3D matrix with the Perfect Grid size
    aligned = nan(nNeurons, nFrames, nTrialsInGrp);

    for ti = 1:nTrialsInGrp
        trialID = trIdxs(ti);

        % skip bad trials 
        if response(1).badTrialMask(trialID), continue; end

        fMask = response(1).responseFrameIdx{trialID};
        rawRelTimes = response(1).responseFrameRelTimes{trialID}; % jittered times
        rawSignal = signalMatrix(:, fMask);

        % interpolateandbaseline
        for n = 1:nNeurons
            % Ignore NaNs (bad frames) and interpolate to exactly 0.0
            valid = ~isnan(rawRelTimes);
            if sum(valid) > 10
                % Re-align to the timebase
                trace = interp1(rawRelTimes(valid), rawSignal(n, valid), Time, 'linear', 'extrap');

                % Subtract this trial's own baseline (-2s to 0s)
                baselineVal = nanmean(trace(Time < 0));
                aligned(n, :, ti) = trace - baselineVal;
            end
        end
    end

    % Compute Stats
    mResp = nanmean(aligned, 3);
    sResp = nanstd(aligned, 0, 3);
    semR  = sResp ./ sqrt(sum(~isnan(aligned), 3));

    % Store
    response(1).psthData(g).stimValue        = grp.value;
    response(1).psthData(g).alignedResponses = aligned;
    response(1).psthData(g).meanResponse     = mResp;
    response(1).psthData(g).stdResponse      = sResp;
    response(1).psthData(g).semResponse      = semR;
    response(1).psthData(g).timeVector       = Time;
    response(1).psthData(g).responseType     = 'Interpolated & Baseline-Subtracted';
end

save(sessionFileInfo.stimFiles(iStim).Response, 'response');
fprintf('PSTH complete with sub-frame alignment for %s.\n', stimName);
end