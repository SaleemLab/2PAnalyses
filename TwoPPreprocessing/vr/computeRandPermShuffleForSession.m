function [response, sessionFileInfo] = computeRandPermShuffleForSession(sessionFileInfo, response, vrStimNamesToStitch, useZScoredProcessedSignals)
%% 1. Setup and Input Handling
if nargin < 4, useZScoredProcessedSignals = true; end
signalVarName = iif(useZScoredProcessedSignals, 'zScoredProcessedSignals', 'processedSignals');

[numROIs, ~, numBins] = size(response.lapPositionActivity.dFF);
numShifts = 1000; 

% Set random seed 
shuffleSeed = 1; 
stream = RandStream('mt19937ar', 'Seed', shuffleSeed);

%% stitch raw signal
fprintf('\n--- Starting RandPerm Shuffle Analysis ---\n');
disp('Stitching RAW signals...');
stitchedRawSignalMatrices = struct();
totalCombinedFrames = 0;

for iStim = 1:length(vrStimNamesToStitch)
    stimName = vrStimNamesToStitch{iStim};
    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}), 1);
    if isempty(stimIdx), error(['File for ', stimName, ' not found.']); end
    
    filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;
    dataChunk = load(filePath, signalVarName);
    signals = dataChunk.(signalVarName);
    signalNames = fieldnames(signals);
    
    for iSignal = 1:length(signalNames)
        sName = signalNames{iSignal};
        if iStim == 1
            stitchedRawSignalMatrices.(sName) = signals.(sName);
        else
            stitchedRawSignalMatrices.(sName) = [stitchedRawSignalMatrices.(sName), signals.(sName)];
        end
    end
    totalCombinedFrames = size(stitchedRawSignalMatrices.(sName), 2);
end

%% create flattened frame map
lapPosIdx = response.lapPosition2PFrameIdx;
frameToBinVector = nan(1, totalCombinedFrames);
for r = 1:size(lapPosIdx, 1)
    for c = 1:size(lapPosIdx, 2)
        frames = lapPosIdx{r, c};
        if ~isempty(frames), frameToBinVector(frames) = c; end
    end
end

validFrameIdx = find(~isnan(frameToBinVector));
validBins = frameToBinVector(validFrameIdx);

%% 
lapPositionActivity_RandPerm = struct();
totalTasks = length(signalNames) * numShifts;
taskCounter = 0;

% Initialize Waitbar
hWait = waitbar(0, 'Initializing Shuffles...', 'Name', 'Permutation Shuffle Progress');
tl = tic;

for iSignal = 1:length(signalNames)
    sName = signalNames{iSignal};
    currentMatrix = stitchedRawSignalMatrices.(sName);
    shuffleCube = nan(numROIs, numBins, numShifts, 'single');
    
    fprintf('Processing Signal: %s (%d ROIs)\n', sName, numROIs);
    
    for iShift = 1:numShifts
        taskCounter = taskCounter + 1;
        
        % Shuffling
        shuffledIndices = randperm(stream, length(validBins));
        shuffledBins = validBins(shuffledIndices);
        
        tempMatrix = nan(numROIs, numBins, 'single');
        for iBin = 1:numBins
            framesInThisBin = validFrameIdx(shuffledBins == iBin);
            if ~isempty(framesInThisBin)
                tempMatrix(:, iBin) = median(currentMatrix(:, framesInThisBin), 2, 'omitnan');
            end
        end
        shuffleCube(:, :, iShift) = tempMatrix;
        
        % Update every 100 shuffles to avoid slowing down the loop
        if mod(iShift, 100) == 0
            prog = taskCounter / totalTasks;
            elapsed = toc(tl);
            estTotal = elapsed / prog;
            remTime = estTotal - elapsed;
            
            waitbar(prog, hWait, sprintf('Signal: %s | Shuffle %d/%d | Rem: %.1fs', ...
                sName, iShift, numShifts, remTime));
            
            fprintf('  -> %s: %d/%d shuffles complete (Time remaining: %.1fs)\n', ...
                sName, iShift, numShifts, remTime);
        end
    end
    lapPositionActivity_RandPerm.(sName) = shuffleCube;
end

close(hWait); % Close waitbar when done
fprintf('--- Shuffle Computation Complete (Total Time: %.2fs) ---\n', toc(tl));

% Add metadata
lapPositionActivity_RandPerm.metadata.seed = shuffleSeed;
lapPositionActivity_RandPerm.metadata.type = 'RandomPermutation';
lapPositionActivity_RandPerm.metadata.dateComputed = char(datetime('now'));

%% 5. Save and Append
response.lapPositionActivity_RandPermMatrix = lapPositionActivity_RandPerm;
stimIdx = find(strcmp(response.stimName, {sessionFileInfo.stimFiles.name}), 1);
savePath = sessionFileInfo.stimFiles(stimIdx).Response;

disp(['Appending results to: ', savePath]);
save(savePath, '-struct', 'response', '-append');

end

function out = iif(condition, true_value, false_value)
    if condition, out = true_value; else, out = false_value; end
end