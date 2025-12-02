function [response, sessionFileInfo] = getLapPositionActivity_withShuffle(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs, useZScoredProcessedSignals)
% Calculates and saves binned lap activity for four signal types, saving
% the full Neuron x Position Bin x Shuffle matrix, where each shuffle is the 
% mean across all laps.
%
% Aman and Sonali February 2025
% Modified Oct 2025 - Optimised for Speed
% Modified Nov 2025 - With shuffle - NEW IMPLEMENTATION
%% Handle optional inputs
if nargin < 3, overwrite = true; end % Default overwrite to false
if nargin < 4, onlyIncludeROIs = false; end
if nargin < 5, useZScoredProcessedSignals = true; end
%% Find VR stimulus and load data
stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
disp('Loading processedTwoPData and Response structs...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');
%% Overwrite check
if overwrite && isfield(response, 'lapPositionActivity')
    disp('Overwrite is true. Removing old analysis fields...');
    % Change the field name to store the new shuffle matrix
    fieldsToRemove = {'lapPositionActivity', 'lapPositionActivity_meanShift', 'lapPositionActivity_ShuffleMatrix', 'cellROIs'};
    response = rmfield(response, intersect(fieldsToRemove, fieldnames(response)));
end
%% Get cell ROIs if needed
if onlyIncludeROIs
    ROIs = find(processedTwoPData.iscell(:, 1));
else
    ROIs = 1:size(processedTwoPData.iscell, 1);
end
numCells = length(ROIs);
%% Select appropriate signal to use
if useZScoredProcessedSignals
    disp('Using zScored dFF and dFFNeuropilCorrected for spatial tuning curves..')
    signals = processedTwoPData.zScoredProcessedSignals; 
    response.signalsZScored = true; 
else  
    signals = processedTwoPData.processedSignals; 
    response.signalsZScored = false; 
    disp('Using dFF and dFFNeuropilCorrected (without zScoring) for spatial tuning curves..')
end 
%% Binning parameters 
if ndims(response.lapPosition2PFrameIdx) == 3
    response.lapPosition2PFrameIdx = squeeze(response.lapPosition2PFrameIdx);
end
totalTimes = length(processedTwoPData.TwoPFrameTime);
numBins = 140; % 1cm bins for a 140cm track
nLaps = size(response.lapPosition2PFrameIdx, 1); % Get nLaps from the 2D cell array
% Shuffle parameters 
maxShift = round(totalTimes./2); % maximum number of elements by which the signal can be circularly shifted (2000)
numShifts = 1000;  % number of times the randomization loop will run (20)
% Generates a vector of numShifts random integers.
rng(1)
randShifts = randi(maxShift,[1 numShifts]); 
lapPositionActivity = struct(); % To store the real (unshuffled) activity (Neurons x Laps x Bins)
lapPositionActivity_ShuffleMatrix = struct(); % To store the shuffle results (Neurons x Bins x Shuffles)
signalNames = fieldnames(signals);

% --- NEW OPTIMIZATION STEP: Map Frames to Lap/Bin Indices ---
disp('Mapping all 2P frame indices to their corresponding Lap and Bin...');
% Pre-allocate a vector to hold the [Lap, Bin] index for every 2P frame
totalFrames = size(signals.(signalNames{1}), 2);
% We will use two parallel vectors: one for the Lap index and one for the Bin index
frameToLapMap = nan(1, totalFrames);
frameToBinMap = nan(1, totalFrames);

for thisLap = 1:nLaps
    for thisBin = 1:numBins
        frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
        if ~isempty(frameIdx)
            frameToLapMap(frameIdx) = thisLap;
            frameToBinMap(frameIdx) = thisBin;
        end
    end
end
% --- END NEW OPTIMIZATION STEP ---

%% 
disp('Calculating Real Activity (Neurons x Laps x Bins) and Shuffle Matrix (Neurons x Bins x Shuffles)...');
tl = tic;
for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};
    currentSignalMatrix = signals.(currentSignalName)(ROIs, :); % Select only the ROIs
    
    numROIs = size(currentSignalMatrix, 1);
    
    % Initialize storage for the real activity
    lapPositionActivity.(currentSignalName) = nan(numROIs, nLaps, numBins);
    
    % Initialize storage for the shuffle matrix (Neurons x Bins x Shuffles)
    lapPositionActivity_ShuffleMatrix.(currentSignalName) = nan(numROIs, numBins, numShifts);
    
    % --- Part 1: Calculate Real Activity (Neurons x Laps x Bins) ---
    disp(['Processing Real Activity for: ', currentSignalName]);
    for thisLap = 1:nLaps
        for thisBin = 1:numBins
            frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
            if ~isempty(frameIdx)
                % Mean across frames in this bin/lap (dim 2)
                meanActivityReal = mean(currentSignalMatrix(:, frameIdx), 2, 'omitnan');
                lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = meanActivityReal;
            end
        end
    end
    
    % --- Part 2: Calculate Shuffle Matrix (Neurons x Bins x Shuffles) ---
    disp(['Processing Shuffle Matrix for: ', currentSignalName]);
    for thisShift = 1:numShifts
        if mod(thisShift, 100) == 0 || thisShift == 1
            % Print status every 100 shifts and the first shift
            fprintf('  -> Processing Shuffle %d of %d... (Shift amount: %d)\n', ...
                thisShift, numShifts, randShifts(thisShift));
        end
        thisShiftAmount = randShifts(thisShift);
        
        % Circularly shift the whole signal matrix (Neurons x Time)
        shiftedSignalMatrix = circshift(currentSignalMatrix, [0 thisShiftAmount]);
        
        % Initialize storage for the mean activity across laps for THIS shift
        % Dimensions: (Neurons x Bins)
        meanLapActivityForShift = nan(numROIs, numBins);
        
        % Loop over position bins
        for thisBin = 1:numBins
            % Find all frames that belong to THIS position bin (across ALL laps)
            allFramesInBin = find(frameToBinMap == thisBin);
            
            if ~isempty(allFramesInBin)
                % Calculate the mean activity for the shifted signal %
                % Changed to median 
                % 1. Extract all activities (across all laps) for this bin:
                activityInBin = shiftedSignalMatrix(:, allFramesInBin);
                
                % 2. Calculate the mean activity across all those frames:
                % (This is the average across all laps for this bin/shuffle)
                medianActivityInBin = median(activityInBin, 2, 'omitnan'); % changed to median 
                
                % Store the result for this bin and shift
                meanLapActivityForShift(:, thisBin) = medianActivityInBin;
            end
        end
        
        % Store the full (Neurons x Bins) result into the shuffle matrix (dim 3)
        lapPositionActivity_ShuffleMatrix.(currentSignalName)(:, :, thisShift) = meanLapActivityForShift;
    end
end
toc(tl)

response.lapPositionActivity = lapPositionActivity;
% Save the new shuffle matrix under a new field name
response.lapPositionActivity_ShuffleMatrix = lapPositionActivity_ShuffleMatrix; 
response.cellROIs = ROIs;

disp(['Saving updated Response struct to ', sessionFileInfo.stimFiles(stimIdx).Response]);
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-append');


end

