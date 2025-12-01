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
numBins = 140; % 1cm bins for a 140cm track
nLaps = size(response.lapPosition2PFrameIdx, 1); % Get nLaps from the 2D cell array
% Shuffle parameters 
maxShift = 2000; % maximum number of elements by which the signal can be circularly shifted (2000)
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
                % Calculate the mean activity for the shifted signal
                % 1. Extract all activities (across all laps) for this bin:
                activityInBin = shiftedSignalMatrix(:, allFramesInBin);
                
                % 2. Calculate the mean activity across all those frames:
                % (This is the average across all laps for this bin/shuffle)
                meanActivity = mean(activityInBin, 2, 'omitnan');
                
                % Store the result for this bin and shift
                meanLapActivityForShift(:, thisBin) = meanActivity;
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

% function [response, sessionFileInfo] = getLapPositionActivity_withShuffle(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs, useZScoredProcessedSignals)
% % Calculates and saves binned lap activity for four signal types, saving
% % ONLY the mean of the shifted distribution. Optimised by pre-calculating
% % all circularly shifted signal matrices.
% %
% % Aman and Sonali February 2025
% % Modified Oct 2025 - Optimised for Speed
% % Modified Nov 2025 - With shuffle 
% %% Handle optional inputs
% if nargin < 3, overwrite = true; end % Default overwrite to false
% if nargin < 4, onlyIncludeROIs = false; end
% if nargin < 5, useZScoredProcessedSignals = true; end
% %% Find VR stimulus and load data
% stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
% if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
% disp('Loading processedTwoPData and Response structs...');
% load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
% load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');
% %% Overwrite check
% if overwrite && isfield(response, 'lapPositionActivity')
%     disp('Overwrite is true. Removing old analysis fields...');
%     fieldsToRemove = {'lapPositionActivity', 'lapPositionActivity_meanShift', 'cellROIs'};
%     response = rmfield(response, intersect(fieldsToRemove, fieldnames(response)));
% end
% %% Get cell ROIs if needed
% if onlyIncludeROIs
%     ROIs = find(processedTwoPData.iscell(:, 1));
% else
%     ROIs = 1:size(processedTwoPData.iscell, 1);
% end
% numCells = length(ROIs);
% %% Select appropriate signal to use
% if useZScoredProcessedSignals
%     disp('Using zScored dFF and dFFNeuropilCorrected for spatial tuning curves..')
%     signals = processedTwoPData.zScoredProcessedSignals; 
%     response.signalsZScored = true; 
% else  
%     signals = processedTwoPData.processedSignals; 
%     response.signalsZScored = false; 
%     disp('Using dFF and dFFNeuropilCorrected (without zScoring) for spatial tuning curves..')
% end 
% %% Binning parameters 
% if ndims(response.lapPosition2PFrameIdx) == 3
%     response.lapPosition2PFrameIdx = squeeze(response.lapPosition2PFrameIdx);
% end
% numBins = 140; % 1cm bins for a 140cm track
% nLaps = size(response.lapPosition2PFrameIdx, 1); % Get nLaps from the 2D cell array
% % Shuffle parameters 
% maxShift = 2000; % maximum number of elements by which the signal can be circularly shifted (2000)
% numShifts = 1000;  % number of times the randomization loop will run (20)
% % Generates a vector of numShifts random integers.
% rng(1)
% randShifts = randi(maxShift,[1 numShifts]); 
% 
% lapPositionActivity = struct();
% lapPositionActivity_meanShift = struct(); 
% signalNames = fieldnames(signals);
% 
% % --- OPTIMIZATION STEP 1: Pre-calculate all shifted signal matrices ---
% disp('Pre-calculating all circularly shifted signal matrices...');
% shiftedSignals = struct();
% for iSignal = 1:length(signalNames)
%     currentSignalName = signalNames{iSignal};
%     currentSignalMatrix = signals.(currentSignalName);
% 
%     % Dimensions: Neurons x Time x Shifts
% 
%     % Store the shifted matrices in a temporary 3D array
%     [~, nTime] = size(currentSignalMatrix);
%     shifted_3D = nan(numCells, nTime, numShifts); 
% 
%     for thisShift = 1:numShifts
%         thisShiftAmount = randShifts(thisShift);
%         shifted_3D(:, :, thisShift) = circshift(currentSignalMatrix, [0 thisShiftAmount]);
%     end
%     shiftedSignals.(currentSignalName) = shifted_3D;
% 
%     % Initialise storage for the results
%     lapPositionActivity.(currentSignalName) = nan(numCells, nLaps, numBins);
%     lapPositionActivity_meanShift.(currentSignalName) = nan(numCells, nLaps, numBins);
% end
% 
% %% 
% disp('Binning lap-position-activity for all signals, and calculating mean of random shifts...');
% % Loop over laps and bins first
% tl = tic;
% for thisLap = 1:nLaps
%     fprintf('Processing Lap: %d / %d\n', thisLap, nLaps);
%     for thisBin = 1:numBins % This position bin 
%         % Get frame indices for this lap and bin.
%         frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
% 
%         if ~isempty(frameIdx)
%             for iSignal = 1:length(signalNames)
%                 currentSignalName = signalNames{iSignal};
%                 currentSignalMatrix = signals.(currentSignalName);
%                 shifted_3D = shiftedSignals.(currentSignalName); % Load the pre-shifted data
% 
%                 % Real activity 
%                 meanActivityReal = mean(currentSignalMatrix(ROIs, frameIdx), 2, 'omitnan');
%                 lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = meanActivityReal;
% 
%                 % Calculate mean of shifted activity (Vectorized) ---
% 
%                 % Extract all shifted activities at once (Gemini speedy
%                 % version)
%                 % Dimensions: (Neurons x Time_in_Bin x Shifts)
%                 allShiftsActivity = shifted_3D(ROIs, frameIdx, :);
% 
%                 % Calculate mean across the Time dimension (dim 2) for all shifts
%                 % Dimensions: (Neurons x 1 x Shifts)
%                 meanActivity_temp_3D = mean(allShiftsActivity, 2, 'omitnan');
% 
%                 % Calculate the FINAL mean across the Shifts dimension (dim 3)
%                 % Dimensions: (Neurons x 1)
%                 meanActivity_overShifts = mean(meanActivity_temp_3D, 3, 'omitnan');
% 
%                 % Store the resulting 3D vector (Mean over shifts)
%                 lapPositionActivity_meanShift.(currentSignalName)(:, thisLap, thisBin) = squeeze(meanActivity_overShifts);
% 
%             end
%         end
%     end
% end
% toc(tl) 
% 
% % if plotFlag 
% % 
% % 
% % end 
% 
% %% Save results to the response struct
% response.lapPositionActivity = lapPositionActivity;
% % Save the structure containing only the mean activity from the shifts
% response.lapPositionActivity_meanShift = lapPositionActivity_meanShift; 
% 
% if onlyIncludeROIs
%     response.cellROIs = ROIs;
% end
% disp('Saving response with updated lapPositionActivities and mean random shift distributions...');
% save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
% end

% function [response, sessionFileInfo] = getLapPositionActivity_withShuffle(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs, useZScoredProcessedSignals)
% % Calculates and saves binned lap activity for four signal types, including a
% % randomized-shift null distribution. This version saves ONLY the mean of the
% % shifted distribution, not the full 4D dataset.
% %
% % Aman and Sonali February 2025
% % Modified Oct 2025
% %% Handle optional inputs
% if nargin < 3, overwrite = false; end % Default overwrite to false
% if nargin < 4, onlyIncludeROIs = false; end
% if nargin < 5, useZScoredProcessedSignals = true; end
% %% Find VR stimulus and load data
% stimIdx = find(strcmp(VRStimName, {sessionFileInfo.stimFiles.name}));
% if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
% disp('Loading processedTwoPData and Response structs...');
% load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');
% load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');
% %% Overwrite check
% if overwrite && isfield(response, 'lapPositionActivity')
%     disp('Overwrite is true. Removing old analysis fields...');
%     % Fields updated to include the new mean shift variable name
%     fieldsToRemove = {'lapPositionActivity', 'lapPositionActivity_meanShift', 'cellROIs'};
%     response = rmfield(response, intersect(fieldsToRemove, fieldnames(response)));
% end
% %% Get cell ROIs if needed
% if onlyIncludeROIs
%     ROIs = find(processedTwoPData.iscell(:, 1));
% else
%     ROIs = 1:size(processedTwoPData.iscell, 1);
% end
% numCells = length(ROIs);
% %% Select appropriate signal to use
% if useZScoredProcessedSignals
%     disp('Using zScored dFF and dFFNeuropilCorrected for spatial tuning curves..')
%     signals = processedTwoPData.zScoredProcessedSignals; 
%     response.signalsZScored = true; 
% else  
%     signals = processedTwoPData.processedSignals; 
%     response.signalsZScored = false; 
%     disp('Using dFF and dFFNeuropilCorrected (without zScoring) for spatial tuning curves..')
% end 
% %% Binning parameters 
% if ndims(response.lapPosition2PFrameIdx) == 3
%     response.lapPosition2PFrameIdx = squeeze(response.lapPosition2PFrameIdx);
% end
% numBins = 140; % 1cm bins for a 140cm track
% nLaps = size(response.lapPosition2PFrameIdx, 1); % Get nLaps from the 2D cell array
% % Shuffle parameters 
% maxShift = 2000; % maximum number of elements by which the signal can be circularly shifted 
% numShifts = 20;  % number of times the randomization loop will run 
% % Generates a vector of numShifts random integers.
% randShifts = randi(maxShift,[1 numShifts]); 
% 
% lapPositionActivity = struct();
% lapPositionActivity_meanShift = struct(); % New structure for storing the mean across shifts
% 
% % storage for all signal types
% signalNames = fieldnames(signals);
% for thisSignal = 1:length(signalNames)
%     currentSignalName = signalNames{thisSignal};
% 
%     % Storage for the REAL activity (Neurons x Laps x Bins)
%     lapPositionActivity.(currentSignalName) = nan(numCells, nLaps, numBins);
% 
%     % Storage for the MEAN SHIFTED activity (Neurons x Laps x Bins)
%     lapPositionActivity_meanShift.(currentSignalName) = nan(numCells, nLaps, numBins);
% end
% %% 
% disp('Binning lap-position-activity for all signals, and calculating mean of random shifts...');
% % Loop over laps and bins first
% tl = tic;
% for thisLap = 1:nLaps
%     fprintf('Circular shifting for Lap: %d....\n', thisLap)
%     for thisBin = 1:numBins % This position bin 
%         % Get frame indices for this lap and bin.
%         frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
% 
%         % If frames exist in this bin, process all signals and cells
%         if ~isempty(frameIdx)
% 
%             for iSignal = 1:length(signalNames)
%                 currentSignalName = signalNames{iSignal};
%                 currentSignalMatrix = signals.(currentSignalName);
% 
%                 % Initialize a temporary matrix to hold all shifts for the current bin/lap
%                 % Neurons x Shifts. This prevents having to save the full 4D matrix.
%                 currentShiftDistribution = nan(numCells, numShifts); 
% 
%                 % Real activity 
%                 % Vectorised calculation for ALL cells at once (mean across time dim 2).
%                 meanActivityReal = mean(currentSignalMatrix(ROIs, frameIdx), 2, 'omitnan');
% 
%                 % Store the resulting REAL activity vector. 
%                 lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = meanActivityReal;
% 
%                 % fprintf('Circular shifting for Lap: %d....\n', thisLap)
%                 % Shifting activity :SSSSS
%                 for thisShift = 1:numShifts
%                     % Generates randomized shift amount with a maximum of
%                     % 2000
% 
%                     thisShiftAmount = randShifts(thisShift);
% 
%                     % Circularly shift the signal matrix in the time dimension (dim 2)
%                     % [0 thisShiftAmount] 0 row shift (neurons), X column shift (transients across time)
%                     newSignal = circshift(currentSignalMatrix, [0 thisShiftAmount]);
% 
%                     % Calculate mean activity from the shifted signal
%                     meanActivityShuf = mean(newSignal(ROIs, frameIdx), 2, 'omitnan');
% 
%                     % Store the result of the current shift in the temporary matrix
%                     currentShiftDistribution(:, thisShift) = meanActivityShuf;
%                 end
% 
%                 %  MEAN across the shifts
%                 % Takes the mean across the 2nd dimension (the shifts dimension)
%                 meanActivity_overShifts = mean(currentShiftDistribution, 2, 'omitnan');
% 
%                 % Mean over shifts
%                 lapPositionActivity_meanShift.(currentSignalName)(:, thisLap, thisBin) = meanActivity_overShifts;
% 
%             end
%         end
%     end
% end
% toc(tl) 
% %% Save results to the response struct
% response.lapPositionActivity = lapPositionActivity;
% % Save the structure containing only the mean activity from the shifts
% response.lapPositionActivity_Shifted = lapPositionActivity_meanShift; 
% 
% if onlyIncludeROIs
%     response.cellROIs = ROIs;
% end
% 
% % disp('Saving response with updated lapPositionActivities and random shift distributions...');
% % save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
% end