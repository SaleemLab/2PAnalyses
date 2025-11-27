function [response, sessionFileInfo] = getLapPositionActivity_withShuffle(sessionFileInfo, VRStimName, overwrite, onlyIncludeROIs, useZScoredProcessedSignals)
% Calculates and saves binned lap activity for four signal types, saving
% ONLY the mean of the shifted distribution. Optimized by pre-calculating
% all circularly shifted signal matrices.
%
% Aman and Sonali February 2025
% Modified Oct 2025 - Optimized for Speed
%% Handle optional inputs
if nargin < 3, overwrite = false; end % Default overwrite to false
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
    fieldsToRemove = {'lapPositionActivity', 'lapPositionActivity_meanShift', 'cellROIs'};
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
maxShift = 2000; % maximum number of elements by which the signal can be circularly shifted 
numShifts = 20;  % number of times the randomization loop will run 
% Generates a vector of numShifts random integers.
randShifts = randi(maxShift,[1 numShifts]); 

lapPositionActivity = struct();
lapPositionActivity_meanShift = struct(); 
signalNames = fieldnames(signals);

% --- OPTIMIZATION STEP 1: Pre-calculate all shifted signal matrices ---
disp('Pre-calculating all circularly shifted signal matrices...');
shiftedSignals = struct();
for iSignal = 1:length(signalNames)
    currentSignalName = signalNames{iSignal};
    currentSignalMatrix = signals.(currentSignalName);
    
    % Dimensions: Neurons x Time x Shifts
    
    % Store the shifted matrices in a temporary 3D array
    [~, nTime] = size(currentSignalMatrix);
    shifted_3D = nan(numCells, nTime, numShifts); 
    
    for thisShift = 1:numShifts
        thisShiftAmount = randShifts(thisShift);
        shifted_3D(:, :, thisShift) = circshift(currentSignalMatrix, [0 thisShiftAmount]);
    end
    shiftedSignals.(currentSignalName) = shifted_3D;
    
    % Initialise storage for the results
    lapPositionActivity.(currentSignalName) = nan(numCells, nLaps, numBins);
    lapPositionActivity_meanShift.(currentSignalName) = nan(numCells, nLaps, numBins);
end

%% 
disp('Binning lap-position-activity for all signals, and calculating mean of random shifts...');
% Loop over laps and bins first
tl = tic;
for thisLap = 1:nLaps
    fprintf('Processing Lap: %d / %d\n', thisLap, nLaps);
    for thisBin = 1:numBins % This position bin 
        % Get frame indices for this lap and bin.
        frameIdx = response.lapPosition2PFrameIdx{thisLap, thisBin};
        
        if ~isempty(frameIdx)
            for iSignal = 1:length(signalNames)
                currentSignalName = signalNames{iSignal};
                currentSignalMatrix = signals.(currentSignalName);
                shifted_3D = shiftedSignals.(currentSignalName); % Load the pre-shifted data
                
                % Real activity 
                meanActivityReal = mean(currentSignalMatrix(ROIs, frameIdx), 2, 'omitnan');
                lapPositionActivity.(currentSignalName)(:, thisLap, thisBin) = meanActivityReal;
                
                % Calculate mean of shifted activity (Vectorized) ---
                
                % Extract all shifted activities at once (Gemini speedy
                % version)
                % Dimensions: (Neurons x Time_in_Bin x Shifts)
                allShiftsActivity = shifted_3D(ROIs, frameIdx, :);
                
                % Calculate mean across the Time dimension (dim 2) for all shifts
                % Dimensions: (Neurons x 1 x Shifts)
                meanActivity_temp_3D = mean(allShiftsActivity, 2, 'omitnan');
                
                % Calculate the FINAL mean across the Shifts dimension (dim 3)
                % Dimensions: (Neurons x 1)
                meanActivity_overShifts = mean(meanActivity_temp_3D, 3, 'omitnan');
                
                % Store the resulting 3D vector (Mean over shifts)
                lapPositionActivity_meanShift.(currentSignalName)(:, thisLap, thisBin) = squeeze(meanActivity_overShifts);

            end
        end
    end
end
toc(tl) 


%% Save results to the response struct
response.lapPositionActivity = lapPositionActivity;
% Save the structure containing only the mean activity from the shifts
response.lapPositionActivity_meanShift = lapPositionActivity_meanShift; 

if onlyIncludeROIs
    response.cellROIs = ROIs;
end
disp('Saving response with updated lapPositionActivities and mean random shift distributions...');
save(sessionFileInfo.stimFiles(stimIdx).Response, 'response', '-v7.3');
end

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