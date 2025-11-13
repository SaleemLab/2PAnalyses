function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF_Flatten(sessionFileInfo, stimName, applyTemporalSmoothing, prctlF, windowSize, overwrite)
% Set default arguments if not provided
if nargin < 3, applyTemporalSmoothing = true; end
if nargin < 4, prctlF = 8; end % The percentile from which to take F0 (baseline F).
if nargin < 5, windowSize = 60; end % The rolling window over which to calculate F0.
if nargin < 6, overwrite = true; end % Overwrite and rerun 

stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

filePath = sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData;

%% Load data (Robust method)
disp('Loading processedTwoPData...');
% Load all variables in the file into a temporary structure
loadedData = load(filePath);

% Check if it's the old nested format or the new flattened format
if isfield(loadedData, 'processedTwoPData')
    % Old format: data is inside the 'processedTwoPData' field
    processedTwoPData = loadedData.processedTwoPData;
else
    % New format: the file content itself IS the data structure
    processedTwoPData = loadedData;
end
clear loadedData; % Clean up memory

%% Overwrite 
if overwrite && isfield(processedTwoPData, 'processedSignals')
    disp('Overwrite is true. Removing processedSignals field...');
    fieldsToRemove = {'processedSignals'}; 
    processedTwoPData = rmfield(processedTwoPData, intersect(fieldsToRemove, fieldnames(processedTwoPData)));
end

%% Pre-processing and Smoothing
F = processedTwoPData.F;    % Interpolated F (ROI x time)
Fneu = processedTwoPData.Fneu; % Interpolated Fneu (ROI x time)

if applyTemporalSmoothing
    disp('Applying temporal smoothing by removing and re-inserting NaNs...');
    w = gausswin(9); % 100ms smoothning  
    w = w / sum(w);
    
    numROIs = size(F, 1);
    fSmoothed = nan(size(F));
    fneuSmoothed = nan(size(Fneu));
    
    % Loop through each ROI to apply the filter
    for i = 1:numROIs
        %Process F signal for the current ROI 
        originalFRoi = F(i, :);
        nanLocationsF = isnan(originalFRoi);
        fWithoutNans = originalFRoi(~nanLocationsF);
        
        if length(fWithoutNans) > length(w)
            smoothedFClean = filtfilt(w, 1, fWithoutNans);
            fSmoothed(i, ~nanLocationsF) = smoothedFClean;
        else
            fSmoothed(i, ~nanLocationsF) = fWithoutNans;
        end
        
        % --- Repeat for the Fneu signal ---
        originalFneuRoi = Fneu(i, :);
        nanLocationsFneu = isnan(originalFneuRoi);
        fneuWithoutNans = originalFneuRoi(~nanLocationsFneu);
        
        if length(fneuWithoutNans) > length(w)
            smoothedFneuClean = filtfilt(w, 1, fneuWithoutNans);
            fneuSmoothed(i, ~nanLocationsFneu) = smoothedFneuClean;
        else
            fneuSmoothed(i, ~nanLocationsFneu) = fneuWithoutNans;
        end
    end
else
    disp('Skipping temporal smoothing.');
    fSmoothed = F;
    fneuSmoothed = Fneu;
end

%% Prepare all four signal matrices
disp('Preparing all four signal types...');
fs = processedTwoPData.ops{1}.fs;  % sampling rate

% Neuropil-Corrected F (Fc)
disp('Computing neuropil correction...');
[Fc, ~, ~, ~] = correct_neuropil(fSmoothed', fneuSmoothed', fs);

% Delta F/F on Raw F
disp('Computing delta f/f...');
f0Raw = get_F0(fSmoothed', fs, prctlF, windowSize)';
processedTwoPData.processedSignals.dFF = get_delta_F_over_F(fSmoothed, f0Raw);

% Delta F/F on Neuropil-Corrected F
disp('Computing delta f/f on neuropil corrected f...');
f0c = get_F0(Fc, fs, prctlF, windowSize)'; 
processedTwoPData.processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0c);

%% Save data (Flattened)
disp('Signal matrices created. Saving individual fields...');
% Use '-struct' to save the fields of processedTwoPData as individual variables
save(filePath, '-struct', 'processedTwoPData', '-v7.3');
end