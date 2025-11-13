function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, stimName, zScoreProcessedSignals ,applyTemporalSmoothing, prctlF, windowSize, overwrite)
% Set default arguments if not provided
if nargin < 3, zScoreProcessedSignals = true; end 
if nargin < 4, applyTemporalSmoothing = true; end
if nargin < 5, prctlF = 8; end % The percentile from which to take F0 (baseline F).
if nargin < 6, windowSize = 60; end % The rolling window over which to calculate F0.
if nargin < 7, overwrite = true; end % Overwrite and rerun 

stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

%% Load data
disp('Loading processedTwoPData...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');

%% Overwrite 
if overwrite && isfield(processedTwoPData, 'processedSignals')
    disp('Overwrite is true. Removing processedSignals field...');
    fieldsToRemove = {'processedSignals'}; % Add any other related fields
    processedTwoPData = rmfield(processedTwoPData, intersect(fieldsToRemove, fieldnames(processedTwoPData)));
end

%% Pre-processing and Smoothing
F = processedTwoPData.F;    % Interpolated F (ROI x time)
Fneu = processedTwoPData.Fneu; % Interpolated Fneu (ROI x time)

if applyTemporalSmoothing
    disp('Applying temporal smoothing by removing and re-inserting NaNs...');
    w = gausswin(9); % 150ms (16.66 ms per frame)
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

if zScoreProcessedSignals 
    processedTwoPData.zScoredProcessedSignals.dFFNeuropilCorrected = zscore(processedTwoPData.processedSignals.dFFNeuropilCorrected, 0,2); %zscore exactly equal to mean and not 1 STD above the mean 
    processedTwoPData.zScoredProcessedSignals.dFF = zscore(processedTwoPData.processedSignals.dFF, 0,2);
    disp('Zscoring Dff and DffNeuropilCorrected signals..')
end 


disp('Signal matrices created and saved.');
save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData', '-v7.3');
end