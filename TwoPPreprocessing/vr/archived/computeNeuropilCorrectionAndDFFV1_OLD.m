function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFFV1_OLD(sessionFileInfo, stimName, applyTemporalSmoothing, prctl_F, window_size)

if nargin < 3, applyTemporalSmoothing = true; end
if nargin < 4, prctl_F = 8; end % The percentile from which to take F0 (baseline F).
if nargin < 5, window_size = 60; end % The rolling window over which to calculate F0.
    
stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end

%% Load data
disp('Loading processedTwoPData...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData');


%% Pre-processing: Optional Temporal Smoothing (recommended for VRStim)

F = processedTwoPData.F;    % Original F (ROI x time)
Fneu = processedTwoPData.Fneu; % Original Fneu (ROI x time)


% Check for and fix non-finite values (NaNs or Infs) before processing.
disp('Checking for non-finite values in F and Fneu...');
if any(~isfinite(F(:))) || any(~isfinite(Fneu(:)))
    fprintf('Found non-finite values. Fixing by linear interpolation...\n');
    % Use 'fillmissing' to interpolate across gaps for each ROI.
    % The 'linear' method connects the points before and after the gap.
    % We specify 'SamplePoints' to ensure interpolation happens along the time dimension (columns).
    F = fillmissing(F, 'linear', 2, 'EndValues', 'nearest');
    Fneu = fillmissing(Fneu, 'linear', 2, 'EndValues', 'nearest');
    fprintf('Fix complete.\n');
else
    fprintf('No non-finite values found. Continuing.\n');
end


if applyTemporalSmoothing
    disp('Applying temporal smoothing to F and Fneu time-series...');
    w = gausswin(9);
    w = w / sum(w);
    
    numROIs = size(F, 1);
    % Using fSmoothed and fneuSmoothed >>>
    fSmoothed = zeros(size(F));
    fneuSmoothed = zeros(size(Fneu));
    
    % Loop through each ROI to apply the filter along the time dimension
    for i = 1:numROIs
        fSmoothed(i, :) = filtfilt(w, 1, F(i, :));
        fneuSmoothed(i, :) = filtfilt(w, 1, Fneu(i, :));
    end
else
    disp('Skipping temporal smoothing.');
    % Might be good to change variable name here to something more neutral 
    fSmoothed = F;
    fneuSmoothed = Fneu;
end

%% Prepare all four signal matrices
disp('Preparing all four signal types from data: fRaw, fNeuropilCorrected, dFF and dFFNeuropilCorrected');
fs = processedTwoPData.ops{1}.fs;  % sampling rate

% Raw F
% processedTwoPData.processedSignals.f = fSmoothed;
% processedTwoPData.temporalSmoothed = applyTemporalSmoothing;

% Neuropil-Corrected F (Fc)
disp('Computing neuropil correction...');
[Fc, ~, ~, ~] = correct_neuropil(fSmoothed', fneuSmoothed', fs);
% processedTwoPData.processedSignals.fNeuropilCorrected = Fc'; % Transpose back to ROI x time

% Delta F/F on Raw F
disp('Computing delta f/f...');
f0Raw = get_F0(fSmoothed', fs, prctl_F, window_size)';
processedTwoPData.processedSignals.dFF = get_delta_F_over_F(fSmoothed, f0Raw);

% Delta F/F on Neuropil-Corrected F
disp('Computing delta f/f on neuropil corrected f...');
F0_c = get_F0(Fc, fs, prctl_F, window_size)'; 
processedTwoPData.processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', F0_c); % Transpose back to ROI x time
disp('Signal matrices created: fRaw, fNeuropilCorrected, dFF & dFFNeuropilCorrected');


disp('Saving processed data files to processedTwoPData..');
save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData', '-v7.3');

end