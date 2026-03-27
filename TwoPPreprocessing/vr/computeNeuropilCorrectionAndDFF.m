function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, stimName, zScoreProcessedSignals, applyTemporalSmoothing, prctlF, windowSize)
    % Set default arguments
    if nargin < 3, zScoreProcessedSignals = true; end 
    if nargin < 4, applyTemporalSmoothing = true; end
    if nargin < 5, prctlF = 8; end 
    if nargin < 6, windowSize = 60; end 

    % Hardcoded Absolute Zero (from your ImageJ vasculature measurement)
    absZero = 23; % measured dark spots in many sessions; primarily over large vasculatures 

    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
    if isempty(stimIdx), error('Specified VRStimName not found.'); end

    %% Load data
    disp('Loading F, FNeu, ops...');
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');

    %% Absolute zero correction (this is new @Aman?) 
    %disp(['Subtracting Absolute Zero Offset: ', num2str(absZero)]);
    %F = F - absZero;
    %Fneu = Fneu - absZero;

    %% Smoothning (f0 will additionally be smoothed) 
    if applyTemporalSmoothing
        disp('Applying temporal smoothing (gausswin 15)...');
        w = gausswin(15); w = w / sum(w);
        % Vectorised smoothing along time dimension
        fSmoothed = filtfilt(w, 1, F')'; 
        fneuSmoothed = filtfilt(w, 1, Fneu')';
    else
        fSmoothed = F; fneuSmoothed = Fneu;
    end

    %% Neuropil Correction (Using sylvia's functions)
    disp('Computing neuropil correction...');
    % handles F0 internally and returns [frames x ROIs]
    [Fc_frames, ~, ~] = correct_neuropil(fSmoothed', fneuSmoothed', prctlF, windowSize);
    Fc = Fc_frames'; % Transpose back to [ROIs x frames]

    %%  Delta F/F calculation
    disp('Computing delta F/F signals...');
    
    % dF/F on Raw F (using fSmoothed)
    f0Raw = get_F0(fSmoothed', prctlF, windowSize);
    processedSignals.dFF = get_delta_F_over_F(fSmoothed', f0Raw)';

    % dF/F on Neuropil-Corrected F
    f0c = get_F0(Fc', prctlF, windowSize);
    processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0c)';

    %% Z-Scoring
    if zScoreProcessedSignals 
        disp('Z-scoring signals...');
        zScoredProcessedSignals.dFFNeuropilCorrected = zscore(processedSignals.dFFNeuropilCorrected, 0, 2);
        zScoredProcessedSignals.dFF = zscore(processedSignals.dFF, 0, 2);
    end 

    %% Saving
    disp('Saving processed data...');
    save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, ...
        'processedSignals', 'zScoredProcessedSignals', '-append');

    processedTwoPData.processedSignals = processedSignals;
    processedTwoPData.zScoredProcessedSignals = zScoredProcessedSignals;
    processedTwoPData.F = F; 
    processedTwoPData.ops = ops;
end

%% Older version 
% function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, stimName, zScoreProcessedSignals ,applyTemporalSmoothing, prctlF, windowSize)
% % Set default arguments if not provided
% if nargin < 3, zScoreProcessedSignals = true; end 
% if nargin < 4, applyTemporalSmoothing = true; end
% if nargin < 5, prctlF = 8; end % The percentile from which to take F0 (baseline F).
% if nargin < 6, windowSize = 60; end % The rolling window over which to calculate F0.
% 
% 
% stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
% if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
% 
% %% Load data
% disp('Loading F, FNeu, ops...');
% load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');
% 
% 
% %% Pre-processing and Smoothing
% if applyTemporalSmoothing
%     disp('Applying temporal smoothing to F and Fneu time-series...');
%     w = gausswin(9); % 144ms smoothning / SD is 25.6ms
%     w = w / sum(w);
% 
%     numROIs = size(F, 1);
%     % Using fSmoothed and fneuSmoothed >>>
%     fSmoothed = zeros(size(F));
%     fneuSmoothed = zeros(size(Fneu));
% 
%     % Loop through each ROI to apply the filter along the time dimension
%     for i = 1:numROIs
%         fSmoothed(i, :) = filtfilt(w, 1, F(i, :));
%         fneuSmoothed(i, :) = filtfilt(w, 1, Fneu(i, :));
%     end
% else
%     disp('Skipping temporal smoothing.');
%     % Might be good to change variable name here to something more neutral 
%     fSmoothed = F;
%     fneuSmoothed = Fneu;
% end
% 
% %% Prepare all four signal matrices
% disp('Preparing all four signal types...');
% fs = ops{1}.fs;  % frame rate
% 
% % Neuropil-Corrected F (Fc)
% disp('Computing neuropil correction...');
% [Fc, ~, ~, ~] = correct_neuropil(fSmoothed', fneuSmoothed', fs);
% 
% % Delta F/F on Raw F
% disp('Computing delta f/f...');
% f0Raw = get_F0(fSmoothed', fs, prctlF, windowSize)';
% processedSignals.dFF = get_delta_F_over_F(fSmoothed, f0Raw);
% 
% % Delta F/F on Neuropil-Corrected F
% disp('Computing delta f/f on neuropil corrected f...');
% f0c = get_F0(Fc, fs, prctlF, windowSize)'; 
% processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0c);
% 
% if zScoreProcessedSignals 
%     zScoredProcessedSignals.dFFNeuropilCorrected = zscore(processedSignals.dFFNeuropilCorrected, 0,2); %zscore exactly equal to mean and not 1 STD above the mean 
%     zScoredProcessedSignals.dFF = zscore(processedSignals.dFF, 0,2);
%     disp('Zscoring Dff and DffNeuropilCorrected signals..')
% end 
% 
% %% Saving (Maintaining the struct as a root variable)
% disp('Appending processedSignals to the file...');
% save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedSignals', 'zScoredProcessedSignals', '-append');
% 
% % Also return for function output (temp) 
% processedTwoPData.processedSignals = processedSignals;
% processedTwoPData.zScoredProcessedSignals = zScoredProcessedSignals;
% processedTwoPData.F = F; 
% processedTwoPData.ops = ops;
% 
% 
% % save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'processedTwoPData', '-v7.3');
% end