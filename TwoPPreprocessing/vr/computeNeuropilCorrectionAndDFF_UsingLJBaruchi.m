function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, stimName, zScoreProcessedSignals, applyTemporalSmoothing, prctl_F0,prctl_F, windowSize, smoothW)
    % Set default arguments
    % Computes Dff on traces interpolated to 60Hz. 
    if nargin < 3, zScoreProcessedSignals = true; end 
    if nargin < 4, applyTemporalSmoothing = false; end 
    if nargin < 5, prctl_F0 = 5; end 
    if nargin < 6, prctl_F = 5; end 
    if nargin < 7, windowSize = 60; end 
    if nargin <8, smoothW = 15; end 

    % This value represents the absolute zero signal and was obtained by averaging the darkest frame over many imaging sessions. 
    % It is important to note that the absolute zero value is arbitrary and depends on the voltage range of the PMTs.
    % Below is the absolute zero value for our B-scope PMTs; Sylvia's scope
    % absolute zero value is 19520 [and thats probably why she implements this]. 
    absZero = -23; 

    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
    if isempty(stimIdx), error('Specified VRStimName not found.'); end

    %% Load data
    disp('Loading F, FNeu, ops...');
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');

    %% Absolute zero correction (@Aman? From Sylvia's pipeline and after a conversation with her) 
    disp(['Subtracting Absolute Zero: ', num2str(absZero)]);
    F = F - absZero;
    Fneu = Fneu - absZero;

    %% Smoothning (f0 will additionally be smoothed) 
    if applyTemporalSmoothing
        fprintf('Applying smoothing (gausswin %d)...', smoothW);
        w = gausswin(smoothW); 
        w = w / sum(w); 
        % Vectorised smoothing along time dimension
        fSmoothed = filtfilt(w, 1, F')'; 
        fneuSmoothed = filtfilt(w, 1, Fneu')';
    else
        fSmoothed = F; fneuSmoothed = Fneu;
    end

    %% Neuropil Correction (Using sylvia's functions)
    disp('Computing neuropil correction...');
    % handles F0 internally and returns [frames x ROIs]
    % [signalTrace,neuropCorrPars]=estimateNeuropil(cellRoiTrace,neuropRoiTrace,opt);
    [Fc_frames, ~, ~] = correct_neuropil(fSmoothed', fneuSmoothed', prctl_F0, prctl_F, windowSize); 
    Fc = Fc_frames'; % Transpose back to [ROIs x frames]

    %%  Delta F/F calculation
    disp('Computing delta F/F signals...');
    
    % dF/F on Raw F (using fSmoothed)
    f0Raw = get_F0(fSmoothed', prctl_F0, windowSize);
    processedSignals.dFF = get_delta_F_over_F(fSmoothed', f0Raw)';

    % dF/F on Neuropil-Corrected F/ Swapping back to using the F0raw @Aman 
    f0c = get_F0(Fc', prctl_F0, windowSize);
    %processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0c)';
    processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0Raw)';

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

