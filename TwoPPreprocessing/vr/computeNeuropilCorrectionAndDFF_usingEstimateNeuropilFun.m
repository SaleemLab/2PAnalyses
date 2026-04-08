function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF_usingEstimateNeuropilFun(sessionFileInfo, stimName, zScoreProcessedSignals, applyTemporalSmoothing, prctl_F0, prctl_F, windowSize, smoothW, numN, minNp, maxNp)
    %% Set defaults
    if nargin < 3 || isempty(zScoreProcessedSignals), zScoreProcessedSignals = true; end 
    if nargin < 4 || isempty(applyTemporalSmoothing), applyTemporalSmoothing = false; end 
    if nargin < 5 || isempty(prctl_F0), prctl_F0 = 5; end 
    if nargin < 6 || isempty(prctl_F), prctl_F = 5; end 
    if nargin < 7 || isempty(windowSize), windowSize = 60; end 
    if nargin < 8 || isempty(smoothW), smoothW = 15; end 
    if nargin < 9 || isempty(numN), numN = 20; end
    if nargin < 10 || isempty(minNp), minNp = 10; end
    if nargin < 11 || isempty(maxNp), maxNp = 90; end
    
    % PMT Offset logic
    % This value represents the absolute zero signal and was obtained by averaging the darkest frame over many imaging sessions.
    % Aabsolute zero value is arbitrary and depends on the voltage range of the PMTs.
    % Below is the absolute zero value for our B-scope PMTs; Sylvia's scope
    % absolute zero value is 19520 [and thats probably why she implements this].
    absZero = -23; 
    
    stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
    if isempty(stimIdx), error('Specified stimName not found.'); end

    %% Load data
    disp('Loading F, FNeu, ops...');
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');

    %% Absolute zero subtraction and smoothning 
    F = F - absZero;
    Fneu = Fneu - absZero;

    if applyTemporalSmoothing
        fprintf('Applying smoothing (gausswin %d)...\n', smoothW);
        w = gausswin(smoothW); w = w / sum(w); 
        fSmoothed = filtfilt(w, 1, F')';  
        fneuSmoothed = filtfilt(w, 1, Fneu')';
    else
        fSmoothed = F; fneuSmoothed = Fneu;
    end

    %% Baseline calculation & centering
    disp('Calculating slow-drift baselines...');
    % We calculate these once and reuse them for, correct_neuropil dF/F calculation 
    % get_F0 expects [frames x ROIs]
    f0_F = get_F0(fSmoothed', prctl_F0, windowSize)'; 
    f0_N = get_F0(fneuSmoothed', prctl_F0, windowSize)'; 
    
    % Centering signals for estimateNeuropil regression
    fCentered = fSmoothed - f0_F;
    fneuCentered = fneuSmoothed - f0_N;

    %% 5. Neuropil Correction
    disp('Computing neuropil correction using estimateNeuropil...');
    opt.numN = numN; 
    opt.minNp = minNp; 
    opt.maxNp = maxNp;
    opt.pCell = prctl_F; 
    opt.noNeg = 1; 
    opt.constrainedFit = 0; 
    opt.window = Inf;     
    
    % estimateNeuropil output (signalTrace) is centered at 0
    [signalTraceCentered, neuropCorrPars] = estimateNeuropil(fCentered, fneuCentered, opt);
    
    % RESTORE BASELINE: Add f0_F back to the corrected centered signal
    % This is the "Corrected Raw Fluorescence" (Fc) required for dF/F 
    % From Sylvia's depth from 2p @ signal[:, iROI] = iF - (b * iN + a) +
    % F0[:, iROI]; estimate neuropil does not do this.. 
    Fc = signalTraceCentered + f0_F; 

    %% 6. Delta F/F Calculation
    disp('Computing delta F/F signals...');
    
    % Raw dF/F (using smoothed F and the F0 we already computed)
    processedSignals.dFF = get_delta_F_over_F(fSmoothed', f0_F')';
    
    % Neuropil Corrected dF/F 
    % Per Sylvia: Normalize corrected signal (Fc) by the original Raw F0 (f0_F)
    processedSignals.dFFNeuropilCorrected = get_delta_F_over_F(Fc', f0_F')';
 
    %% 7. Z-Scoring
    if zScoreProcessedSignals 
        disp('Z-scoring signals...');
        zScoredProcessedSignals.dFFNeuropilCorrected = zscore(processedSignals.dFFNeuropilCorrected, 0, 2);
        zScoredProcessedSignals.dFF = zscore(processedSignals.dFF, 0, 2);
    end 

    %% 8. Sanity Check Plot
    % figure('Name', 'Neuropil Correction Sanity Check', 'Color', 'w', 'Position', [100 100 1200 800]);
    % t = (0:size(F,2)-1) / 60; % Assuming 60Hz
    % roisToPlot = 1:min(3, size(F,1));
    % 
    % for i = 1:length(roisToPlot)
    %     roiIdx = roisToPlot(i);
    %     subplot(length(roisToPlot), 1, i);
    %     hold on;
    %     plot(t, fSmoothed(roiIdx,:), 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw F');
    %     plot(t, fneuSmoothed(roiIdx,:), 'Color', [0.4 0.6 0.8], 'LineStyle', ':', 'DisplayName', 'Neuropil');
    %     plot(t, f0_F(roiIdx,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'F0 Baseline');
    %     plot(t, Fc(roiIdx,:), 'r', 'LineWidth', 1, 'DisplayName', 'Corrected Fc');
    % 
    %     ylabel(['ROI ' num2str(roiIdx) ' (a.u.)']);
    %     if i == 1, title(['Neuropil Correction (Slope r = ' num2str(neuropCorrPars(roiIdx,2), '%.3f') ')']); end
    %     if i == length(roisToPlot), xlabel('Time (s)'); legend('Location', 'northeastoutside'); end
    %     grid on; axis tight;
    % end

    %% 9. Saving Results
    disp('Saving processed data...');
    save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, ...
        'processedSignals', 'zScoredProcessedSignals', 'neuropCorrPars', '-append');
    
    processedTwoPData.processedSignals = processedSignals;
    processedTwoPData.zScoredProcessedSignals = zScoredProcessedSignals;
    processedTwoPData.neuropCorrPars = neuropCorrPars;
    processedTwoPData.F = F; 
    processedTwoPData.ops = ops;
end