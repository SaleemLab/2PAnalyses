function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, stimName, zScoreProcessedSignals, applyTemporalSmoothing, applyRedChannelCorrection, prctl_F0, prctl_F, windowSize, smoothW, numN, minNp, maxNp, plotFlag)
% estimates neuropil correction, dff and dff without any neuropil
% correction on 60hz interpolated traces from processedTwoPData 
% Parameters:
%   sessionFileInfo        : [struct]  Contains directory paths; must have .stimFiles subfield
%   stimName               : [string]  Name of the stimulus to process (e.g., 'Gratings')
%   zScoreProcessedSignals : [boolean] If true, standardises final dFF traces (mean=0, std=1)
%   applyTemporalSmoothing : [boolean] If true, applies Gaussian filtering
%                                      to raw F and Neu
%   applyRedChannelCorrection : [boolean] If true, regresses out motion using F_chan2
%   prctl_F0               : [Int]  Percentile for slow-drift baseline (e.g., 8th percentile)
%   prctl_F                : [Int]  Percentile of the measured signal that will be matched to neuropil.
%                                      The default is 5. (e.g., 5th)
%   windowSize             : [Int] Length of moving window in seconds for baseline estimation
%   smoothW                : [Int] Length of Gaussian smoothing window in SAMPLES/BINS
%   numN                   : [Int] Number of bins used to partition the distribution of neuropil values.
%                                  Each bin will be associated with a mean neuropil value and a mean signal value.
%   minNp                  : [Int] Minimum values of neuropil considered, expressed in percentile.
%   maxNp                  : [Int] Maximum values of neuropil considered, expressed in percentile.

%% Set Defaults
if nargin < 3 || isempty(zScoreProcessedSignals), zScoreProcessedSignals = true; end
if nargin < 4 || isempty(applyTemporalSmoothing), applyTemporalSmoothing = true; end
if nargin < 5 || isempty(applyRedChannelCorrection), applyRedChannelCorrection = false; end
if nargin < 6 || isempty(prctl_F0), prctl_F0 = 8; end
if nargin < 7 || isempty(prctl_F), prctl_F = 5; end
if nargin < 8 || isempty(windowSize), windowSize = 60; end % Aman: what are units? Time or bins, good to clarify
if nargin < 9 || isempty(smoothW), smoothW = 5; end
if nargin < 10 || isempty(numN), numN = 20; end
if nargin < 11 || isempty(minNp), minNp = 10; end
if nargin < 12 || isempty(maxNp), maxNp = 90; end
if nargin < 13, plotFlag = false; end

interpolatedFrameRate = 60;
% PMT Offset logic
% Absolute zero value 
% was obtained by averaging the darkest frame over many imaging sessions.
% The  absolute zero value is arbitrary and depends on the voltage range of the PMTs. 
absoluteZeroValue = -23;

stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified stimName not found.'); end

%% Load Data
disp('Loading F, FNeu, ops...');
% Check for Red Channel (F_chan2) in file
varsInFile = whos('-file', sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData);
hasRed = any(strcmp({varsInFile.name}, 'F_chan2'));

if hasRed
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'F_chan2', 'ops');
else
    load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');
    applyRedChannelCorrection = false;
end

%% Absolute zero subtraction & Smoothing
F = F - absoluteZeroValue; % currently not saving the absolute zero subtracted version anywhere; might be good to? (Aman: just need to save absZero)
Fneu = Fneu - absoluteZeroValue;

if applyTemporalSmoothing
    fprintf('Applying smoothing (gausswin %d)...\n', smoothW);
    w = gausswin(smoothW); w = w / sum(w);
    fSmoothed = filtfilt(w, 1, F')';
    fneuSmoothed = filtfilt(w, 1, Fneu')';
    if applyRedChannelCorrection
        % Transpose F_chan2 first to apply filter along frames
        fRedSmoothed = filtfilt(w, 1, (F_chan2 - absoluteZeroValue)')';
    end
else
    fSmoothed = F;
    fneuSmoothed = Fneu;
    if applyRedChannelCorrection, fRedSmoothed = F_chan2 - absoluteZeroValue; end
end

%% Transpose once to use in all the functions below: [ROIs x Frames] to [Frames x ROIs]
fSmoothed = fSmoothed';
fneuSmoothed = fneuSmoothed';
if applyRedChannelCorrection, fRedSmoothed = fRedSmoothed'; end % Ensure [Frames x ROIs]

%% Baseline correction and centering
disp('Calculating (slow-drift) baseline...');
tic;
% get_F0 expects [frames x ROIs]
f0_F = get_F0(fSmoothed, prctl_F0, windowSize, interpolatedFrameRate);

%% Neuropil Correction
disp('Computing neuropil correction...');
% This is the "Corrected Raw Fluorescence" (Fc) required for dF/F
% From Sylvia's depth from 2p @ signal[:, iROI] = iF - (b * iN + a) +
% F0[:, iROI];
[Fc, regPars, ~, ~] = correct_neuropil(fSmoothed, f0_F, fneuSmoothed, interpolatedFrameRate, prctl_F0, prctl_F, windowSize, numN, minNp, maxNp);

%% Delta F/F 
disp('Computing delta F/F signals...');
% Raw dF/F (using smoothed F and the F0 we already computed)
dFF = get_delta_F_over_F(fSmoothed, f0_F);
% Neuropil Corrected dF/F
dFF_NC = get_delta_F_over_F(Fc, f0_F);
toc

%% Red Channel Motion Regression (Additional Step)
if applyRedChannelCorrection
    disp('Regressing out red channel motion artifacts...');
    % Both should now be [Frames x ROIs]
    f0_Red = get_F0(fRedSmoothed, prctl_F0, windowSize, interpolatedFrameRate);
    fRed_f = fRedSmoothed - f0_Red;
    
    medWin = round(10 * interpolatedFrameRate);
    if mod(medWin, 2) == 0, medWin = medWin + 1; end
    fRed_clean = movmedian(fRed_f, medWin, 1);
    
    dFF_final = zeros(size(dFF_NC));
    for thisROI = 1:size(dFF_NC, 2)
        y = dFF_NC(:, thisROI); 
        x = fRed_clean(:, thisROI);
        valid = ~isnan(y) & ~isnan(x);
        
        if sum(valid) > 10
            p = polyfit(x(valid), y(valid), 1);
            % % Subtract ONLY the slope component, 
            dFF_final(:, thisROI) = y - (p(1) * x); %+ p(2)
        else
            dFF_final(:, thisROI) = y;
        end

        % figure;
        % subplot(1,2,1);
        % scatter(x(valid), y(valid), 5, 'filled', 'MarkerFaceAlpha', 0.2);
        % hold on; 
        % plot(x(valid), polyval(p, x(valid)), 'r', 'LineWidth', 2);
        % title('Before: dFF vs Red Channel');
        % xlabel('F chan2 (smoothed and baseline subtracted)'); ylabel('dFF');
        % grid on;
        % 
        % subplot(1,2,2);
        % y_corr = dFF_final(valid, thisROI);
        % scatter(x(valid), y_corr, 5, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerFaceColor', [0 .7 0]);
        % title('After: Corrected dFF vs Red Channel');
        % xlabel('F chan2 (smoothed and baseline subtracted)'); ylabel('Corrected dFF');
        % grid on;
    end
    dFF_NC = dFF_final; 
end

%% Transpose back: convert [Frames x ROIs] to [ROIs x Frames]
processedSignals.dFF = dFF';
processedSignals.dFFNeuropilCorrected = dFF_NC';

% these are for plotting
Fc_plot = Fc';
f0_F_plot = f0_F';
fSmoothed_plot = fSmoothed';
fneuSmoothed_plot = fneuSmoothed';

%% Z-Scoring [for boutons]
if zScoreProcessedSignals
    disp('Z-scoring signals...');
    zScoredProcessedSignals.dFFNeuropilCorrected = zscore(processedSignals.dFFNeuropilCorrected, 0, 2);
    zScoredProcessedSignals.dFF = zscore(processedSignals.dFF, 0, 2);
end

%% Sanity Check Plot
if plotFlag
    figure('Name', 'Neuropil & Motion Correction Sanity Check', 'Color', 'w', 'Position', [100 100 1200 800]);
    t = (0:size(F,2)-1) / interpolatedFrameRate;
    roisToPlot = 1:min(3, size(F,1));
    for thisROI = 1:length(roisToPlot)
        roiIdx = roisToPlot(thisROI);
        subplot(length(roisToPlot), 1, thisROI);
        hold on;
        % Plot using the ROIs x Frames versions
        plot(t, fSmoothed_plot(roiIdx,:), 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw F');
        plot(t, fneuSmoothed_plot(roiIdx,:), 'Color', [0.4 0.6 0.8], 'LineStyle', ':', 'DisplayName', 'Neuropil');
        plot(t, f0_F_plot(roiIdx,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'F0 Baseline');
        plot(t, Fc_plot(roiIdx,:), 'm', 'LineWidth', 0.5, 'DisplayName', 'Neuropil-Corr (Fc)');
       % plot(t, processedSignals.dFFNeuropilCorrected(roiIdx,:), 'g', 'LineWidth', 1, 'DisplayName', 'Final dFF (NC+MC)');
        
        ylabel(['ROI ' num2str(roiIdx)]);
        title(['Neuropil & Motion Correction (Slope r = ' num2str(regPars(2,roiIdx), '%.3f') ')']);
        if thisROI == length(roisToPlot), xlabel('Time (s)'); legend('Location', 'northeastoutside'); end
        grid on; axis tight;
    end
end
%% Saving results
processedTwoPData.processedSignals = processedSignals;
processedTwoPData.zScoredProcessedSignals = zScoredProcessedSignals;
processedTwoPData.neuropCorrPars = regPars;
processedTwoPData.F = F;
processedTwoPData.ops = ops;
processedTwoPData.absoluteZeroValue = absoluteZeroValue; 

save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, ...
    'processedSignals', 'zScoredProcessedSignals', 'absoluteZeroValue', ...
    'regPars', 'F', '-append');
end