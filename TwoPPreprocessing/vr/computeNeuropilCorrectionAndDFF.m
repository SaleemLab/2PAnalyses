function [processedTwoPData, sessionFileInfo] = computeNeuropilCorrectionAndDFF(sessionFileInfo, stimName, zScoreProcessedSignals, applyTemporalSmoothing, prctl_F0, prctl_F, windowSize, smoothW, numN, minNp, maxNp)
% estimates neuropil correction, dff and dff without any neuropil
% correction on 60hz interpolated traces from processedTwoPData 
% 
%% Set Defaults
if nargin < 3 || isempty(zScoreProcessedSignals), zScoreProcessedSignals = true; end
if nargin < 4 || isempty(applyTemporalSmoothing), applyTemporalSmoothing = true; end
if nargin < 5 || isempty(prctl_F0), prctl_F0 = 8; end
if nargin < 6 || isempty(prctl_F), prctl_F = 5; end
if nargin < 7 || isempty(windowSize), windowSize = 60; end
if nargin < 8 || isempty(smoothW), smoothW = 5; end
if nargin < 9 || isempty(numN), numN = 20; end
if nargin < 10 || isempty(minNp), minNp = 10; end
if nargin < 11 || isempty(maxNp), maxNp = 90; end

interpolatedFrameRate = 60;
% PMT Offset logic
absZero = -23;

stimIdx = find(strcmp(stimName, {sessionFileInfo.stimFiles.name}));
if isempty(stimIdx), error('Specified stimName not found.'); end

%% Load Data
disp('Loading F, FNeu, ops...');
load(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, 'F', 'Fneu', 'ops');

%% Absolute zero subtraction & Smoothing
F = F - absZero; % currently not saving the absolute zero subtracted version anywhere; might be good to? 
Fneu = Fneu - absZero;

if applyTemporalSmoothing
    fprintf('Applying smoothing (gausswin %d)...\n', smoothW);
    w = gausswin(smoothW); w = w / sum(w);
    fSmoothed = filtfilt(w, 1, F')';
    fneuSmoothed = filtfilt(w, 1, Fneu')';
else
    fSmoothed = F;
    fneuSmoothed = Fneu;
end

%% Transpose once to use in all the functions below: [ROIs x Frames] to [Frames x ROIs]
fSmoothed = fSmoothed';
fneuSmoothed = fneuSmoothed';

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
% This function removes slow-drifts [F0,N0] from traces and aligns both
% traces at a zero-baseline. It might be important for neuropil
% correction to ensure the correction factor (r) is estimated based on
% common high-frequency fluctuations rather than absolute offsets.
% F0_F is added back to the corrected trace [Fc]; Performs a first-degree polynomial fit
[Fc, regPars, ~, ~] = correct_neuropil(fSmoothed, f0_F, fneuSmoothed, interpolatedFrameRate, prctl_F0, prctl_F, windowSize, numN, minNp, maxNp);

%% Delta F/F 
disp('Computing delta F/F signals...');

% Raw dF/F (using smoothed F and the F0 we already computed)
% Inputs are [Frames x ROIs]
dFF = get_delta_F_over_F(fSmoothed, f0_F);

% Neuropil Corrected dF/F
% Per Sylvia: Normalise corrected signal (Fc) by the original Raw F0 (f0_F)
dFF_NC = get_delta_F_over_F(Fc, f0_F);
toc

%% Transpose back: convert [Frames x ROIs] to [ROIs x Frames]
processedSignals.dFF = dFF';
processedSignals.dFFNeuropilCorrected = dFF_NC';
% these are for plotting
Fc = Fc';
f0_F = f0_F';
fSmoothed=fSmoothed';
fneuSmoothed=fneuSmoothed';

%% Z-Scoring [for boutons]
if zScoreProcessedSignals
    disp('Z-scoring signals...');
    % Standard z-score along the time dimension (dim 2 for [ROIs x Time])
    zScoredProcessedSignals.dFFNeuropilCorrected = zscore(processedSignals.dFFNeuropilCorrected, 0, 2);
    zScoredProcessedSignals.dFF = zscore(processedSignals.dFF, 0, 2);
end

%% Sanity Check Plot
figure('Name', 'Neuropil Correction Sanity Check', 'Color', 'w', 'Position', [100 100 1200 800]);
t = (0:size(F,2)-1) / interpolatedFrameRate;
roisToPlot = 1:min(3, size(F,1));
for i = 1:length(roisToPlot)
    roiIdx = roisToPlot(i);
    subplot(length(roisToPlot), 1, i);
    hold on;
    plot(t, fSmoothed(roiIdx,:), 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw F');
    plot(t, fneuSmoothed(roiIdx,:), 'Color', [0.4 0.6 0.8], 'LineStyle', ':', 'DisplayName', 'Neuropil');
    plot(t, f0_F(roiIdx,:), 'b--', 'LineWidth', 1.2, 'DisplayName', 'F0 Baseline');
    plot(t, Fc(roiIdx,:), 'r', 'LineWidth', 1, 'DisplayName', 'Corrected Fc');
    ylabel(['ROI ' num2str(roiIdx) ' (a.u.)']);
    % regPars is [2 x nROIs], row 2 is the slope (r)
    title(['Neuropil Correction (Slope r = ' num2str(regPars(2,roiIdx), '%.3f') ')']);
    if i == length(roisToPlot), xlabel('Time (s)'); legend('Location', 'northeastoutside'); end
    grid on; axis tight;
end

%% Saving results
disp('Saving processed data...');
save(sessionFileInfo.stimFiles(stimIdx).processedMergedBonsaiSuite2pData, ...
    'processedSignals', 'zScoredProcessedSignals', '-append');

processedTwoPData.processedSignals = processedSignals;
processedTwoPData.zScoredProcessedSignals = zScoredProcessedSignals;
processedTwoPData.neuropCorrPars = regPars; % new
processedTwoPData.F = F;
processedTwoPData.ops = ops;
end