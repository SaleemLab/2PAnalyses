function [decodedPos, truePos, decodeError, confusionMat] = bayesianDecoder_LapByLap(RegionData, signalName)
% Bayesian position decoder trained on odd laps, tested on even laps
% 
% Input:
%   RegionData  - struct array with ConditionData.Baseline.LapActivity and FilteredROIs
%   signalName  - string, e.g. 'dFFNeuropilCorrected'
%
% Output:
%   decodedPos  - decoded position for each test bin
%   truePos     - true position for each test bin
%   decodeError - absolute error in cm
%   confusionMat - numBins x numBins confusion matrix

if nargin < 2, signalName = 'dFFNeuropilCorrected'; end

fprintf('Running Bayesian decoder for %s...\n', signalName);

allDecodedPos  = [];
allTruePos     = [];
allDecodeError = [];
confusionMat   = zeros(200, 200);

for s = 1:length(RegionData)
    sess = RegionData(s);
    
    if ~isfield(sess, 'ConditionData') || ~isfield(sess.ConditionData, 'Baseline') || ...
       ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
        continue;
    end
    
    lapActivity = sess.ConditionData.Baseline.LapActivity;
    [numROIsTotal, numLaps, numBins] = size(lapActivity);
    
    % smooth
    w = gausswin(15); w = w / sum(w);
    smoothedActivity = lapActivity;
    for iCell = 1:numROIsTotal
        for iLap = 1:numLaps
            trace = squeeze(lapActivity(iCell, iLap, :));
            if all(isnan(trace)), continue; end
            nanMask = isnan(trace); trace(nanMask) = 0;
            smoothed = filtfilt(w, 1, trace);
            smoothed(nanMask) = NaN;
            smoothedActivity(iCell, iLap, :) = smoothed;
        end
    end
    
    % use filtered ROIs only
    roiActivity = smoothedActivity(sess.FilteredROIs, :, :);  % nCells x nLaps x nBins
    nCells = length(sess.FilteredROIs);
    
    if nCells == 0, continue; end
    
    % split odd/even laps
    oddLaps  = 1:2:numLaps;
    evenLaps = 2:2:numLaps;
    
    % --- TRAINING: compute mean tuning curve per cell on odd laps ---
    % tuningCurves: nCells x nBins
    tuningCurves = squeeze(mean(roiActivity(:, oddLaps, :), 2, 'omitnan'));
    if nCells == 1, tuningCurves = tuningCurves'; end
    
    % --- TESTING: decode position on each even lap ---
    for iLap = 1:length(evenLaps)
        lapIdx = evenLaps(iLap);
        
        % activity on this test lap: nCells x nBins
        testActivity = squeeze(roiActivity(:, lapIdx, :));
        if nCells == 1, testActivity = testActivity'; end
        
        % only decode bins 30-170 (exclude unreliable edges)
        for iBin = 30:170
            % population activity vector at this position bin
            popVec = testActivity(:, iBin);  % nCells x 1
            
            if all(isnan(popVec)), continue; end
            
            % --- Bayesian decoding ---
            % compute log likelihood of each position given population activity
            % using Gaussian likelihood: P(r|x) ∝ exp(-0.5 * sum((r - f(x))^2))
            logLikelihood = NaN(numBins, 1);
            for xBin = 1:numBins
                expected = tuningCurves(:, xBin);  % expected activity at position xBin
                diff = popVec - expected;
                % ignore NaN cells
                validCells = ~isnan(diff);
                if sum(validCells) < 3, continue; end
                logLikelihood(xBin) = -0.5 * sum(diff(validCells).^2);
            end
            
            % flat prior (uniform occupancy assumption)
            % posterior ∝ likelihood (since prior is flat)
            [~, decodedBin] = max(logLikelihood);
            
            % store results
            allDecodedPos  = [allDecodedPos; decodedBin];
            allTruePos     = [allTruePos; iBin];
            allDecodeError = [allDecodeError; abs(decodedBin - iBin)];
            
            % update confusion matrix
            confusionMat(iBin, decodedBin) = confusionMat(iBin, decodedBin) + 1;
        end
    end
    
    fprintf('Session %d: %d cells, %d test laps decoded\n', s, nCells, length(evenLaps));
end

decodedPos  = allDecodedPos;
truePos     = allTruePos;
decodeError = allDecodeError;

fprintf('\n=== Decoding Results ===\n');
fprintf('Mean error:   %.1f cm\n', mean(decodeError, 'omitnan'));
fprintf('Median error: %.1f cm\n', median(decodeError, 'omitnan'));

%% --- Plot ---
figure('Color', 'w', 'Position', [100 100 900 400]);

% confusion matrix
subplot(1, 2, 1);
imagesc(1:200, 1:200, confusionMat);
colormap(gca, flipud(gray));
hold on;
xline(40,  'r--', 'LineWidth', 1.2);
xline(80,  'b--', 'LineWidth', 1.2);
xline(120, 'r--', 'LineWidth', 1.2);
xline(160, 'b--', 'LineWidth', 1.2);
yline(40,  'r--', 'LineWidth', 1.2);
yline(80,  'b--', 'LineWidth', 1.2);
yline(120, 'r--', 'LineWidth', 1.2);
yline(160, 'b--', 'LineWidth', 1.2);
xlabel('Decoded position (cm)');
ylabel('True position (cm)');
title('Confusion matrix');
axis square;
colorbar;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 11);

% decoding error distribution
subplot(1, 2, 2);
histogram(decodeError, 'BinEdges', 0:5:200, ...
    'Normalization', 'probability', ...
    'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xline(40,  'r--', '40cm', 'LineWidth', 1.2);
xline(80,  'b--', '80cm', 'LineWidth', 1.2);
xline(120, 'r--', '120cm', 'LineWidth', 1.2);
xlabel('Absolute decoding error (cm)');
ylabel('Probability');
title(sprintf('Decoding error\nMedian = %.1f cm', median(decodeError, 'omitnan')));
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 11);

sgtitle(sprintf('Bayesian Position Decoding — %s', signalName), 'FontName', 'Arial');

%% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4\BayesianDecoding';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(gcf, fullfile(outputDir, sprintf('BayesianDecoding_%s', signalName)));

end