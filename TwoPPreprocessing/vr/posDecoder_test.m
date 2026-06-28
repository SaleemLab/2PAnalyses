% response = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_BaselineCorridor_20260322_00002.mat")
% procTwoP = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_processed2PData_M26003_BaselineCorridor_20260322_00002.mat")
% pairs = struct; 
% pairs.M26003 = ['20260322']

% 
% response = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_BaselineCorridor_20260226_00002.mat")
% procTwoP = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_processed2PData_M25132_BaselineCorridor_20260226_00002.mat")
% pairs = struct; 
% pairs.M25132 = ['20260226']

response = load("Z:\ibn-vision\DATA\SUBJECTS\M25133\Analysis\20260220\M25133_20260220_Response_M25133_BaselineCorridor_20260220_00001.mat")
procTwoP = load("Z:\ibn-vision\DATA\SUBJECTS\M25133\Analysis\20260220\M25133_20260220_processed2PData_M25133_BaselineCorridor_20260220_00001.mat")
pairs = struct; 
pairs.M25133 = ['20260220']



%%
testmouse = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0);

testmousedata = getTuningDataByCondition(testmouse);

testmousedata = appendFilteredROIs(testmousedata,'UseExpVar_SigNullDist', true,'ExpVarSigThreshold', 0.01, 'UseExpVar', true, 'cvExpvarThreshold', 0.1, 'FilterEdgeSMI', true);

FilteredROIs = testmousedata(1).FilteredROIs;
timeVec = procTwoP.TwoPFrameTime;
wheelSpeed = response.wheelSpeed;
mousePos = response.mouseVirtualPosition;
abortedLapIdx = response.abortedLaps_AbsoluteIdx;
lapStartTime = response.completedStartTimes;
lapEndTime = response.completedEndTimes;

% --- Build frame-by-frame lap ID vector ---
nFrames = length(timeVec);
frameLapID = zeros(1, nFrames);

nLaps = length(lapStartTime);
for iLap = 1:nLaps
    lapMask = timeVec >= lapStartTime(iLap) & timeVec <= lapEndTime(iLap);
    frameLapID(lapMask) = iLap;
end

% exclude aborted laps
if ~isempty(abortedLapIdx)
    for iAbort = 1:length(abortedLapIdx)
        frameLapID(frameLapID == abortedLapIdx(iAbort)) = 0;
    end
end

% force everything to column vectors
timeVec     = timeVec(:);
wheelSpeed  = wheelSpeed(:);
mousePos    = mousePos(:);
frameLapID  = frameLapID(:);

% running mask: speed > 1 cm/s AND in a valid lap
runningMask = wheelSpeed > 1 & frameLapID > 0;

fprintf('Total frames: %d\n', nFrames);
fprintf('Running frames in valid laps: %d\n', sum(runningMask));
fprintf('Number of valid laps: %d\n', length(unique(frameLapID(frameLapID > 0))));
fprintf('runningMask size: %d x %d\n', size(runningMask,1), size(runningMask,2));
fprintf('Running frames: %d\n', sum(runningMask));

% --- Extract running frame data ---
framePos_running    = mousePos(runningMask);
frameLapID_running  = frameLapID(runningMask);
frameSpeed_running  = wheelSpeed(runningMask);

% extract fluorescence for filtered ROIs only
fluorTrace = procTwoP.processedSignals.dFFNeuropilCorrected(FilteredROIs, :);
fluorTrace_running = fluorTrace(:, runningMask);

fprintf('fluorTrace size: %d ROIs x %d frames\n', size(fluorTrace_running, 1), size(fluorTrace_running, 2));
fprintf('framePos_running range: %.1f to %.1f cm\n', min(framePos_running), max(framePos_running));
fprintf('Unique laps: %d\n', length(unique(frameLapID_running)));

%% --- Bayesian Position Decoder ---
nROIs      = size(fluorTrace_running, 1);
nFrames    = size(fluorTrace_running, 2);
nBins      = 200;  % 1 cm bins
binEdges   = 0:1:200;
binCentres = (binEdges(1:end-1) + binEdges(2:end)) / 2;  % FIX: added binCentres

% get unique valid laps
validLaps = unique(frameLapID_running);
oddLaps   = validLaps(1:2:end);
evenLaps  = validLaps(2:2:end);

%% --- TRAINING: build tuning curves on odd laps ---
tuningCurves = NaN(nROIs, nBins);
oddMask = ismember(frameLapID_running, oddLaps);

for iBin = 1:nBins
    binMask = oddMask & framePos_running >= binEdges(iBin) & framePos_running < binEdges(iBin+1);
    if sum(binMask) < 3, continue; end
    tuningCurves(:, iBin) = mean(fluorTrace_running(:, binMask), 2, 'omitnan');
end

% smooth tuning curves with 5 cm Gaussian window (Saleem et al.)
w = gausswin(15); w = w/sum(w);  % FIX: gausswin(10) = 5cm with 1cm bins
for iROI = 1:nROIs
    tc = tuningCurves(iROI, :);
    if all(isnan(tc)), continue; end
    nm = isnan(tc); tc(nm) = 0;
    tc = filtfilt(w, 1, tc);
    tc(nm) = NaN;
    tuningCurves(iROI, :) = tc;
end

fprintf('Tuning curves computed for %d ROIs\n', nROIs);

%% --- TESTING: decode position on even laps ---
evenMask   = ismember(frameLapID_running, evenLaps);
testFrames = find(evenMask);

allDecoded = NaN(length(testFrames), 1);
allTrue    = framePos_running(testFrames);

for f = 1:length(testFrames)
    frameIdx = testFrames(f);
    popVec = fluorTrace_running(:, frameIdx);

    if sum(~isnan(popVec)) < 3, continue; end

    % Gaussian log likelihood
    logLik = NaN(nBins, 1);
    for xBin = 1:nBins
        expected = tuningCurves(:, xBin);
        d = popVec - expected;
        valid = ~isnan(d);
        if sum(valid) < 3, continue; end
        logLik(xBin) = -0.5 * sum(d(valid).^2);
    end

    [~, bestBin]   = max(logLik);
    allDecoded(f)  = binCentres(bestBin);  % FIX: store position in cm not bin index
end

%% --- Compute error ---
decodeError = abs(allDecoded - allTrue);

fprintf('Median decoding error: %.1f cm\n', median(decodeError, 'omitnan'));
fprintf('Mean decoding error: %.1f cm\n', mean(decodeError, 'omitnan'));


%% --- Plot ---
binSize = 2;
binEdges_plot = 0:binSize:200;
nBins_plot = length(binEdges_plot) - 1;

% build density map
densityMap = zeros(nBins_plot, nBins_plot);
for f = 1:length(testFrames)
    if isnan(allDecoded(f)), continue; end
    trueBin   = min(nBins_plot, max(1, ceil(allTrue(f) / binSize)));
    decodeBin = min(nBins_plot, max(1, ceil(allDecoded(f) / binSize)));
    densityMap(trueBin, decodeBin) = densityMap(trueBin, decodeBin) + 1;
end

% normalise each row by occupancy
for iBin = 1:nBins_plot
    rowSum = sum(densityMap(iBin, :));
    if rowSum > 0
        densityMap(iBin, :) = densityMap(iBin, :) / rowSum;
    end
end

% normalise to chance and smooth
chanceLevel              = 1 / nBins_plot;
densityMap_chance        = densityMap / chanceLevel;
densityMap_chance_smooth = imgaussfilt(densityMap_chance, 2);

figure('Color', 'w', 'Position', [100 100 500 480]);
imagesc(binEdges_plot(1:end-1), binEdges_plot(1:end-1), log2(densityMap_chance_smooth'));
colormap(redWhiteBlue(-1, 1, 256));

cb = colorbar;
cb.Label.String = 'Prob. density/chance';
caxis(log2([0.25 4.0]));
cb.Ticks = log2([0.5 1.0 2.0]);
cb.TickLabels = {'0.5', '1.0', '2.0'};
cb.TickDirection = 'out';
cb.Box = 'off';
cb.FontName = 'Arial';
cb.FontSize = 10;

hold on;
plot([0 200], [0 200], 'k-', 'LineWidth', 1.5);
xline(40,  'k--', 'LineWidth', 1.2);
xline(80,  'k--', 'LineWidth', 1.2);
xline(120, 'k--', 'LineWidth', 1.2);
xline(160, 'k--', 'LineWidth', 1.2);
yline(40,  'k--', 'LineWidth', 1.2);
yline(80,  'k--', 'LineWidth', 1.2);
yline(120, 'k--', 'LineWidth', 1.2);
yline(160, 'k--', 'LineWidth', 1.2);

xlabel('True position (cm)', 'FontName', 'Arial', 'FontSize', 12);
ylabel('Decoded position (cm)', 'FontName', 'Arial', 'FontSize', 12);

set(gca, 'XTick', [40 80 120 160], ...
    'XTickLabel', {'40', '80', '120', '160'}, ...
    'YTick', [40 80 120 160], ...
    'YTickLabel', {'40', '80', '120', '160'});

title(sprintf('RSP position decoding (n=%d ROIs)\nMedian error = %.1fcm', ...
    nROIs, median(decodeError, 'omitnan')), ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'normal');

axis square;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 11);



%% --- Save ---
% %% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\positionDecoding\BayesianDecoding_m25133_20260220_day2';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(gcf, fullfile(outputDir, 'RSP_BayesianDecoding_DensityMap_redo'));