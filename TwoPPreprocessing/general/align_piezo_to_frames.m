%% align_piezo_to_frames.m
% Aligns piezo (z-actuator) data to imaging frame times using sync pulses,
% then plots the average piezo trajectory (in microns) across one full
% volume cycle, with each plane's acquisition window overlaid.
%
% This plot is used to VISUALLY IDENTIFY which planes are fly-back planes,
% based on where the piezo trace shows a fast, non-linear return sweep -
% it does not assume the fly-back planes in advance.
%
% Requires: align2PSyncPulses.m (SGS)

%%
%piezoFile = "\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\PiezoReadTest\PRT-piezo\EyeTracking\20260203\PRT_PRT8P2S_20260203_00001_Piezo2026-02-03T18_46_45.csv";
%frameFile = "\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\PiezoReadTest\PRT-frametimes\Bonsai\20260203\PRT_PRT8P2S_20260203_00001_2P2026-02-03T18_47_10.csv";
% 




%% make_flyback_figure.m
% Generates the final figure showing the average piezo trajectory per
% plane, with fly-back planes labelled based on trace direction: planes
% 0, 1, and 2 show a descending trajectory (consistent with the fly-back
% return sweep, including its slower tail near the turnaround point),
% while planes 3-7 show the ascending ramp.

%%
piezoFile = "Z:\ibn-vision\DATA\SUBJECTS\M25133\EyeTracking\20260219\M25133_Darkness_20260219_00001_Piezo2026-02-19T17_43_57.csv";
frameFile = "Z:\ibn-vision\DATA\SUBJECTS\M25133\Bonsai\20260219\M25133_Darkness_20260219_00001_2P2026-02-19T17_44_16.csv";

%% make_flyback_figure.m
% Generates the final figure showing the average piezo trajectory per
% plane, with fly-back planes labelled based on trace direction: planes
% 0, 1, and 2 show a descending trajectory (consistent with the fly-back
% return sweep, including its slower tail near the turnaround point),
% while planes 3-7 show the ascending ramp.


numPlanes = 8;
piezoADCRange = 1023;
piezoMicronRange = 400;
steadyStateUpperBound = 100;   % microns, for trimming a startup transient if present
bufferMs = 200;

% Final result: planes 0, 1, 2 identified as fly-back (confirmed
% directly from the acquisition pipeline's own ignore_flyback setting).
flybackPlanes = [0, 1, 2];

%% oad data 
piezoData = readtable(piezoFile);
frameData = readtable(frameFile);

piezoData.Microns = (double(piezoData.Piezo) / piezoADCRange) * piezoMicronRange;

%% - Trim to steady-state window 
% Look at the raw piezo trace (piezoData.ArduinoTime vs piezoData.Microns)
% once, find the clean steady-state region by eye, and set the two
% boundaries below directly - simpler and more reliable than automatic
% detection.
manualTrimStartMs = 7010000; 
manualTrimEndMs    = 7320000;  

piezoData = piezoData(piezoData.ArduinoTime > manualTrimStartMs & piezoData.ArduinoTime < manualTrimEndMs, :);
fprintf('Trimmed piezo data to manual window: [%d, %d] ms\n', manualTrimStartMs, manualTrimEndMs);

% Sanity check plot
figure('Name', 'Piezo trace after manual trimming', 'Color', 'w');
plot(piezoData.ArduinoTime, piezoData.Microns);
xlabel('Arduino time (ms)');
ylabel('Piezo position (microns)');
title('Piezo trace after trimming - check both ends look clean');
grid on;

%% Align piezo time onto frame timebase using sync pulses 
% NOTE: matching is done using the FULL sync pulse trains from both logs
% (not pre-truncated to the shorter log's count before searching for the
% best lag), since pre-truncation can discard genuinely matching pulses
% and bias the result.
syncPiezo = unique(piezoData.LastSyncPulsetime);
syncFrame = unique(frameData.LastSyncPulseTime);

d1 = diff(syncPiezo);
d2 = diff(syncFrame);

maxLag = 100;
bestLag = NaN; bestCorr = -Inf; bestLen = 0;
for lag = -maxLag:maxLag
    if lag < 0
        a = d1(-lag+1:end);
        b = d2(1:min(length(a), length(d2)));
    else
        b = d2(lag+1:end);
        a = d1(1:min(length(b), length(d1)));
    end
    L = min(length(a), length(b));
    if L < 50
        continue
    end
    a = a(1:L); b = b(1:L);
    R = corrcoef(a, b);
    corrVal = R(1, 2);
    if corrVal > bestCorr
        bestCorr = corrVal;
        bestLag = lag;
        bestLen = L;
    end
end
fprintf('Sync alignment: best_lag=%d, corr=%.6f, matched_len=%d\n', bestLag, bestCorr, bestLen);

if bestLag < 0
    t1 = syncPiezo(-bestLag+1 : -bestLag+bestLen+1);
    t2 = syncFrame(1 : bestLen+1);
else
    t1 = syncPiezo(1 : bestLen+1);
    t2 = syncFrame(bestLag+1 : bestLag+bestLen+1);
end

piezoData.AlignedTime = interp1(t1, t2, piezoData.ArduinoTime, 'linear', 'extrap');
piezoData = sortrows(piezoData, 'AlignedTime');

%% Estimate true period via autocorrelation 
t = piezoData.AlignedTime;
y = piezoData.Microns;
dt = median(diff(t));

sig = y - mean(y, 'omitnan');
sig(isnan(sig)) = 0;

maxLagMs = 300;
maxLagSamples = round(maxLagMs / dt);
[acFull, lags] = xcorr(sig, maxLagSamples, 'coeff');
acPositive = acFull(lags >= 0);

[~, pkLocs] = findpeaks(acPositive(2:end));
if isempty(pkLocs)
    error('No clear periodicity found.');
end
periodSamples = pkLocs(1) + 1;
periodMs = periodSamples * dt;
fprintf('Estimated true volume period: %.2f ms (%.2f Hz)\n', periodMs, 1000/periodMs);

%% Find troughs (true start-of-ramp points) 
[~, troughLocs] = findpeaks(-sig, 'MinPeakDistance', round(periodSamples * 0.8));
troughTimes = t(troughLocs);
fprintf('Found %d candidate cycle troughs.\n', numel(troughTimes));

%% Build average cycle shape (to locate ramp peak / fly-back boundary) 
timeX = 0:0.1:periodMs;
allCycles = nan(numel(timeX), numel(troughLocs) - 1);
for c = 1:numel(troughLocs) - 1
    idxStart = troughLocs(c);
    idxEnd   = troughLocs(c) + periodSamples;
    if idxEnd > height(piezoData), continue; end
    tSeg = t(idxStart:idxEnd) - t(idxStart);
    ySeg = y(idxStart:idxEnd);
    if numel(tSeg) > 1
        allCycles(:, c) = interp1(tSeg, ySeg, timeX, 'linear', 'extrap');
    end
end
avgCycle = mean(allCycles, 2, 'omitnan');
[~, peakIdx] = max(avgCycle);
peakTimeMs = timeX(peakIdx);
fprintf('Ramp peak at %.2f ms; fly-back window = [%.2f, %.2f] ms\n', peakTimeMs, peakTimeMs, periodMs);

%% Assign each frame its plane index and phase within the cycle 
numFrames = height(frameData);
planeOffset = 1;   % shift applied so planes 0,1,2 are correctly assigned as fly-back TODO: need to double check how piezo is logged and fix this
frameData.PlaneIdx = mod((0:numFrames-1) + planeOffset, numPlanes)';

frameTimes = frameData.TwoPFrameTime;
troughTimesSorted = sort(troughTimes);
phase = nan(numFrames, 1);
for i = 1:numFrames
    ft = frameTimes(i);
    idx = find(troughTimesSorted <= ft, 1, 'last');
    if ~isempty(idx)
        phase(i) = ft - troughTimesSorted(idx);
    end
end
frameData.PhaseMs = phase;

%%  Compute per-plane average trajectory (needed for boundary test + final figure) ----
frameDurationMs60 = 1000 / 60;
timeXFrame = 0:0.1:frameDurationMs60;
avgShapePerPlaneRaw = nan(numel(timeXFrame), numPlanes);

for p = 0:numPlanes - 1
    planeFrameIdx = find(frameData.PlaneIdx == p);
    segments = nan(numel(timeXFrame), numel(planeFrameIdx));

    for i = 1:numel(planeFrameIdx)
        tStart = frameData.TwoPFrameTime(planeFrameIdx(i));
        tEnd   = tStart + frameDurationMs60;

        mask = piezoData.AlignedTime >= tStart & piezoData.AlignedTime < tEnd;
        seg     = piezoData.Microns(mask);
        segTime = piezoData.AlignedTime(mask) - tStart;

        if numel(seg) > 1
            segments(:, i) = interp1(segTime, seg, timeXFrame, 'linear', 'extrap');
        end
    end

    avgShapePerPlaneRaw(:, p + 1) = mean(segments, 2, 'omitnan');
end

%% Final figure: per-plane average trajectory, labelled with the final fly-back result 
avgShapePerPlane = avgShapePerPlaneRaw;

colors = turbo(numPlanes);
figA=figure('Name', 'Average piezo trajectory per plane', 'Color', 'w');
hold on;

for p = 0:numPlanes - 1
    labelSuffix = '';
    if ismember(p, flybackPlanes)
        labelSuffix = ' (fly-back)';
    end
    plot(timeXFrame, avgShapePerPlane(:, p + 1), ...
        'Color', colors(p + 1, :), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('Plane %d%s', p, labelSuffix));
end

xlabel('Time within frame (ms, 60 Hz)');
ylabel('Piezo position (microns)');
title('Average piezo trajectory per plane');
legend('show', 'Location', 'eastoutside');
grid on;
box on;
defaultAxesProperties(gca, 1);

%%
baseFileName = sprintf('piezo_flybackexample_m25133_20260219');
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\Methods\Fig2.5_zmotioncCorrection\piezo_flyback';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = fullfile(outputDir, baseFileName);
saveFigureFormats(figA, fullSavePath);