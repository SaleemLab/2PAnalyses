
%% 
load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_PeripheralData_M25132_GrayScreen_20260226_00001.mat")
t = peripheralData.Wheel.sampleTimes;
v = peripheralData.Wheel.Value;

fs = 60; 
tickToCm = 3.1415 * 20 / 1024; 
wheelSpeed = [0; diff(peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(peripheralData.Wheel.sampleTimes)];
wheelSpeed(abs(wheelSpeed) > 150) = NaN; 

mIdx = find(wheelSpeed > 1 & ~isnan(wheelSpeed)); 
isRunning = zeros(size(v));
isRunning(mIdx) = 1;

figure('Color', 'w', 'Name', 'Running vs Stationary Check', 'Position', [50 100 1200 400]);
hold on;

runStarts = find(diff([0; isRunning]) == 1);
runEnds   = find(diff([isRunning; 0]) == -1);

for i = 1:length(runStarts)
    patch([t(runStarts(i)) t(runEnds(i)) t(runEnds(i)) t(runStarts(i))], ...
          [-10 -10 150 150], [0.8 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

%
patch([t(1) t(end) t(end) t(1)], [-10 -10 150 150], [1 0.8 0.8], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

plot(t, v, 'k', 'LineWidth', 1, 'DisplayName', 'Wheel Speed');
yline(1, 'r--', 'Threshold (1 cm/s)', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Speed (cm/s)');
title('Green Shading = Moving (mIdx) | Red Shading = Stationary (sIdx)');
legend('Running Bouts');
xlim([t(1) t(1)+120]); % Zoom into the first 2 minutes for a clear view
ylim([-2 40]);         % Adjust Y-axis to see the low-speed transitions
defaultAxesProperties(gca, true);

%% Including running and stationary (manually change files from darkness to gray)
load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_GrayScreen_20260226_00001.mat")
sigIdx = find(response.tuningCurve.dFFNeuropilCorrected.isSignificant_999); % significant using running and stationanry bins 
nSig = length(sigIdx);
nTotal = length(response.tuningCurve.dFFNeuropilCorrected.statMean);
stat  = response.tuningCurve.dFFNeuropilCorrected.statMean(sigIdx);
move  = response.tuningCurve.dFFNeuropilCorrected.moveMean(sigIdx, :);

allBins = [stat, move]; 
maxR = max(allBins, [], 2);
minR = min(allBins, [], 2);

modulation_ratio = 2 * (maxR - minR) ./ (maxR + minR);
modulation = modulation_ratio * 100;

figure('Color', 'w', 'Position', [100 100 500 400]);
h = histogram(modulation, 0:10:200, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'w');
xlabel('Modulation (%)');
ylabel('Number of ROIs');
title(sprintf('Session Modulation (n=%d/%d Sig ROIs)', nSig, nTotal));
defaultAxesProperties(gca, true);

targets = [30, 100, 180]; 
labels = {'Low modulation', 'Mid modulation', 'High modulation (Switch)'};

figure('Color', 'w', 'Position', [100 100 1000 350], 'Name', 'ROI Examples');
tlo = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

for i = 1:length(targets)
    [~, bestMatch] = min(abs(modulation - targets(i)));
    roiIdx = sigIdx(bestMatch);
    currentMod = modulation(bestMatch);
    
    yStat = response.tuningCurve.dFFNeuropilCorrected.statMean(roiIdx);
    yMove = response.tuningCurve.dFFNeuropilCorrected.moveMean(roiIdx, :);
    
    edges = response.tuningCurve.speedBins;
    xMove = edges(1:end-1) + diff(edges)/2; 
    xStat = 0;
    
    ax = nexttile(tlo);
    hold(ax, 'on');
    
    plot(ax, xMove, yMove, 'o-k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    plot(ax, xStat, yStat, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    title(ax, sprintf('%s\nROI %d: %.1f%% Mod', labels{i}, roiIdx, currentMod));
    xlabel(ax, 'Speed (cm/s)');
    ylabel(ax, 'dFF (Neuropil Corr)');
    
    %defaultAxesProperties(ax, true);
    yAll = [yStat, yMove];
    ylim(ax, [min(0, min(yAll)*1.1), max(yAll)*1.2]);
end
%% Only including running bins: (manually change files from darkness to gray) 
load("\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\DATA\SUBJECTS\M26004\Analysis\20260321\M26004_20260321_Response_M26004_GrayScreen_20260321_00001.mat")
sigIdx = find(response.tuningCurve.dFFNeuropilCorrected.isSignificantMoving_999); % significant using running bins only
nSig = length(sigIdx);
nTotal = length(response.tuningCurve.dFFNeuropilCorrected.statMean);
move  = response.tuningCurve.dFFNeuropilCorrected.moveMean(sigIdx, :);


maxR = max(move, [], 2);
minR = min(move, [], 2);

modulation_ratio = 2 * (maxR - minR) ./ (maxR + minR);
modulation = modulation_ratio * 100;

figure('Color', 'w', 'Position', [100 100 500 400]);
h = histogram(modulation, 0:5:100, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'w');
xlabel('% Modulation');
ylabel('Number of ROIs');
title(sprintf('Session Modulation (n=%d/%d Sig ROIs)', nSig, nTotal));

mVal = median(modulation);

defaultAxesProperties(gca, true);

targets = [60, 85, 100]; 
labels = {'Low modulation', 'Mid modulation', 'High modulation (Switch)'};

figure('Color', 'w', 'Position', [100 100 1000 350], 'Name', 'ROI Examples');
tlo = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

for i = 1:length(targets)
    [~, bestMatch] = min(abs(modulation - targets(i)));
    roiIdx = sigIdx(bestMatch);
    currentMod = modulation(bestMatch);
    
    yStat = response.tuningCurve.dFFNeuropilCorrected.statMean(roiIdx);
    yMove = response.tuningCurve.dFFNeuropilCorrected.moveMean(roiIdx, :);
    
    edges = response.tuningCurve.speedBins;
    xMove = edges(1:end-1) + diff(edges)/2; 
    xStat = 0;
    
    ax = nexttile(tlo);
    hold(ax, 'on');
    
    plot(ax, xMove, yMove, 'o-k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    % plot(ax, xStat, yStat, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    title(ax, sprintf('%s\nROI %d: %.1f%% Mod', labels{i}, roiIdx, currentMod));
    xlabel(ax, 'Speed (cm/s)');
    ylabel(ax, 'dFF (Neuropil Corr)');
    
    %defaultAxesProperties(ax, true);
    yAll = [yStat, yMove];
    ylim(ax, [min(0, min(yAll)*1.1), max(yAll)*1.2]);
end

%% % modulation gray vs darkness 
% Load
res = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_GrayScreen_20260226_00001.mat"); 
respGray = res.response; 
res2 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_Darkness_20260226_00001.mat");  
respDark = res2.response; 

sigDark = find(respDark.tuningCurve.dFFNeuropilCorrected.isSignificant_999);
sigGray = find(respGray.tuningCurve.dFFNeuropilCorrected.isSignificant_999);
allSigROIs = union(sigDark, sigGray);
bothIdx = intersect(sigDark, sigGray);
onlyDarkIdx = setdiff(sigDark, sigGray);
onlyGrayIdx = setdiff(sigGray, sigDark);

% Function uses the 2 * (max-min)/(max+min) ratio and converts to % (0-200%)
getMod = @(stat, move) 100 * (2 * (max([stat, move],[],2) - min([stat, move],[],2)) ./ ...
                                 (max([stat, move],[],2) + min([stat, move],[],2)));

figure('Color', 'w', 'Name', 'Speed Modulation Comparison');
hold on;

% Both Significant
mD_both = getMod(respDark.tuningCurve.dFFNeuropilCorrected.statMean(bothIdx), respDark.tuningCurve.dFFNeuropilCorrected.moveMean(bothIdx,:));
mG_both = getMod(respGray.tuningCurve.dFFNeuropilCorrected.statMean(bothIdx), respGray.tuningCurve.dFFNeuropilCorrected.moveMean(bothIdx,:));
scatter(mD_both, mG_both, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.7);

% Dark Only
mD_onlyD = getMod(respDark.tuningCurve.dFFNeuropilCorrected.statMean(onlyDarkIdx), respDark.tuningCurve.dFFNeuropilCorrected.moveMean(onlyDarkIdx,:));
mG_onlyD = getMod(respGray.tuningCurve.dFFNeuropilCorrected.statMean(onlyDarkIdx), respGray.tuningCurve.dFFNeuropilCorrected.moveMean(onlyDarkIdx,:));
scatter(mD_onlyD, mG_onlyD, 40, [0.5 0.5 0.5], 'LineWidth', 1);

% Gray Only
mD_onlyG = getMod(respDark.tuningCurve.dFFNeuropilCorrected.statMean(onlyGrayIdx), respDark.tuningCurve.dFFNeuropilCorrected.moveMean(onlyGrayIdx,:));
mG_onlyG = getMod(respGray.tuningCurve.dFFNeuropilCorrected.statMean(onlyGrayIdx), respGray.tuningCurve.dFFNeuropilCorrected.moveMean(onlyGrayIdx,:));
scatter(mD_onlyG, mG_onlyG, 40, 'b', 'LineWidth', 1);

% 
plot([0 200], [0 200], 'k--', 'HandleVisibility', 'off'); 
axis square; 
xlim([0 200]); ylim([0 200]);

xlabel('Modulation % (Darkness)', 'FontSize', 14); 
ylabel('Modulation % (Gray Screen)', 'FontSize', 14);
legend('Significant in Both', 'Significant in Dark Only', 'Significant in Gray Only', 'Location', 'southeast');
legend boxoff;
grid off;

defaultAxesProperties(gca, true)
set(gca, 'FontSize', 12);
title('Session Modulation: Gray vs darkness ');


%% RUnning only 
 
% Load
res = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_GrayScreen_20260226_00001.mat"); 
respGray = res.response; 
res2 = load("Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260226\M25132_20260226_Response_M25132_Darkness_20260226_00001.mat");  
respDark = res2.response; 

sigDark = find(respDark.tuningCurve.dFFNeuropilCorrected.isSignificantMoving_999);
sigGray = find(respGray.tuningCurve.dFFNeuropilCorrected.isSignificantMoving_999);
allSigROIs = union(sigDark, sigGray);
bothIdx = intersect(sigDark, sigGray);
onlyDarkIdx = setdiff(sigDark, sigGray);
onlyGrayIdx = setdiff(sigGray, sigDark);

% Function uses the 2 * (max-min)/(max+min) ratio and converts to % (0-200%)
getMod_run = @(move) 100 * (2 * (max(move,[],2) - min(move,[],2)) ./ ...
                                 (max(move,[],2) + min(move,[],2)));

figure('Color', 'w', 'Position', [100 100 600 600]);
hold on;

mD_both = getMod_run(respDark.tuningCurve.dFFNeuropilCorrected.moveMean(bothIdx,:));
mG_both = getMod_run(respGray.tuningCurve.dFFNeuropilCorrected.moveMean(bothIdx,:));
scatter(mD_both, mG_both, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.6);

mD_onlyD = getMod_run(respDark.tuningCurve.dFFNeuropilCorrected.moveMean(onlyDarkIdx,:));
mG_onlyD = getMod_run(respGray.tuningCurve.dFFNeuropilCorrected.moveMean(onlyDarkIdx,:));
scatter(mD_onlyD, mG_onlyD, 40, [0.6 0.6 0.6], 'LineWidth', 1);

mD_onlyG = getMod_run(respDark.tuningCurve.dFFNeuropilCorrected.moveMean(onlyGrayIdx,:));
mG_onlyG = getMod_run(respGray.tuningCurve.dFFNeuropilCorrected.moveMean(onlyGrayIdx,:));
scatter(mD_onlyG, mG_onlyG, 40, 'b', 'LineWidth', 1);

plot([0 200], [0 200], 'k--', 'HandleVisibility', 'off'); 
axis square; 
xlim([0 100]); ylim([0 100]);

xlabel('Modulation % (Darkness)'); 
ylabel('Modulation % (Gray Screen)');
title('Session Modulation: Gray vs darkness');
legend('Significant Both', 'Dark Only', 'Gray Only', 'Location', 'southeast');
legend boxoff;
defaultAxesProperties(gca, true);
set(gca, 'FontSize', 12);

% --- Example Picking ---
pVals = [25, 50, 75]; 
targetsG = prctile(mG_both, pVals);
titles = {'25th Percentile (Low Mod)', '50th Percentile (Median)', '75th Percentile (High Mod)'};

figure('Color', 'w', 'Position', [100 100 1200 600]);
tlo = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'tight');

for i = 1:length(pVals)
    [~, bestMatchIdx] = min(abs(mG_both - targetsG(i)));
    roiIdx = bothIdx(bestMatchIdx); 
    
    yD = respDark.tuningCurve.dFFNeuropilCorrected.moveMean(roiIdx, :);
    yG = respGray.tuningCurve.dFFNeuropilCorrected.moveMean(roiIdx, :);
    
    edgesD = respDark.tuningCurve.speedBins;
    edgesG = respGray.tuningCurve.speedBins;
    xD = edgesD(1:end-1) + diff(edgesD)/2;
    xG = edgesG(1:end-1) + diff(edgesG)/2;
    
    axD = nexttile(tlo, i);
    plot(axD, xD, yD, 'o-k', 'MarkerFaceColor', 'k', 'LineWidth', 1.2);
    title(axD, sprintf('%s (ROI %d)\nDark Mod: %.1f%%', titles{i}, roiIdx, mD_both(bestMatchIdx)));
    ylabel(axD, 'Dark dFFNeu');
    set(axD, 'XScale', 'log');
    % defaultAxesProperties(axD, true);
    
    axG = nexttile(tlo, i+3);
    plot(axG, xG, yG, 'o-b', 'MarkerFaceColor', 'b', 'LineWidth', 1.2);
    title(axG, sprintf('Gray Mod: %.1f%%', mG_both(bestMatchIdx)));
    ylabel(axG, 'Gray dFFNeu');
    xlabel(axG, 'Speed (cm/s)');
    set(axG, 'XScale', 'log');
    % defaultAxesProperties(axG, true);
    
    allX = [xD, xG];
    xlim(axD, [0 max(allX)]);
    xlim(axG, [0 max(allX)]);
end