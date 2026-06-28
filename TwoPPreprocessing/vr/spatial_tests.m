% trying SI 
baseTrials = response.trialIndicesByCondition.Baseline;
numLaps = size(response.lapPosition2PFrameIdx, 1);
numBins = size(response.lapPosition2PFrameIdx, 2);

occupancyMatrix = zeros(numLaps, numBins);
for iLap = 1:numLaps
    for iBin = 1:numBins
        occupancyMatrix(iLap, iBin) = numel(response.lapPosition2PFrameIdx{iLap, iBin});
    end
end

% then sum across baseline laps
occupancy = sum(occupancyMatrix(baseTrials, :), 1);  % 1 x 200
p = occupancy / sum(occupancy);


% mean activity across baseline laps per cell (numROIs x numBins)
baseLapActivity = response.lapPositionActivity.dFFNeuropilCorrected(:, baseTrials,:); % (numROIs x laps x position) 
meanRate = squeeze(mean(baseLapActivity, 2, 'omitnan'));  % numROIs x 200

numROIs= size(baseLapActivity, 1)

SI = NaN(numROIs, 1);
for iCell = 1:numROIs
    r = meanRate(iCell, :);
    if all(isnan(r)), continue; end
    
    % clamp negative values to 0 (log2 of negative = complex)
    r = max(0, r);
    
    R = sum(p .* r, 'omitnan');
    if R == 0, continue; end
    
    terms = p .* (r / R) .* log2(r / R + eps);
    SI(iCell) = real(sum(terms, 'omitnan'));  % take real part as safety net
end

histogram(SI, 50);
xlabel('Spatial Information (bits/spike)');
ylabel('Count');
title('SI distribution');

% --- Define filters (all must be 477 x 1 logical) ---
excludeEdge = response.SMI_Metrics.dFFNeuropilCorrected.ExcludeEdgePeakCells;
highSI = SI > 0.5;
stableExpVar = crossValExpVar.dFFNeuropilCorrected.meanExpVar > 0.3;
sigExpVar = crossValExpVar.dFFNeuropilCorrected.pValues < 0.01;
sigOddEven = lapCorr_OddEven.rho < 0.8;

% check sizes match
fprintf('excludeEdge: %d, highSI: %d, stableExpVar: %d, sigExpVar: %d\n', ...
    numel(excludeEdge), numel(highSI), numel(stableExpVar), numel(sigExpVar));

% --- Combine all filters ---
keepCells = ~excludeEdge & highSI & stableExpVar & sigExpVar;

% --- Apply to full SMI vector (477 x 1) ---
allSMI = response.SMI_Metrics.dFFNeuropilCorrected.SMI;  % full 477 x 1
smiFiltered = allSMI(keepCells);

fprintf('Total cells: %d\n', numel(allSMI));
fprintf('Cells passing all filters: %d\n', sum(keepCells));
fprintf('Median SMI (filtered): %.3f\n', median(smiFiltered, 'omitnan'));

figure;

%All cells (no filtering)


% --- Filter without SI threshold ---
keepCells_noSI = ~excludeEdge & stableExpVar & sigExpVar;
allSMI_noSI = response.SMI_Metrics.dFFNeuropilCorrected.SMI(keepCells_noSI);
medianNoSI = median(allSMI_noSI, 'omitnan');

subplot(1,2,1);
histogram(allSMI_noSI, 30);
xlabel('SMI'); ylabel('Count');
title(sprintf('SMI - no SI filter (n=%d)', sum(keepCells_noSI)));
xline(medianNoSI, 'r-', sprintf('median=%.3f', medianNoSI), 'LineWidth', 2);

%Filtered cells 
medianFiltered = median(smiFiltered, 'omitnan');

subplot(1,2,2);
histogram(smiFiltered, 30);
xlabel('SMI'); ylabel('Count');
title(sprintf('SMI - filtered cells (n=%d)', sum(keepCells)));
xline(medianFiltered, 'r-', sprintf('median=%.3f', medianFiltered), 'LineWidth', 2);

%
prefBins = response.SMI_Metrics.dFFNeuropilCorrected.RpBin(keepCells);

% plot distribution of preferred bins
figure;
histogram(prefBins, 'BinEdges', 0:5:200);
xlabel('Position (cm)');
ylabel('Count');
title('Preferred landmark bin distribution (filtered cells)');
xline(40, 'r--', '40cm'); xline(80, 'r--', '80cm');
xline(120, 'r--', '120cm'); xline(160, 'r--', '160cm');



%
% get indices of putative position coders
positionCoders = keepCells_noSI & (response.SMI_Metrics.dFFNeuropilCorrected.SMI > 0.3);
positionCoderIdx = find(positionCoders);


meanOdd  = squeeze(mean(baseLapActivity(:, 1:2:end, :), 2, 'omitnan'));  % numROIs x numBins
meanEven = squeeze(mean(baseLapActivity(:, 2:2:end, :), 2, 'omitnan'));  % numROIs x numBins

figure;
nCols = 7;
nRows = ceil(numel(positionCoderIdx) / nCols);

for i = 1:numel(positionCoderIdx)
    iCell = positionCoderIdx(i);
    
    trainTrace = meanOdd(iCell, :);
    testTrace = meanEven(iCell, :);
    
    subplot(nRows, nCols, i); hold on;
    plot(trainTrace, 'b-', 'LineWidth', 1.5);
    plot(testTrace, 'r-', 'LineWidth', 1.5);
    
    % mark landmark positions
    xline(40, 'k--'); xline(80, 'k--');
    xline(120, 'k--'); xline(160, 'k--');
    
    smi = response.SMI_Metrics.dFFNeuropilCorrected.SMI(iCell);
    title(sprintf('ROI %d | SMI=%.2f', iCell, smi), 'FontSize', 7);
    
    if i == 1
        legend('Odd', 'Even', 'Location', 'best', 'FontSize', 6);
    end
end

sgtitle('Putative position coders (SMI > 0.3)');


%
% autocorrelation - include all 4 landmarks (40, 80, 120, 160cm)
binRange = 30:170;  % includes all landmarks, excludes unreliable edges

allACF = NaN(numFiltered, 2*length(binRange)-1);

for i = 1:numFiltered
    iCell = filteredIdx(i);
    meanTrace = meanTuning(iCell, binRange);
    if all(isnan(meanTrace)), continue; end
%     demeaned = meanTrace - mean(meanTrace, 'omitnan');
    allACF(i, :) = xcorr(meanTrace, 'normalized');
end

[~, lags] = xcorr(zeros(1, length(binRange)), 'normalized');

meanACF = mean(allACF, 1, 'omitnan');

figure;
plot(lags, meanACF, 'k-', 'LineWidth', 2);
xline(0, 'k--');
xline(40,  'r--', '40cm');  xline(-40,  'r--');
xline(80,  'b--', '80cm');  xline(-80,  'b--');
xline(120, 'r--', '120cm'); xline(-120, 'r--');
xline(160, 'b--', '160cm'); xline(-160, 'b--');
xlabel('Lag (position bins)');
ylabel('Autocorrelation');
title(sprintf('Mean ACF - 30 to 170cm (n=%d)', numFiltered));
xlim([-200 200]);
grid on;


%%
figure;
subplot(1,3,1);
plot(lags, mean(allACF(prefersAdjacent,:), 1, 'omitnan'), 'k-', 'LineWidth', 2);
xline(0, 'k--');
xline(40,  'r--', '40cm');  xline(-40,  'r--');
xline(80,  'b--', '80cm');  xline(-80,  'b--');
xline(120, 'r--', '120cm'); xline(-120, 'r--');
ylim([-0.5 0.4]); xlim([-140 140]);
xlabel('Lag (bins)'); ylabel('ACF');
title(sprintf('Fires every landmark\n(n=%d)', sum(prefersAdjacent)));
grid on;

subplot(1,3,2);
plot(lags, mean(allACF(prefersSameType,:), 1, 'omitnan'), 'k-', 'LineWidth', 2);
xline(0, 'k--');
xline(40,  'r--', '40cm');  xline(-40,  'r--');
xline(80,  'b--', '80cm');  xline(-80,  'b--');
xline(120, 'r--', '120cm'); xline(-120, 'r--');
ylim([-0.5 0.4]); xlim([-140 140]);
xlabel('Lag (bins)'); ylabel('ACF');
title(sprintf('Fires every other landmark\n(n=%d)', sum(prefersSameType)));
grid on;

subplot(1,3,3);
plot(lags, mean(allACF(ambiguous,:), 1, 'omitnan'), 'k-', 'LineWidth', 2);
xline(0, 'k--');
xline(40,  'r--', '40cm');  xline(-40,  'r--');
xline(80,  'b--', '80cm');  xline(-80,  'b--');
xline(120, 'r--', '120cm'); xline(-120, 'r--');
ylim([-0.5 0.4]); xlim([-140 140]);
xlabel('Lag (bins)'); ylabel('ACF');
title(sprintf('Ambiguous\n(n=%d)', sum(ambiguous)));
grid on;




%%
% get indices of every other landmark cells
everyOtherIdx = filteredIdx(prefersSameType);

figure;
nCols = 4;
nRows = ceil(numel(everyOtherIdx) / nCols);

for i = 1:numel(everyOtherIdx)
    iCell = everyOtherIdx(i);
    
    subplot(nRows, nCols, i); hold on;
    
    % plot mean tuning curve
    plot(meanTuning(iCell, :), 'k-', 'LineWidth', 2);
    
    % mark landmark positions
    xline(40,  'r--', 'A'); 
    xline(120, 'r--', 'A');
    xline(80,  'b--', 'B'); 
    xline(160, 'b--', 'B');
    
    % mark excluded edges
    xline(30,  'k:', 'start');
    xline(170, 'k:', 'end');
    
    acf80 = allACF(prefersSameType(i), lags == 80);
    acf40 = allACF(prefersSameType(i), lags == 40);
    
    title(sprintf('ROI %d\nACF40=%.2f ACF80=%.2f', iCell, acf40, acf80), 'FontSize', 7);
    xlabel('Position (bins)'); 
    ylabel('Activity');
    xlim([1 200]);
end

sgtitle('Every other landmark cells (80cm periodic, n=13)');

%
sgtitle('Mean ACF by spatial periodicity group');
% mean response in window around each landmark
lm1 = mean(meanTuning(filteredIdx, 30:50), 2, 'omitnan');   % 40cm  - type A
lm2 = mean(meanTuning(filteredIdx, 70:90), 2, 'omitnan');   % 80cm  - type B
lm3 = mean(meanTuning(filteredIdx, 110:130), 2, 'omitnan'); % 120cm - type A
lm4 = mean(meanTuning(filteredIdx, 150:170), 2, 'omitnan'); % 160cm - type B

% pairwise correlations across cells
fprintf('=== Same type ===\n');
fprintf('A-A (lm1 vs lm3, 40 vs 120cm): r = %.3f\n', corr(lm1, lm3));
fprintf('B-B (lm2 vs lm4, 80 vs 160cm): r = %.3f\n', corr(lm2, lm4));

fprintf('=== Different type ===\n');
fprintf('A-B (lm1 vs lm2, 40 vs 80cm):  r = %.3f\n', corr(lm1, lm2));
fprintf('A-B (lm1 vs lm4, 40 vs 160cm): r = %.3f\n', corr(lm1, lm4));
fprintf('A-B (lm2 vs lm3, 80 vs 120cm): r = %.3f\n', corr(lm2, lm3));
fprintf('A-B (lm3 vs lm4, 120 vs 160cm):r = %.3f\n', corr(lm3, lm4));



%%
lag40idx = lags == 40;
lag80idx = lags == 80;
lag120idx = lags == 120;

% extract ACF at landmark lags for each cell
acf_lag40  = allACF(:, lag40idx);   % adjacent landmark
acf_lag80  = allACF(:, lag80idx);   % same type landmark
acf_lag120 = allACF(:, lag120idx);  % 3rd landmark

% test if each is significantly above 0 across population
[p40,  ~, stats40]  = signrank(acf_lag40);
[p80,  ~, stats80]  = signrank(acf_lag80);
[p120, ~, stats120] = signrank(acf_lag120);

fprintf('ACF at lag 40:  median=%.3f, p=%.4f\n', median(acf_lag40),  p40);
fprintf('ACF at lag 80:  median=%.3f, p=%.4f\n', median(acf_lag80),  p80);
fprintf('ACF at lag 120: median=%.3f, p=%.4f\n', median(acf_lag120), p120);
%%
figure;
data = [acf_lag40, acf_lag80, acf_lag120];
labels = {'Lag 40\newline(adjacent)', 'Lag 80\newline(same type)', 'Lag 120'};

boxplot(data, 'Labels', labels);
hold on;
yline(0, 'k--');

% add significance markers
plot(1, max(acf_lag40)+0.05, 'k*');
plot(2, max(acf_lag80)+0.05, 'k*');
text(3, max(acf_lag120)+0.05, 'n.s.', 'HorizontalAlignment', 'center');

ylabel('Spatial autocorrelation');
title('ACF at landmark lags (n=107 cells)');
grid on;