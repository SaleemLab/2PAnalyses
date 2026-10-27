%% RSP bouton running-speed tuning in darkness (Edd's template classification)
%
% Selection:        circular-shift test on moving frames (pValMoving; Saleem et al., 2013)
% Shape:            fitGaussianTemplates_tuning (Horrocks et al., 2024), unchanged from
%                   the visual speed analysis (customOffset = 0.5), fitted to the mean
%                   response in each equal-occupancy speed bin. No R2 threshold.
% Preferred speed:  fitted template centre (band-pass/trough only), bin rank -> cm/s
%
% Requires: fitGaussianTemplates_tuning.m, plotAllSpeedTuningCategories_Edd.m

clear; clc;

%% Parameters
pairs = struct();
pairs.M25132 = {'20260219','20260223','20260226','20260228','20260303','20260313','20260306'};
pairs.M25133 = {'20260219','20260223','20260221'};
pairs.M26003 = {'20260316','20260322','20260324','20260325'};

targetStruct          = 'tuningCurve';
useField              = 'dFFNeuropilCorrected';
pval_shuffleThreshold = 0.01;   % set to 0.05 to match the visual selection threshold
minBins               = 5;
customOffset          = 0.5;    % same as the visual speed analysis
onlyBoutons           = true;
binFraction           = 0.07;   % Saleem et al. (2013) occupancy per bin

saveToResponseFiles = false;    % true = store classification_Edd in each session's response file
makeQCpdf           = false;    % true = write a PDF of all speed-tuned tuning curves

saveFolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section3_Fig4_6\';

% Edd's codes: 1 = low-pass, 2 = high-pass, 3 = band-pass, 4 = trough
codeNames = {'lowpass', 'highpass', 'bandpass', 'trough_inverted'};

colors.highpass        = [0.89, 0.10, 0.11];   % red
colors.bandpass        = [0.99, 0.60, 0.17];   % orange
colors.lowpass         = [0.12, 0.47, 0.71];   % blue
colors.trough_inverted = [0.20, 0.63, 0.17];   % green

tuningTypesPlot = {'highpass', 'bandpass', 'lowpass', 'trough_inverted'};
faceColors_fig  = {colors.highpass, colors.bandpass, colors.lowpass, colors.trough_inverted};
labelNames      = {'High-pass', 'Band-pass', 'Low-pass', 'Trough-inv'};

%% Filter and load session info
filteredTable = filterMasterTable_usingNameSessionPairs('MousePairs', pairs, 'Exclude', 0, ...
    'HasStimulus', {'Darkness'}); % , 'GrayScreen'
mouseInfo = sessionsToProcess(filteredTable);

%% Accumulators
allSessionPcts_fig = [];
allTypes      = {};
allPrefSpeeds = [];
allR2s        = [];
allBinDurations = [];
sessionResults  = struct('mouse', {}, 'session', {}, 'E', {}, 'moveMean', {});

sessionLog = struct('Mouse', {}, 'Session', {}, 'Type', {}, 'Status', {}, ...
    'RecDur_s', {}, 'RunTime_s', {}, 'RunFrac', {}, 'nBins', {}, ...
    'BinMean_s', {}, 'BinMin_s', {}, 'BinMax_s', {}, 'Target7pct_s', {}, 'nROIs', {}, 'nSig', {});
newEntry = @(m, s, t, status) struct('Mouse', string(m), 'Session', string(s), ...
    'Type', string(t), 'Status', string(status), 'RecDur_s', NaN, 'RunTime_s', NaN, ...
    'RunFrac', NaN, 'nBins', NaN, 'BinMean_s', NaN, 'BinMin_s', NaN, 'BinMax_s', NaN, ...
    'Target7pct_s', NaN, 'nROIs', NaN, 'nSig', NaN);

%% Load and classify each session
for thisMouse = 1:size(mouseInfo, 1)
    mousenumber  = mouseInfo{thisMouse, 1};
    sessionNames = mouseInfo{thisMouse, 2};

    for thisSession = 1:numel(sessionNames)
        sessionName = sessionNames{thisSession};
        typeImaged  = "";

        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if isempty(infoPath) || ~isfile(infoPath)
            sessionLog(end+1) = newEntry(mousenumber, sessionName, "", "skipped: no session info file");
            continue;
        end
        matchRow = find(strcmp(filteredTable.MouseID, mousenumber) & strcmp(filteredTable.Session, sessionName), 1);
        if isempty(matchRow)
            sessionLog(end+1) = newEntry(mousenumber, sessionName, "", "skipped: not in metadata table"); 
            continue;
        end
        typeImaged = filteredTable.TypeImaged{matchRow};
        if onlyBoutons && ~strcmpi(typeImaged, 'Boutons')
            sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: not a bouton session"); 
            continue;
        end

        try
            sessionFileInfo = load(infoPath, 'sessionFileInfo').sessionFileInfo;
            stimNames = {sessionFileInfo.stimFiles.name};
            targetIdx = find(contains(stimNames, 'Darkness'));
            if isempty(targetIdx)
                sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: no darkness stimulus"); 
                continue;
            end

            for thisStim = 1:numel(targetIdx)
                thisStimName = stimNames{targetIdx(thisStim)};
                fileFullPath = fullfile(sessionFileInfo.Directories.save_folder, ...
                    sprintf('%s_%s_Response_%s.mat', mousenumber, sessionName, thisStimName));
                if ~isfile(fileFullPath)
                    sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: no response file");
                    continue;
                end

                resp = load(fileFullPath, 'response').response;
                if ~isfield(resp, targetStruct) || ~isfield(resp.(targetStruct), useField)
                    sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "skipped: no tuning curves"); 
                    continue;
                end
                tc = resp.(targetStruct);

                % --- Occupancy and session exclusion (before any pooling) ---
                occ   = tc.occupancy.moving;
                statT = tc.occupancy.stationary;
                runT  = sum(occ);
                recDur = statT + runT;   % excludes artefact frames (>150 cm/s)
                nBins  = numel(occ);

                entry = newEntry(mousenumber, sessionName, typeImaged, "");
                entry.RecDur_s     = recDur;
                entry.RunTime_s    = runT;
                entry.RunFrac      = runT / recDur;
                entry.nBins        = nBins;
                entry.BinMean_s    = mean(occ);
                entry.BinMin_s     = min(occ);
                entry.BinMax_s     = max(occ);
                entry.Target7pct_s = binFraction * recDur;

                if nBins < minBins
                    entry.Status = sprintf("excluded: only %d bin(s)", nBins);
                    sessionLog(end+1) = entry; %#ok<SAGROW>
                    continue;
                end

                % --- Selection + Edd's classification ---
                E = classifySessionEdd(tc, useField, pval_shuffleThreshold, customOffset);

                entry.Status = "included";
                entry.nROIs  = E.nROIs;
                entry.nSig   = numel(E.sigIdx);
                sessionLog(end+1) = entry; %#ok<SAGROW>
                allBinDurations = [allBinDurations; occ(:)]; %#ok<AGROW>

                sessionResults(end+1) = struct('mouse', mousenumber, 'session', sessionName, ...
                    'E', E, 'moveMean', tc.(useField).moveMean); %#ok<SAGROW>

                % Per-session percentages (of speed-tuned boutons), in tuningTypesPlot order
                sigTypes = cell(numel(E.code), 1);
                valid    = ~isnan(E.code);
                sigTypes(valid)  = codeNames(E.code(valid));
                sigTypes(~valid) = {'untuned'};
                if ~isempty(E.sigIdx)
                    pcts = cellfun(@(t) 100 * sum(strcmp(sigTypes, t)) / numel(E.sigIdx), tuningTypesPlot);
                    allSessionPcts_fig = [allSessionPcts_fig; pcts]; %#ok<AGROW>
                end

                % Pool: speed-tuned boutons with their shape, the rest as 'untuned'
                nUntunedHere = E.nROIs - numel(E.sigIdx);
                allTypes      = [allTypes; sigTypes; repmat({'untuned'}, nUntunedHere, 1)]; %#ok<AGROW>
                allPrefSpeeds = [allPrefSpeeds; E.pref; nan(nUntunedHere, 1)]; %#ok<AGROW>
                allR2s        = [allR2s; E.R2; nan(nUntunedHere, 1)]; %#ok<AGROW>

                % --- Optional: store classification in the response file ---
                if saveToResponseFiles
                    resp.(targetStruct).(useField).classification_Edd = E;
                    resp.(targetStruct).(useField).classification_Edd.settings = ...
                        struct('pval_shuffleThreshold', pval_shuffleThreshold, 'customOffset', customOffset);
                    response = resp; %#ok<NASGU>
                    save(fileFullPath, 'response', '-append');
                end
            end
        catch ME
            sessionLog(end+1) = newEntry(mousenumber, sessionName, typeImaged, "error: " + string(ME.message)); %#ok<SAGROW>
        end
    end
end

%% Session report
logTable = struct2table(sessionLog);
fprintf('\n===== SESSION LOG =====\n');
disp(logTable(:, {'Mouse','Session','Type','Status','RecDur_s','RunFrac','nBins','BinMean_s','BinMin_s','BinMax_s','Target7pct_s','nROIs','nSig'}));

inc = logTable.Status == "included";
fprintf('\nSessions listed: %d | included: %d | excluded (<%d bins): %d | skipped: %d\n', ...
    height(logTable), sum(inc), minBins, sum(startsWith(logTable.Status, "excluded")), ...
    sum(startsWith(logTable.Status, "skipped") | startsWith(logTable.Status, "error")));
fprintf('Mice included: %d (%s)\n', numel(unique(logTable.Mouse(inc))), strjoin(unique(logTable.Mouse(inc)), ', '));

fprintf('\n===== INCLUDED SESSIONS ONLY =====\n');
fprintf('Recording duration: mean %.1f s | range %.1f to %.1f s\n', ...
    mean(logTable.RecDur_s(inc)), min(logTable.RecDur_s(inc)), max(logTable.RecDur_s(inc)));
fprintf('Fraction of time running: mean %.2f | range %.2f to %.2f\n', ...
    mean(logTable.RunFrac(inc)), min(logTable.RunFrac(inc)), max(logTable.RunFrac(inc)));
fprintf('Speed bins per session: range %d to %d\n', min(logTable.nBins(inc)), max(logTable.nBins(inc)));
fprintf('Running time per bin: mean %.1f s | range %.1f to %.1f s\n', ...
    mean(allBinDurations), min(allBinDurations), max(allBinDurations));

%% Bouton summary
tuningTypes  = {'lowpass', 'trough_inverted', 'bandpass', 'highpass'};
tuningLabels = {'Low-pass', 'Trough-inverted', 'Band-pass', 'High-pass'};
totalROIsCounted = numel(allTypes);
nClassified      = sum(~strcmp(allTypes, 'untuned'));   % = all speed-tuned boutons
nUntuned         = sum(strcmp(allTypes, 'untuned'));

fprintf('\n===== BOUTONS =====\n');
fprintf('Total: %d\n', totalROIsCounted);
fprintf('Speed-tuned (p <= %.3f): %d (%.1f%%)\n', pval_shuffleThreshold, nClassified, 100 * nClassified / totalROIsCounted);
fprintf('Untuned: %d (%.1f%%)\n', nUntuned, 100 * nUntuned / totalROIsCounted);
for k = 1:numel(tuningTypes)
    n = sum(strcmp(allTypes, tuningTypes{k}));
    fprintf('  %-16s %d (%.1f%% of speed-tuned)\n', tuningLabels{k}, n, 100 * n / nClassified);
end
fprintf('\nPreferred speed (fitted centre; band-pass and trough only):\n');
for k = [2 3]
    sp = allPrefSpeeds(strcmp(allTypes, tuningTypes{k}));
    if ~isempty(sp)
        fprintf('  %-16s median %.1f | range %.1f to %.1f cm/s\n', tuningLabels{k}, median(sp), min(sp), max(sp));
    end
end

%% Save pooled results
if ~exist(saveFolder, 'dir'), mkdir(saveFolder); end
settings = struct('pval_shuffleThreshold', pval_shuffleThreshold, 'minBins', minBins, ...
    'customOffset', customOffset, 'useField', useField, 'targetStruct', targetStruct);
save(fullfile(saveFolder, 'runningSpeedTuning_darkness_Edd.mat'), 'allTypes', 'allPrefSpeeds', 'allR2s', ...
    'allSessionPcts_fig', 'logTable', 'sessionResults', 'settings', '-v7.3');

%% Figure A: Running speed and dF/F example
exampleBase  = 'Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\';
proc2PData   = load(fullfile(exampleBase, 'M26003_20260322_processed2PData_M26003_Darkness_20260322_00001.mat'));
peripheraleg = load(fullfile(exampleBase, 'M26003_20260322_PeripheralData_M26003_Darkness_20260322_00001.mat'));
respeg       = load(fullfile(exampleBase, 'M26003_20260322_Response_M26003_Darkness_20260322_00001.mat'));

% Same selection criterion as the analysis
isSigShuffle = find(respeg.response.tuningCurve.(useField).pValMoving <= pval_shuffleThreshold);
dFF          = proc2PData.processedSignals.(useField)(isSigShuffle, :);
timeVec      = proc2PData.TwoPFrameTime;

tickToCm   = pi * 20 / 1024;
wheelTimes = peripheraleg.peripheralData.Wheel.sampleTimes;
runSpeed   = [0; diff(peripheraleg.peripheralData.Wheel.Value * tickToCm)] ./ [1; diff(wheelTimes)];
runSpeed(abs(runSpeed) > 150) = NaN;

tStart = 47; tEnd = 360;
keep_dff   = timeVec >= tStart & timeVec <= tEnd;
timeVec    = timeVec(keep_dff);
dFF        = dFF(:, keep_dff);
keep_wheel = wheelTimes >= tStart & wheelTimes <= tEnd;
wheelTimes = wheelTimes(keep_wheel);
runSpeed   = runSpeed(keep_wheel);

[~, nearestIdx]  = min(abs(wheelTimes - timeVec'), [], 1);
runSpeed_matched = runSpeed(nearestIdx);
runSpeed_matched = runSpeed_matched(:);

corrs = corr(dFF', runSpeed_matched, 'rows', 'complete');
[~, sortOrder] = sort(corrs, 'ascend');
dFF = dFF(sortOrder, :);

dFF_norm = dFF - min(dFF, [], 2);
dFF_norm = dFF_norm ./ max(dFF_norm, [], 2);

figA = figure('Color', 'w', 'Position', [100 100 900 400]);
ax1  = subplot(2, 1, 1);
plot(wheelTimes, runSpeed, 'k', 'LineWidth', 0.8);
ylabel('Speed (cm s^{-1})', 'FontSize', 9);
title('RSP bouton activity in darkness', 'FontSize', 10);
set(ax1, 'Box', 'off', 'TickDir', 'out');

ax2 = subplot(2, 1, 2);
imagesc(timeVec, 1:size(dFF_norm, 1), dFF_norm);
colormap(ax2, flipud(gray));
ylabel('Boutons', 'FontSize', 9); xlabel('Time (s)', 'FontSize', 9);
set(ax2, 'Box', 'off', 'TickDir', 'out', 'YDir', 'normal');

cb = colorbar(ax2, 'eastoutside');
cb.Label.String = '\DeltaF/F (norm.)';

linkaxes([ax1 ax2], 'x');
drawnow;
pos1 = get(ax1, 'Position'); pos2 = get(ax2, 'Position');
set(ax2, 'Position', [pos1(1), pos2(2), pos1(3), pos2(4)]);
set(cb,  'Position', [pos1(1) + pos1(3) + 0.01, pos2(2), 0.02, pos2(4)]);

outputDir = fullfile(saveFolder, 'exampleWheel_dFF');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figA, fullfile(outputDir, 'wheelSpeed_dff_m26003_day5'));

%% Figure B: Gaussian fit examples (Edd's classification)
outputDir = fullfile(saveFolder, 'exampleFits');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

exampleFiles = {
    'Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260223\M25132_20260223_Response_M25132_Darkness_20260223_00001.mat', 'M25132_20260223_takebandpass';
    'Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260220\M25132_20260220_Response_M25132_Darkness_20260220_00001.mat', 'M25132_20260220_taketrough';
    'Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260322\M26003_20260322_Response_M26003_Darkness_20260322_00001.mat', 'M26003_20260322_takelowpass';
    'Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260316\M26003_20260316_Response_M26003_Darkness_20260316_00001.mat', 'M26003_20260316_takehighpass'};

for e = 1:size(exampleFiles, 1)
    if ~isfile(exampleFiles{e, 1})
        warning('Example file not found: %s', exampleFiles{e, 1});
        continue;
    end
    fits = load(exampleFiles{e, 1}, 'response');
    figB = plotAllSpeedTuningCategories_Edd(fits.response, useField, pval_shuffleThreshold, customOffset);
    saveFigureFormats(figB, fullfile(outputDir, exampleFiles{e, 2}));
end

%% Figure C: Category proportions with per-session dots
meanPct = mean(allSessionPcts_fig, 1);
semPct  = std(allSessionPcts_fig, 0, 1) / sqrt(size(allSessionPcts_fig, 1));

figC1 = figure('Position', [100, 100, 380, 360], 'Color', 'w');
ax_c1 = gca; hold on;
for k = 1:numel(tuningTypesPlot)
    bar(k, meanPct(k), 'FaceColor', faceColors_fig{k}, 'EdgeColor', 'none', 'BarWidth', 0.6);
end
errorbar(1:numel(tuningTypesPlot), meanPct, semPct, 'k.', 'LineWidth', 1, 'CapSize', 4);
for k = 1:numel(tuningTypesPlot)
    xJitter = k + (rand(size(allSessionPcts_fig, 1), 1) - 0.5) * 0.25;
    scatter(xJitter, allSessionPcts_fig(:, k), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
end

ax_c1.XTick              = 1:numel(tuningTypesPlot);
ax_c1.XTickLabel         = labelNames;
ax_c1.XTickLabelRotation = 30;
ax_c1.FontSize           = 8;
ax_c1.FontName           = 'Arial';
ax_c1.TickDir            = 'out';
ax_c1.LineWidth          = 0.75;
ax_c1.Position           = [0.18, 0.28, 0.75, 0.62];
box off;

ylabel('% of speed-tuned boutons', 'FontSize', 9, 'FontName', 'Arial');
title('RSP', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('%d speed-tuned of %d boutons | n = %d sessions', ...
    nClassified, totalROIsCounted, size(allSessionPcts_fig, 1)), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');
defaultAxesProperties(ax_c1, 0);
offsetAxes(ax_c1);

outputDir = fullfile(saveFolder, 'speedModulated_histogram_darkness');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figC1, fullfile(outputDir, 'category_counts_darkness'));

%% Figure D: Band-pass preferred speed histogram (fitted centre)
prefSpeedBins = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5];
binCentres    = [5, 10, 15, 20, 25, 30];
bp_speeds     = allPrefSpeeds(strcmp(allTypes, 'bandpass'));
histBP        = histcounts(bp_speeds, prefSpeedBins);

figC2 = figure('Position', [100, 100, 380, 360], 'Color', 'w');
ax_c2 = gca; hold on;
bar(binCentres, histBP, 'FaceColor', colors.bandpass, 'EdgeColor', 'none', 'BarWidth', 0.5);

xlim([prefSpeedBins(1), prefSpeedBins(end)]);
ax_c2.XTick      = binCentres;
ax_c2.XTickLabel = {'5', '10', '15', '20', '25', '30'};
ylim([0, max([histBP, 1]) * 1.2]);
ax_c2.FontSize   = 8;
ax_c2.FontName   = 'Arial';
ax_c2.TickDir    = 'out';
ax_c2.LineWidth  = 0.75;
ax_c2.Layer      = 'top';
ax_c2.Position   = [0.18, 0.22, 0.75, 0.65];
box off;

xlabel('Preferred speed (cm s^{-1})', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Number of boutons', 'FontSize', 9, 'FontName', 'Arial');
title('RSP: band-pass', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('n = %d boutons', numel(bp_speeds)), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');
defaultAxesProperties(ax_c2, 0);
offsetAxes(ax_c2);
saveFigureFormats(figC2, fullfile(outputDir, 'preferred_speed_bandpass_darkness'));

%% Figure E: Preferred speed, all speed-tuned boutons (Saleem et al., 2013, Fig. 2d style)
% Low-pass and high-pass grouped at the ends; band-pass at their fitted centre.
% Trough-shaped boutons excluded (their centre marks minimum activity).
nLow  = sum(strcmp(allTypes, 'lowpass'));
nHigh = sum(strcmp(allTypes, 'highpass'));
nTr   = sum(strcmp(allTypes, 'trough_inverted'));
countsE = [nLow, histBP, nHigh];
labelsE = [{'Low'}, arrayfun(@(a, b) sprintf('%g-%g', a, b), prefSpeedBins(1:end-1), prefSpeedBins(2:end), ...
           'UniformOutput', false), {'High'}];
barCols = [colors.lowpass; repmat(colors.bandpass, numel(histBP), 1); colors.highpass];

figE = figure('Position', [100, 100, 480, 360], 'Color', 'w');
ax_e = gca; hold on;
bE = bar(1:numel(countsE), countsE, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.7);
bE.CData = barCols;

ax_e.XTick              = 1:numel(labelsE);
ax_e.XTickLabel         = labelsE;
ax_e.XTickLabelRotation = 30;
ax_e.FontSize           = 8;
ax_e.FontName           = 'Arial';
ax_e.TickDir            = 'out';
ax_e.LineWidth          = 0.75;
ax_e.Layer              = 'top';
box off;

xlabel('Preferred running speed (cm s^{-1})', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Number of boutons', 'FontSize', 9, 'FontName', 'Arial');
title('RSP', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial');
subtitle(sprintf('n = %d speed-tuned boutons (trough-shaped excluded, n = %d)', nClassified - nTr, nTr), ...
    'FontSize', 8, 'Color', [0.45 0.45 0.45], 'FontName', 'Arial');
defaultAxesProperties(ax_e, 0);
offsetAxes(ax_e);
saveFigureFormats(figE, fullfile(outputDir, 'preferred_speed_histogram_darkness'));

%% Optional QC: PDF of all speed-tuned tuning curves, with Edd's fits
if makeQCpdf
    pdfSavePath = fullfile(saveFolder, 'exampleFits', 'all_speedTuned_tuningCurves_Edd.pdf');
    if isfile(pdfSavePath), delete(pdfSavePath); end
    typeColors = {colors.lowpass, colors.highpass, colors.bandpass, colors.trough_inverted};

    for s = 1:numel(sessionResults)
        E  = sessionResults(s).E;
        Y  = sessionResults(s).moveMean;
        nP = numel(E.sigIdx);
        if nP == 0, continue; end

        nB = numel(E.centres);
        rDense = linspace(1, nB, 100);
        xDense = interp1(1:nB, E.centres, rDense);

        nCols = 8; nRows = ceil(nP / nCols);
        fig = figure('Color', 'w', 'Position', [100 100 1600 max(nRows * 180, 300)], 'Visible', 'off');
        for k = 1:nP
            r = E.sigIdx(k); c = E.code(k);
            if isnan(c), continue; end
            p = E.params(k, :);
            ax_k = subplot(nRows, nCols, k); hold on;
            plot(E.centres, Y(r, :), 'o', 'Color', typeColors{c}, 'MarkerSize', 3, 'MarkerFaceColor', typeColors{c});
            plot(xDense, p(1) + p(2) * exp(-((rDense - p(3)).^2) / (2 * p(4)^2)), '-', 'Color', typeColors{c}, 'LineWidth', 1);
            title(sprintf('%s %d', codeNames{c}, r), 'FontSize', 5.5, 'Color', typeColors{c}, 'Interpreter', 'none');
            box off; set(ax_k, 'TickDir', 'out', 'FontSize', 5);
        end
        sgtitle(sprintf('%s  %s: %d / %d boutons speed-tuned', sessionResults(s).mouse, ...
            sessionResults(s).session, nP, E.nROIs), 'FontSize', 10, 'FontWeight', 'bold');
        exportgraphics(fig, pdfSavePath, 'Append', true, 'Resolution', 120);
        close(fig);
    end
    fprintf('QC PDF saved: %s\n', pdfSavePath);
end

%% ===================== Local functions =====================

function E = classifySessionEdd(tc, useField, alphaShuffle, customOffset)
% Selects speed-tuned boutons (circular-shift test on moving frames) and
% classifies their tuning shape with Edd's fitGaussianTemplates_tuning.
    Y       = tc.(useField).moveMean;             % nROIs x nBins
    pvals   = tc.(useField).pValMoving(:);
    edges   = tc.speedBins(:)';
    E.centres = edges(1:end-1) + diff(edges) / 2;
    E.nROIs   = size(Y, 1);
    E.sigIdx  = find(pvals <= alphaShuffle);

    n  = numel(E.sigIdx);
    nB = numel(E.centres);
    E.code   = nan(n, 1);
    E.pref   = nan(n, 1);
    E.R2     = nan(n, 1);
    E.params = nan(n, 4);

    for k = 1:n
        y = Y(E.sigIdx(k), :);
        if any(isnan(y)), continue; end
        [p, c, r2] = fitGaussianTemplates_tuning(num2cell(y(:)), customOffset, false);
        E.code(k) = c; E.R2(k) = r2; E.params(k, :) = p;
        if c == 3 || c == 4   % band-pass or trough: fitted centre, bin rank -> cm/s
            E.pref(k) = interp1(1:nB, E.centres, min(max(p(3), 1), nB));
        end
    end
end