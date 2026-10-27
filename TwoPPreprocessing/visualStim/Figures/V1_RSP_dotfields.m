% RSP_vs_V1_SpeedTuning_Comparison.m
%
% Loads the running-only dot-field speed tuning data for RSP and V1 by
% running DotFields_IncTuningCurveAnalysis_RunningOnly_2PData once per area,
% then builds one comparison figure with metrics that are NOT in Fig 4.4:
%   a) proportion tuned per session
%   b) tuning class probability per area (session-level)
%   c) preferred speed CDF (continuous, from the Gaussian fit)
%   d) tuning width CDF (FWHM of the Gaussian, bandpass units)
%   e) mean peak-normalised tuning curve
%   f) fraction of units excited / suppressed at each speed
%   g) speed decoding accuracy vs population size (pseudo-population LDA)
%   h) signal correlation distributions (within RSP, within V1, across)
%
% No edit to the pipeline script is needed. loadArea() reads a copy of the
% pipeline, forces its mouseList to the list for that area, and runs the copy.
% Nothing is cached, so every run recomputes both areas (the cross-validated
% R^2 step is slow).
%
% Requires the Statistics and Machine Learning Toolbox (signrank, ranksum,
% fitcdiscr), plus your saveFigureFormats and defaultAxesProperties helpers.

clear; close all; clc
rng(1);

%% ------------------------------ config ---------------------------------
rspMice = {'M25132', 'M25133', 'M26003'};
v1Mice  = {'M25133', 'M26004', 'M26005', 'M25126'};

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section3_Fig4_4\RSP_vs_V1';

speeds      = [0 16 32 64 128 256];
areaNames   = {'RSP', 'V1'};
areaColors  = [hex2rgb('#7570b3'); hex2rgb('#d95f02')];
classNames  = {'Lowpass', 'Highpass', 'Bandpass', 'Trough'};

alphaBonf          = 0.05 / 6;   % same Bonferroni correction as the exc/supp panel
minUnitsPerSession = 3;          % sessions with fewer valid units are left out of session-level stats
minTrialsDecode    = 4;          % min trials per speed for a unit to enter the decoder
Ns                 = [5 10 20 50 100 200];
nRepDecode         = 30;
nPseudo            = 40;         % pseudo-trials per speed per split
maxCorrUnits       = 400;        % units per area used for signal correlations

%% ------------------------------ load -----------------------------------
A(1) = loadArea(rspMice, 'RSP', alphaBonf);
A(2) = loadArea(v1Mice,  'V1',  alphaBonf);

% The two areas must not be the same dataset.
if isequal(sort(unique(A(1).sessionLabel)), sort(unique(A(2).sessionLabel))) && A(1).n == A(2).n
    error('RSP and V1 contain identical sessions. Check the mouse lists and the master table.');
end

% M25133 is in both mouse lists. If it was filtered by mouse only, the same
% sessions will appear in both areas and contaminate the comparison.
overlap = intersect(unique(A(1).sessionLabel), unique(A(2).sessionLabel));
if ~isempty(overlap)
    warning('These sessions appear in BOTH areas (check the mouse lists / area assignment):');
    disp(overlap);
end

%% ------------------------- per-area metrics ----------------------------
xd = linspace(1, 6, 200);   % dense axis in speed-index units
R  = struct();

for a = 1:2
    S = A(a);
    v = S.valid;
    [us, ~, g] = unique(S.sessionLabel, 'stable');

    isValid = false(S.n, 1); isValid(v) = true;
    nU = accumarray(g, 1);
    nV = accumarray(g, double(isValid));

    R(a).nUnits    = S.n;
    R(a).nValid    = numel(v);
    R(a).nSessions = numel(unique(S.sessionLabel(v)));
    R(a).nMice     = numel(unique(S.mouseID(v)));
    R(a).fracTuned = nV ./ nU;

    % class probability per session (valid units only)
    propMat = nan(numel(us), 4);
    for s = 1:numel(us)
        idx = v(g(v) == s);
        if numel(idx) < minUnitsPerSession, continue; end
        propMat(s, :) = histcounts(S.gaussChar(idx), 0.5:1:4.5) / numel(idx);
    end
    R(a).propMat = propMat;

    % continuous preferred speed and width from the Gaussian fits
    prefCont = nan(numel(v), 1);
    fwhm     = nan(numel(v), 1);
    for k = 1:numel(v)
        i = v(k);
        p = S.gaussParams(i, :);
        if any(isnan(p)), continue; end
        y = gaussCurve(p, xd);
        if S.gaussChar(i) == 4, [~, j] = min(y); else, [~, j] = max(y); end
        prefCont(k) = xd(j);
        fwhm(k)     = 2 * sqrt(2 * log(2)) * abs(p(4));   % index units (~octaves above 16 deg/s)
    end
    R(a).prefCont   = prefCont;
    R(a).fwhm       = fwhm;
    R(a).validSess  = g(v);
    R(a).validClass = S.gaussChar(v);

    % peak-normalised tuning (min-max per unit)
    T  = S.tuning(v, :);
    Tn = (T - min(T, [], 2)) ./ (max(T, [], 2) - min(T, [], 2));
    R(a).tuning     = T;
    R(a).tuningNorm = Tn(all(~isnan(Tn), 2), :);

    % fraction excited / suppressed per speed
    rt = S.respType(v, :);
    R(a).fracExc  = sum(rt == 1,  1) ./ sum(~isnan(rt), 1);
    R(a).fracSupp = sum(rt == -1, 1) ./ sum(~isnan(rt), 1);
end

%% ---------------------------- decoding ---------------------------------
TRall = cell(1, 2);
for a = 1:2
    units = A(a).units;
    el = eligibleUnits(units, minTrialsDecode);
    idx = find(el);
    TR = cell(numel(idx), 1);
    for k = 1:numel(idx)
        at = units(idx(k)).alltraces;
        TR{k} = cellfun(@(x) x(~isnan(x)), at, 'UniformOutput', false);
    end
    TRall{a} = TR;
    fprintf('Decoding: %s has %d eligible units (>= %d trials at every speed).\n', ...
        areaNames{a}, numel(TR), minTrialsDecode);
end
Nmax   = min(cellfun(@numel, TRall));
NsUse  = Ns(Ns <= Nmax);
if isempty(NsUse) || NsUse(end) < Nmax, NsUse = [NsUse Nmax]; end
decAcc = nan(2, numel(NsUse), nRepDecode);
for a = 1:2
    for n = 1:numel(NsUse)
        for r = 1:nRepDecode
            decAcc(a, n, r) = decodeOnce(TRall{a}, NsUse(n), nPseudo);
        end
    end
    fprintf('Decoding done for %s.\n', areaNames{a});
end

%% ------------------------ signal correlations --------------------------
Tc = cell(1, 2);
for a = 1:2
    T = R(a).tuning;
    keep = all(~isnan(T), 2) & std(T, 0, 2) > 0;
    Tc{a} = T(keep, :);
end
nSub = min([cellfun(@(x) size(x, 1), Tc), maxCorrUnits]);
for a = 1:2
    Tc{a} = Tc{a}(randperm(size(Tc{a}, 1), nSub), :);
end
C11 = corr(Tc{1}'); C22 = corr(Tc{2}'); C12 = corr(Tc{1}', Tc{2}');
mask = triu(true(nSub), 1);
corrWithin = {C11(mask), C22(mask)};
corrAcross = C12(:);

%% ------------------------------ stats ----------------------------------
pFrac  = sessionTest(R(1).fracTuned, (1:numel(R(1).fracTuned))', R(2).fracTuned, (1:numel(R(2).fracTuned))', 1, false);
pPref  = sessionTest(R(1).prefCont, R(1).validSess, R(2).prefCont, R(2).validSess, minUnitsPerSession, true);
bp1 = R(1).validClass == 3; bp2 = R(2).validClass == 3;
pWidth = sessionTest(R(1).fwhm(bp1), R(1).validSess(bp1), R(2).fwhm(bp2), R(2).validSess(bp2), 2, true);

fprintf('\n=========== SUMMARY ===========\n');
for a = 1:2
    fprintf('%s: %d units, %d tuned (%.1f%%), %d sessions, %d mice\n', areaNames{a}, ...
        R(a).nUnits, R(a).nValid, 100 * R(a).nValid / R(a).nUnits, R(a).nSessions, R(a).nMice);
    fprintf('   median preferred speed index %.2f | median bandpass FWHM %.2f\n', ...
        median(R(a).prefCont, 'omitnan'), median(R(a).fwhm(R(a).validClass == 3), 'omitnan'));
end
fprintf('Session-level rank-sum p: fraction tuned %.3g | preferred speed %.3g | width %.3g\n', pFrac, pPref, pWidth);
fprintf('Signal correlation medians: RSP-RSP %.3f | V1-V1 %.3f | RSP-V1 %.3f (n = %d units per area)\n', ...
    median(corrWithin{1}, 'omitnan'), median(corrWithin{2}, 'omitnan'), median(corrAcross, 'omitnan'), nSub);
fprintf('Pairs are not independent, so no p-values are given for the correlation panel.\n');

%% ------------------------------ figure ---------------------------------
fig = figure('Color', 'w', 'Position', [50 50 1500 720]);
tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
letters = 'abcdefgh';
axAll = gobjects(1, 8);

% a) proportion tuned per session
ax = nexttile; axAll(1) = ax; hold on
for a = 1:2
    f = R(a).fracTuned;
    bar(a, mean(f), 'FaceColor', areaColors(a, :), 'EdgeColor', 'none', 'BarWidth', 0.6);
    scatter(a + (rand(size(f)) - 0.5) * 0.3, f, 22, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
    errorbar(a, mean(f), std(f) / sqrt(numel(f)), 'k', 'LineWidth', 1.5, 'CapSize', 8);
end
set(ax, 'XTick', 1:2, 'XTickLabel', areaNames, 'XLim', [0.4 2.6]);
ylabel('Proportion tuned (per session)');
title(sprintf('Tuned units (%s)', pStr(pFrac)));
box off

% b) class probability per area
ax = nexttile; axAll(2) = ax; hold on
w = 0.32;
for c = 1:4
    for a = 1:2
        xc = c + (a - 1.5) * (w + 0.04);
        pv = R(a).propMat(:, c); pv = pv(~isnan(pv));
        bar(xc, mean(pv), w, 'FaceColor', areaColors(a, :), 'EdgeColor', 'none');
        scatter(xc + (rand(size(pv)) - 0.5) * 0.12, pv, 14, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
        errorbar(xc, mean(pv), std(pv) / sqrt(numel(pv)), 'k', 'LineWidth', 1.2, 'CapSize', 5);
    end
end
set(ax, 'XTick', 1:4, 'XTickLabel', classNames);
ylabel('Probability (per session)');
title('Tuning class');
h1 = bar(nan, nan, 'FaceColor', areaColors(1, :), 'EdgeColor', 'none');
h2 = bar(nan, nan, 'FaceColor', areaColors(2, :), 'EdgeColor', 'none');
legend([h1 h2], areaNames, 'Location', 'northwest', 'Box', 'off');
box off

% c) preferred speed CDF
ax = nexttile; axAll(3) = ax; hold on
hl = gobjects(1, 2);
for a = 1:2
    hl(a) = plotCDF(R(a).prefCont, areaColors(a, :), '-');
    m = median(R(a).prefCont, 'omitnan');
    plot(m, 1.05, 'v', 'MarkerFaceColor', areaColors(a, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
end
set(ax, 'XTick', 1:6, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false), 'XLim', [1 6], ...
    'YLim', [0 1.1], 'YTick', 0:0.25:1);
xlabel('Preferred visual speed (\circ/s)'); ylabel('Cumulative fraction');
title(sprintf('Preferred speed (%s)', pStr(pPref)));
legend(hl, areaNames, 'Location', 'southeast', 'Box', 'off');
box off

% d) tuning width CDF (bandpass only)
ax = nexttile; axAll(4) = ax; hold on
hl = gobjects(1, 2);
for a = 1:2
    wv = R(a).fwhm(R(a).validClass == 3);
    hl(a) = plotCDF(wv, areaColors(a, :), '-');
    plot(median(wv, 'omitnan'), 1.05, 'v', 'MarkerFaceColor', areaColors(a, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
end
set(ax, 'YLim', [0 1.1], 'YTick', 0:0.25:1);
xlabel('Tuning width, FWHM (octaves)'); ylabel('Cumulative fraction');
title(sprintf('Bandpass width (%s)', pStr(pWidth)));
legend(hl, areaNames, 'Location', 'southeast', 'Box', 'off');
box off

% e) mean normalised tuning curve
ax = nexttile; axAll(5) = ax; hold on
for a = 1:2
    Tn = R(a).tuningNorm;
    m = mean(Tn, 1); s = std(Tn, 0, 1) / sqrt(size(Tn, 1));
    fill([1:6 6:-1:1], [m + s fliplr(m - s)], areaColors(a, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(1:6, m, '-o', 'Color', areaColors(a, :), 'LineWidth', 1.8, 'MarkerFaceColor', areaColors(a, :), 'MarkerSize', 5);
end
set(ax, 'XTick', 1:6, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false));
xlabel('Visual speed (\circ/s)'); ylabel('Normalised response');
title('Mean normalised tuning');
box off

% f) excited / suppressed per speed
ax = nexttile; axAll(6) = ax; hold on
for a = 1:2
    plot(1:6, R(a).fracExc,  '-o', 'Color', areaColors(a, :), 'LineWidth', 1.8, 'MarkerFaceColor', areaColors(a, :), 'MarkerSize', 5);
    plot(1:6, R(a).fracSupp, '--o', 'Color', areaColors(a, :), 'LineWidth', 1.8, 'MarkerFaceColor', 'w', 'MarkerSize', 5);
end
set(ax, 'XTick', 1:6, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false));
xlabel('Visual speed (\circ/s)'); ylabel('Fraction of tuned units');
title('Excited (solid) / suppressed (dashed)');
box off

% g) decoding accuracy vs population size
ax = nexttile; axAll(7) = ax; hold on
for a = 1:2
    d = squeeze(decAcc(a, :, :));       % nN x nRep
    if isvector(d), d = d(:); end
    m = mean(d, 2)'; s = std(d, 0, 2)';
    fill([NsUse fliplr(NsUse)], [m + s fliplr(m - s)], areaColors(a, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(NsUse, m, '-o', 'Color', areaColors(a, :), 'LineWidth', 1.8, 'MarkerFaceColor', areaColors(a, :), 'MarkerSize', 5);
end
yline(1/6, 'k:', 'LineWidth', 1.2);
set(ax, 'XScale', 'log', 'XTick', NsUse, 'XTickLabel', arrayfun(@num2str, NsUse, 'UniformOutput', false));
xlabel('Number of units'); ylabel('Decoding accuracy');
title('Speed decoding (6 classes)');
ylim([0 1]);
box off

% h) signal correlations
ax = nexttile; axAll(8) = ax; hold on
hl = gobjects(1, 3);
hl(1) = plotCDF(corrWithin{1}, areaColors(1, :), '-');
hl(2) = plotCDF(corrWithin{2}, areaColors(2, :), '-');
hl(3) = plotCDF(corrAcross, [0.4 0.4 0.4], '--');
xlabel('Signal correlation'); ylabel('Cumulative fraction');
title(sprintf('Tuning similarity (n = %d units/area)', nSub));
legend(hl, {'RSP-RSP', 'V1-V1', 'RSP-V1'}, 'Location', 'southeast', 'Box', 'off');
box off

% panel letters
for k = 1:8
    text(axAll(k), -0.18, 1.12, letters(k), 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
end

if exist('defaultAxesProperties', 'file')
    for k = 1:8, defaultAxesProperties(axAll(k), true); end
end

if ~exist(outputDir, 'dir'), mkdir(outputDir); end
if exist('saveFigureFormats', 'file')
    saveFigureFormats(fig, fullfile(outputDir, 'RSP_vs_V1_speedTuning_comparison'));
else
    exportgraphics(fig, fullfile(outputDir, 'RSP_vs_V1_speedTuning_comparison.pdf'), 'ContentType', 'vector');
end

%% ============================ local functions ===========================
function S = loadArea(mice, name, alphaBonf)
    % Runs a copy of the pipeline script in this function's workspace, with
    % its mouseList forced to `mice`, so the original file is never edited.
    areaLabel__ = name;
    alphaBonf__ = alphaBonf;

    pipelineName = 'DotFields_IncTuningCurveAnalysis_RunningOnly_2PData';
    srcFile = which(pipelineName);
    if isempty(srcFile)
        error('Cannot find %s on the MATLAB path.', pipelineName);
    end
    srcLines = regexp(fileread(srcFile), '\r?\n', 'split');
    isML = ~cellfun('isempty', regexp(srcLines, '^\s*mouseList\s*=', 'once'));
    if ~any(isML)
        error('No active "mouseList = ..." line found in %s.', srcFile);
    end
    srcLines(isML) = {'mouseList = mice;'};

    tmpName = ['tmpPipeline_' areaLabel__];
    tmpFile = fullfile(tempdir, [tmpName '.m']);
    fid = fopen(tmpFile, 'w');
    fwrite(fid, strjoin(srcLines, newline));
    fclose(fid);
    addpath(tempdir); rehash;
    cleaner = onCleanup(@() cleanupTmp(tmpFile)); %#ok<NASGU>

    fprintf('\n##### Loading %s (mice: %s) #####\n', areaLabel__, strjoin(mice, ', '));
    eval(tmpName);   % runs the patched pipeline in this workspace

    S = packArea(allDotUnits, validIdx_run, areaLabel__, alphaBonf__); %#ok<NODEF>
    loaded = unique({allDotUnits.mouseID});
    fprintf('%s loaded with mice: %s\n', areaLabel__, strjoin(loaded, ', '));
    missing = setdiff(mice, loaded);
    if ~isempty(missing)
        warning('%s: no data found for mice %s (check the master table filters).', areaLabel__, strjoin(missing, ', '));
    end
    close all
end

function cleanupTmp(tmpFile)
    if isfile(tmpFile), delete(tmpFile); end
    rmpath(tempdir);
end

function S = packArea(units, validIdx, name, alphaBonf)
    n = numel(units);
    S.name         = name;
    S.n            = n;
    S.units        = units;
    S.sessionLabel = {units.sessionLabel}';
    S.mouseID      = {units.mouseID}';
    S.gaussChar    = cat(1, units.gaussChar_run);
    S.gaussParams  = nan(n, 4);
    S.tuning       = nan(n, 6);
    S.respType     = nan(n, 6);
    S.has6         = false(n, 1);

    for i = 1:n
        at = units(i).alltraces;
        if numel(at) ~= 6, continue; end
        S.has6(i) = true;
        S.tuning(i, :) = units(i).tuning(:)';
        gp = units(i).gaussParams_run;
        if numel(gp) == 4, S.gaussParams(i, :) = gp(:)'; end

        pre = units(i).preStimTrials;
        for s = 1:6
            ev = at{s}(:); bl = pre{s}(:);
            ok = ~isnan(ev) & ~isnan(bl);
            ev = ev(ok); bl = bl(ok);
            if numel(ev) < 3, continue; end
            if signrank(ev, bl) < alphaBonf
                if median(ev - bl) > 0, S.respType(i, s) = 1; else, S.respType(i, s) = -1; end
            else
                S.respType(i, s) = 0;
            end
        end
    end

    validIdx = validIdx(:);
    S.valid = validIdx(S.has6(validIdx));
end

function y = gaussCurve(p, x)
    y = p(1) + p(2) * exp(-((x - p(3)).^2) / (2 * p(4)^2));
end

function el = eligibleUnits(units, minTr)
    el = false(numel(units), 1);
    for i = 1:numel(units)
        at = units(i).alltraces;
        if numel(at) ~= 6, continue; end
        n = cellfun(@(x) sum(~isnan(x)), at);
        el(i) = all(n >= minTr);
    end
end

function acc = decodeOnce(TR, N, nPseudo)
    % Pseudo-population LDA. Each unit's trials are split into train/test
    % halves, pseudo-trials are resampled independently per unit (so noise
    % correlations are ignored), and accuracy is measured on the test half.
    nS  = 6;
    sel = randperm(numel(TR), N);
    Xtr = nan(nS * nPseudo, N);
    Xte = nan(nS * nPseudo, N);
    y   = repelem((1:nS)', nPseudo);
    for u = 1:N
        for s = 1:nS
            t = TR{sel(u)}{s};
            t = t(randperm(numel(t)));
            h = floor(numel(t) / 2);
            trn = t(1:h); tst = t(h+1:end);
            rows = (s - 1) * nPseudo + (1:nPseudo);
            Xtr(rows, u) = trn(randi(numel(trn), nPseudo, 1));
            Xte(rows, u) = tst(randi(numel(tst), nPseudo, 1));
        end
    end
    mu = mean(Xtr, 1); sd = std(Xtr, 0, 1); sd(sd == 0) = 1;
    Xtr = (Xtr - mu) ./ sd;
    Xte = (Xte - mu) ./ sd;
    mdl = fitcdiscr(Xtr, y, 'DiscrimType', 'pseudoLinear');
    acc = mean(predict(mdl, Xte) == y);
end

function [p, mA, mB] = sessionTest(valsA, sessA, valsB, sessB, minN, doMedian)
    % Collapse to one value per session, then rank-sum across sessions.
    mA = sessionSummary(valsA, sessA, minN, doMedian);
    mB = sessionSummary(valsB, sessB, minN, doMedian);
    if numel(mA) >= 2 && numel(mB) >= 2
        p = ranksum(mA, mB);
    else
        p = nan;
    end
end

function m = sessionSummary(vals, sess, minN, doMedian)
    vals = vals(:); sess = sess(:);
    ok = ~isnan(vals); vals = vals(ok); sess = sess(ok);
    if isempty(vals), m = []; return; end
    [~, ~, g] = unique(sess);
    n = accumarray(g, 1);
    if doMedian
        m = accumarray(g, vals, [], @median);
    else
        m = accumarray(g, vals, [], @mean);
    end
    m(n < minN) = [];
end

function h = plotCDF(v, col, ls)
    v = sort(v(~isnan(v)));
    if isempty(v), h = plot(nan, nan, ls, 'Color', col, 'LineWidth', 1.8); return; end
    f = (1:numel(v)) / numel(v);
    h = plot(v, f, ls, 'Color', col, 'LineWidth', 1.8);
end

function s = pStr(p)
    if isnan(p),      s = 'n/a';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'ns';
    end
end

function rgb = hex2rgb(hexStr)
    hexStr = strrep(hexStr, '#', '');
    rgb = [hex2dec(hexStr(1:2)), hex2dec(hexStr(3:4)), hex2dec(hexStr(5:6))] / 255;
end