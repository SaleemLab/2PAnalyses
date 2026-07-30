% this include data only locomotion trials pooled; 
% what speeds are the boutons
% run this to load all the data from another script; warning this lauches
% many figures; close all ofter inspection 
DotFields_IncTuningCurveAnalysis_RunningOnly_2PData
speeds = [0 16 32 64 128 256];
%% panel a b c d  - gaussian fits; include error bars; highlight with black dotted line the mean response during the blank

boutonsToPlot= [1 1 6 3]; %6
nPerCategory = 5;
nPerCategory_cat3 = 10; 
r2Thresh = 0.1;
cvPvalThresh = 0.05;

category = cat(1, allDotUnits.gaussChar_run);
r2       = cat(1, allDotUnits.runR2);
cvPval   = cat(1, allDotUnits.runR2_pval);
roiIdx   = cat(1, allDotUnits.roiIdx);
boutonIdx = (1:numel(allDotUnits))';
mouseID     = {allDotUnits.mouseID}';
sessionName = {allDotUnits.sessionName}';
sessionKey  = {allDotUnits.sessionLabel}'; % used only for diversity-checking

Category = []; BoutonIdx = []; MouseID = {}; Session = {}; RoiIdx = []; R2 = []; Pval = [];

for c = 1:4
    if c == 3
        thisNPerCategory = nPerCategory_cat3;
    else
        thisNPerCategory = nPerCategory;
    end

    idx = find(category == c & r2 >= r2Thresh & cvPval < cvPvalThresh);
    [~, order] = sort(r2(idx), 'descend');
    idx = idx(order);
    used = {};
    picked = [];
    for k = 1:numel(idx)
        key = sessionKey{idx(k)};
        if c ~= 3 && ismember(key, used), continue; end % diversity check skipped only for category 3
        picked(end+1) = idx(k); 
        used{end+1} = key;
        if numel(picked) >= thisNPerCategory, break; end
    end
    for k = 1:numel(picked)
        Category(end+1,1) = c; 
        BoutonIdx(end+1,1) = boutonIdx(picked(k)); 
        MouseID{end+1,1} = mouseID{picked(k)};
        Session{end+1,1} = sessionName{picked(k)}; 
        RoiIdx(end+1,1) = roiIdx(picked(k)); 
        R2(end+1,1) = r2(picked(k)); 
        Pval(end+1,1) = cvPval(picked(k)); 
    end
end

exampleTable = table(Category, BoutonIdx, MouseID, Session, RoiIdx, R2, Pval);
disp(exampleTable);

% can also use this for the previous figure to select boutons 

%bandpassRows = exampleTable(exampleTable.Category==1, :);
%plotDotFieldsExampleBouton(allDotUnits, bandpassRows.MouseID{1}, bandpassRows.Session{1}, bandpassRows.RoiIdx(1));
% plotDotFieldsExampleBouton(allDotUnits, bandpassRows.MouseID{1}, bandpassRows.Session{1}, bandpassRows.RoiIdx(1));

boutonsToPlot = [1 1 6 3]; 


figA = figure('Color', 'w', 'Position', [100 100 1000 260]);
for c = 1:4
    theseRows = exampleTable(exampleTable.Category == c, :);
    if isempty(theseRows), continue; end

    boutonIdx = theseRows.BoutonIdx(boutonsToPlot(c));
    unit = allDotUnits(boutonIdx);

    meanVals        = unit.tuning;
    semVals         = cellfun(@(x) std(x(~isnan(x)))/sqrt(sum(~isnan(x))), unit.alltraces);
    gaussParams     = unit.gaussParams_run;
    char            = unit.gaussChar_run;
    R2_thisUnit     = unit.runR2;
    pval            = unit.runR2_pval;
    fitR2           = unit.gaussR2_run;
    preStimBaseMean = unit.preStimBaselineMean;

    ax = subplot(1, 4, c);
    plotSpeedTuningGaussianFit(speeds, meanVals, semVals, gaussParams, char, R2_thisUnit, pval, fitR2, preStimBaseMean, ax);
end

baseFileName = 'visSpeedGaussianFits_4Types';
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section3_Fig4_4\gaussianfits';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = fullfile(outputDir, baseFileName);
saveFigureFormats(figA, fullSavePath);
%% Shared category names/colors (matches gaussChar_run: 1=Lowpass, 2=Highpass, 3=Bandpass, 4=Trough)
categoryNames  = {'Lowpass', 'Highpass', 'Bandpass', 'Trough'};
categoryColors = [ ...
    hex2rgb('#3288bd'); ...  % lowpass  - dark blue
    hex2rgb('#a6cee3'); ...  % highpass - light blue
    hex2rgb('#e6ab02'); ...  % bandpass - orange/gold
    hex2rgb('#66a61e')  ...  %  trough   - green
    ];

gaussChar_run = cat(1, allDotUnits.gaussChar_run);
prefSpeed_run = cat(1, allDotUnits.prefSpeed_run);
% validIdx_run should already exist from the main pipeline (R^2>0.1, p<0.05 filter)


%% panel e: probability histogram of tuning category
% for the four speed categories 
catCounts = histcounts(gaussChar_run(validIdx_run), 0.5:1:4.5);
catProbs  = catCounts / sum(catCounts);

figE=figure('Color', 'w');
b = bar(catProbs, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.6);
for i = 1:4
    b.CData(i,:) = categoryColors(i,:);
end
set(gca, 'XTick', 1:4, 'XTickLabel', categoryNames);
ylabel('Probability');
title(sprintf('Tuning class probability (n=%d)', numel(validIdx_run)));
box off;

hold on;
sessionLabels_valid = {allDotUnits(validIdx_run).sessionLabel}';
uniqueSessions = unique(sessionLabels_valid);

for i = 1:4
    perSessionProp = nan(numel(uniqueSessions), 1);
    for s = 1:numel(uniqueSessions)
        thisSessionMask = strcmp(sessionLabels_valid, uniqueSessions{s});
        nThisSession = sum(thisSessionMask);
        if nThisSession == 0, continue; end
        perSessionProp(s) = sum(gaussChar_run(validIdx_run(thisSessionMask)) == i) / nThisSession;
    end
    perSessionProp = perSessionProp(~isnan(perSessionProp));
    jitterX = i + (rand(size(perSessionProp)) - 0.5) * 0.3;
    scatter(jitterX, perSessionProp, 25, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);

    meanProp = mean(perSessionProp);
    semProp  = std(perSessionProp) / sqrt(numel(perSessionProp));
    errorbar(i, meanProp, semProp, 'k', 'LineWidth', 1.5, 'MarkerSize', 12, 'CapSize', 8);
end
hold off;
defaultAxesProperties(gca, true);
baseFileName = 'Tuning_class_probability';
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section3_Fig4_4\tuningCategories_gaussTemp';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = fullfile(outputDir, baseFileName);
saveFigureFormats(figE, fullSavePath);

%% PANEL f equivalent: excitatory vs suppressed, per category (relative to the pre-stimulus period)

% For each bouton, across its 6 conditions: run the paired signed-rank test per condition, 
% get a count of how many were excitatory and how many were suppressed (two numbers per bouton, each between 0 and 6).
% Loop through categories: for each of the 4 shape categories, find every bouton that belongs to it 
% and has already passed the reliability filter (validIdx_run).
% Take the mean: average those per-bouton excitatory counts across all boutons in that category
% (and separately for suppressed counts) — giving one number per category, which becomes the bar height.


% Bonferroni correction: dividing 0.05 by 6 because each bouton is tested
% at 6 separate speed conditions -- without this correction, running 6
% tests per bouton would inflate the false-positive rate well above 5%.

alpha_ExcSupp = 0.05 / 6;

% This loop computes exc/supp classification for EVERY bouton (not just
% validIdx_run) because respType_run is a per-bouton field stored
% once and reused later -- filtering to validIdx_run happens afterward,
% when we actually TALLY the counts, not here during computation.
for thisBouton = 1:numel(allDotUnits)
    alltraces     = allDotUnits(thisBouton).alltraces;
    preStimTrials = allDotUnits(thisBouton).preStimTrials;
    nSpeeds = numel(alltraces);
    respType = nan(nSpeeds, 1); % 1 = excitatory, -1 = suppressed, 0 = not significant

    for s = 1:nSpeeds
        evoked   = alltraces{s}(:);
        baseline = preStimTrials{s}(:);
        % signrank requires PAIRED, non-NaN samples -- this mask keeps
        % only trials where both the evoked value and its own matching
        % pre-stim baseline are valid numbers.
        validPairs = ~isnan(evoked) & ~isnan(baseline);
        evoked = evoked(validPairs);
        baseline = baseline(validPairs);

        % signrank needs a minimum number of pairs to be meaningful at
        % all; below this we simply leave respType as NaN for this speed.
        if numel(evoked) < 3
            continue;
        end

        % paired Wilcoxon signed-rank test: is the per-trial evoked value
        % reliably different from that SAME trial's pre-stim baseline?
        % This is what makes the test robust to each bouton's own
        % arbitrary fluorescence offset (see earlier discussion) -- the
        % offset appears on both sides of the pair and cancels out.
        p = signrank(evoked, baseline);
        if p < alpha_ExcSupp
            % direction of the effect is judged by the MEDIAN of the
            % paired differences, not the mean -- this matches what the
            % signed-rank test itself is actually testing, so the
            % direction call can't disagree with the significance call.
            % tried mean first and that was silly 
            if median(evoked - baseline) > 0
                respType(s) = 1;  % excitatory
            else
                respType(s) = -1; % suppressed
            end
        else
            respType(s) = 0; % not significant
        end
    end
    allDotUnits(thisBouton).respType_run = respType;
end

nExc  = zeros(4,1);
nSupp = zeros(4,1);
nBoutonsPerCat = zeros(4,1);

% NOW we restrict to validIdx_run -- only boutons that already passed the
% cross-validated R^2/p reliability filter contribute to this tally, so
% we're not counting exc/supp responses from boutons whose tuning curve
% itself isn't considered trustworthy in the first place.
for idx = validIdx_run(:)'
    c = allDotUnits(idx).gaussChar_run;
    if isnan(c) || c < 1 || c > 4, continue; end
    respType = allDotUnits(idx).respType_run;
    nExc(c)  = nExc(c)  + sum(respType == 1,  'omitnan');
    nSupp(c) = nSupp(c) + sum(respType == -1, 'omitnan');
    nBoutonsPerCat(c) = nBoutonsPerCat(c) + 1;
end

% mean number of significant responses PER BOUTON, matching Ed's y-axis
% scale (~0-2.5). We deliberately do NOT use raw summed totals (would
% just scale with how many boutons happen to be in each category, making
% categories with more boutons look artificially "more responsive") and
% NOT percentages (tested and looked wrong/uninterpretable against his
% figure) per-bouton mean count is what actually matches his panel.
meanExc  = nExc  ./ nBoutonsPerCat;
meanSupp = nSupp ./ nBoutonsPerCat;

figF=figure('Color', 'w');
bh = bar([meanExc, meanSupp], 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.6);
bh(1).FaceColor = [0.15 0.15 0.15]; % excitatory - dark
bh(2).FaceColor = [0.75 0.75 0.75]; % suppressed - light
set(gca, 'XTick', 1:4, 'XTickLabel', categoryNames);
legend({'Excitatory', 'Suppressed'}, 'Location', 'best');
ylabel('# significant responses (mean per bouton)');
title('Excitation vs suppression by tuning class');
box off;
defaultAxesProperties(gca, true);

baseFileName = 'exc_supp_per_tuningClass';
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section3_Fig4_4\tuningCategories_exc_supp';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = fullfile(outputDir, baseFileName);
saveFigureFormats(figF, fullSavePath);
%% PANEL g equivalent: preferred speed histogram, stacked by category
nSpeedsTotal = numel(speeds);
stackedCounts = zeros(nSpeedsTotal, 4);

% For each category, pull out only the boutons in BOTH validIdx_run AND
% that category, then count how many have each preferred-speed index.
% Building this as a [nSpeeds x 4] matrix (rather than 4 separate
% histograms) is what lets bar(...,'stacked') stack the 4 category
% colors on top of each other at each speed position, matching his
% stacked-bar style exactly.
for c = 1:4
    theseIdx = validIdx_run(gaussChar_run(validIdx_run) == c);
    theseSpeedIdx = prefSpeed_run(theseIdx);
    for s = 1:nSpeedsTotal
        stackedCounts(s, c) = sum(theseSpeedIdx == s);
    end
end

stackedProps = stackedCounts / numel(validIdx_run);

figG=figure('Color', 'w');
bh2 = bar(stackedProps, 'stacked', 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.6);
for c = 1:4
    bh2(c).CData = categoryColors(c,:);
end
set(gca, 'XTick', 1:nSpeedsTotal, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), speeds, 'UniformOutput', false));
xlabel('Preferred visual speed (\circ/s)');
ylabel('Proportion of boutons');
legend(categoryNames, 'Location', 'best');
title(sprintf('Preferred speed distribution (n=%d)', numel(validIdx_run)));
box off;
defaultAxesProperties(gca, true);

fullSavePath = fullfile(outputDir, baseFileName);
saveFigureFormats(figG, fullSavePath);
%% local helper
function rgb = hex2rgb(hexStr)
    hexStr = strrep(hexStr, '#', '');
    rgb = [hex2dec(hexStr(1:2)), hex2dec(hexStr(3:4)), hex2dec(hexStr(5:6))] / 255;
end