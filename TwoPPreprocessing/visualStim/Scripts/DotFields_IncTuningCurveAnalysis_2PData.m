% DotFields_Pooled_AllTrials.m
%
% Pools DotMotion_SpeedTuning boutons across MULTIPLE mice/sessions, with
% NO running/stationary split -- all valid, non-blank trials used
% together regardless of behavioral state.
%
% Adapted from DotFields_TuningCurveAnalysis_compareStates_2PData.m's
% trial classification / calcium-mapping logic, generalized into a
% multi-session pooling loop (same pattern as the DirTuning scripts).
%
% SCOPE: ports trial classification + cross-validated R^2
% (calc_kfold_R2) + one-way ANOVA (speeds+blank) + 2SD-above-BLANK
% threshold (using genuine blank trials, not a pre-stim workaround) +
% Gaussian tuning-curve fitting + dynamic-range/Fano-factor + example
% fit overlays, adapted for the SINGLE-STATE case (no
% stationary-vs-running comparison, since there's only one state here).
%
% NOTE: unlike the DirTuning scripts (no blank stimulus available, so a
% pre-stim baseline window had to substitute), DotFields genuinely has
% blank trials (VelX1==1, numDots1==0). This script keeps them (the
% original script discarded them) and uses them as the literal blank
% reference for the ANOVA+threshold criterion -- a more faithful
% implementation of the originally published method.
%
% Requires calc_kfold_R2.m on MATLAB path.


mouseList = {'M25132', 'M25133', 'M26003'};  % 

stimWindowMask_range = [0 2.5];  % window for extracting per-trial calcium response (matches original script)
stimFramesMask_range = [0 2.0];  % window for wheel-speed behavior evaluation (matches original script; unused here since no state split, kept for reference)

ALPHA = 0.05;

r2opts.kval       = 3;
r2opts.nPerms     = 10;
r2opts.randFlag   = 1;
r2opts.validMeans = 1;
r2opts.nShuffle   = 100;
%% ===================================================

filteredTable = filterMasterTable('MouseID', mouseList, 'HasStimulus', 'DotMotion_SpeedTuning', ...
    'Suite2PPreprocessing', 1, 'Exclude', 0);
allMice    = filteredTable.MouseID;
uniqueMice = unique(allMice, 'stable');

allDotUnits = struct('alltraces', {}, 'blankTrials', {}, 'tuning', {}, 'sessionLabel', {});
haveWarnedShape = false;

for iMouse = 1:length(uniqueMice)
    thisMouse    = uniqueMice{iMouse};
    mouseSessIdx = find(strcmp(allMice, thisMouse));
    fprintf('MOUSE: %s | %d sessions\n', thisMouse, length(mouseSessIdx));

    for iSess = 1:length(mouseSessIdx)
        tableRow        = filteredTable(mouseSessIdx(iSess), :);
        thisSessionName = char(tableRow.Session);
        fprintf('  --- Session: %s ---\n', thisSessionName);

        infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
        if ~isfile(infoPath), warning('    sfi missing -- skipping.'); continue; end
        loadedInfo      = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames       = {sessionFileInfo.stimFiles.name};
        dotIdx = find(contains(stimNames, 'DotMotion_SpeedTuning'), 1);
        if isempty(dotIdx), warning('    No DotMotion_SpeedTuning file -- skipping.'); continue; end

        try
            load(sessionFileInfo.stimFiles(dotIdx).Response, 'response');
            load(sessionFileInfo.stimFiles(dotIdx).BonsaiData, 'bonsaiData');
        catch ME
            warning('    Could not load response/bonsaiData: %s', ME.message); continue;
        end
        if ~isfield(response, 'wheelData') || ~isfield(response, 'psthData')
            warning('    Missing wheelData or psthData -- skipping.'); continue;
        end

        %%build trial list (adapted from original script, no run/stat classification) ---
        nGroups = numel(response.wheelData);
        trialsSpeed2D = struct('VelX1', {}, 'numDots1', {}, 'startTime', {});
        trialCounter = 1;

        for g = 1:nGroups
            grpWheel  = response.wheelData(g);
            grpBonsai = bonsaiData.trialGroups(g);
            nTrialsInGrp = size(grpWheel.alignedResponses, 2);

            for ti = 1:nTrialsInGrp
                trialID = grpBonsai.trials(ti);
                trialsSpeed2D(trialCounter).VelX1 = grpWheel.stimValue;
                if grpWheel.stimValue == 1
                    trialsSpeed2D(trialCounter).numDots1 = 0; % blank
                else
                    trialsSpeed2D(trialCounter).numDots1 = 573;
                end
                trialsSpeed2D(trialCounter).startTime = bonsaiData.onARDTimes(trialID);
                trialsSpeed2D(trialCounter).origGroup = g;
                trialsSpeed2D(trialCounter).origTrialInGroup = ti;
                trialCounter = trialCounter + 1;
            end
        end

        tsd = trialsSpeed2D;
        temp_tsd  = tsd([tsd.numDots1] == 573); % real (non-blank) trials -- used for tuning/R2/gaussian fits
        blank_tsd = tsd([tsd.numDots1] == 0);   % genuine BLANK trials -- used as the baseline reference

        if isempty(temp_tsd)
            warning('    No non-blank trials for %s -- skipping.', thisSessionName);
            continue;
        end
        if isempty(blank_tsd)
            warning('    No blank trials found for %s -- ANOVA+threshold vs blank will be NaN for this session.', thisSessionName);
        end

        uniqueVelocities = unique(abs([temp_tsd.VelX1]));
        nSpeeds = numel(uniqueVelocities);

        nBoutons = size(response.psthData(1).alignedResponses, 1);
        if ~haveWarnedShape
            fprintf('    Detected shape: %d boutons x %d timepoints x %d trials (group 1)\n', ...
                size(response.psthData(1).alignedResponses,1), ...
                size(response.psthData(1).alignedResponses,2), ...
                size(response.psthData(1).alignedResponses,3));
            haveWarnedShape = true;
        end

        timeVec = response.psthData(1).timeVector(:)';
        stimWindowMask = timeVec >= stimWindowMask_range(1) & timeVec <= stimWindowMask_range(2);

        %%  per-bouton extraction (no state split), speeds AND blank 
        for thisROI = 1:nBoutons
            alltraces = cell(nSpeeds, 1);
            for s = 1:nSpeeds
                matchingTrials = find(abs([temp_tsd.VelX1]) == uniqueVelocities(s));
                traceAccumulator = nan(1, numel(matchingTrials));
                for mt = 1:numel(matchingTrials)
                    origGroup = temp_tsd(matchingTrials(mt)).origGroup;
                    origTi    = temp_tsd(matchingTrials(mt)).origTrialInGroup;
                    fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
                    traceAccumulator(mt) = nanmean(fullTrace(stimWindowMask));
                end
                alltraces{s} = traceAccumulator;
            end

            % genuine blank trials, same stimulus-window extraction
            blankAccumulator = nan(1, numel(blank_tsd));
            for mt = 1:numel(blank_tsd)
                origGroup = blank_tsd(mt).origGroup;
                origTi    = blank_tsd(mt).origTrialInGroup;
                fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
                blankAccumulator(mt) = nanmean(fullTrace(stimWindowMask));
            end

            allDotUnits(end+1) = struct( ...
                'alltraces',     {alltraces}, ...
                'blankTrials',   blankAccumulator, ...
                'tuning',        cellfun(@nanmean, alltraces), ...
                'sessionLabel',  sprintf('%s_%s', thisMouse, thisSessionName));
        end
        fprintf('    Added %d boutons (running total: %d).\n', nBoutons, numel(allDotUnits));
    end
end

nBoutonsTotal = numel(allDotUnits);
fprintf('\nDone pooling. %d boutons total across %d mice.\n', nBoutonsTotal, numel(uniqueMice));

%% ANOVA (speeds + blank) + 2SD-above-blank threshold (testing using the same metric as the rf mapping and direction tuning; this is not part of edd's script)
%  DotFields has genuine blank trials, so no
% pre-stim-window workaround is needed. ANOVA is run across all speed
% conditions PLUS the blank condition together (nSpeeds+1 groups), and
% the preferred speed's response is compared against the blank
% trials' mean + NSD*SD (not a pre-stim baseline period as done for the direction tuning stimulus).

alpha = 0.05;
numSD = 2;

anovaP_thresh    = nan(nBoutonsTotal, 1);
blankMean        = nan(nBoutonsTotal, 1);
blankSD          = nan(nBoutonsTotal, 1);
meanPreferredRaw = nan(nBoutonsTotal, 1);
prefSpeedIdx     = nan(nBoutonsTotal, 1);
isResponsive     = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    alltraces   = allDotUnits(b).alltraces;
    blankTrials = allDotUnits(b).blankTrials(:);
    blankTrials = blankTrials(~isnan(blankTrials));
    nSpeeds = numel(alltraces);

    % ANOVA across speeds AND blank (nSpeeds+1 groups)
    y = []; grp = [];
    for s = 1:nSpeeds
        vals = alltraces{s}(:);
        vals = vals(~isnan(vals));
        y   = [y; vals];
        grp = [grp; repmat(s, numel(vals), 1)];
    end
    y   = [y;   blankTrials];
    grp = [grp; repmat(nSpeeds + 1, numel(blankTrials), 1)]; % blank = its own group

    if numel(y) < (nSpeeds + 1) * 2 || isempty(blankTrials)
        continue; % not enough trials, or no blank trials for this bouton/session
    end
    anovaP_thresh(b) = anova1(y, grp, 'off');

    blankMean(b) = mean(blankTrials);
    blankSD(b)   = std(blankTrials);

    [~, prefSpeedIdx(b)] = max(allDotUnits(b).tuning);
    meanPreferredRaw(b) = mean(alltraces{prefSpeedIdx(b)}, 'omitnan');

    isResponsive(b) = (anovaP_thresh(b) < alpha) && ...
                      (meanPreferredRaw(b) > blankMean(b) + numSD * blankSD(b));
end

fprintf('\n%d / %d boutons classified as visually responsive (ANOVA[speeds+blank]+threshold, p<%.2f, >%dSD above BLANK).\n', ...
    sum(isResponsive), nBoutonsTotal, alpha, numSD);

for b = 1:nBoutonsTotal
    allDotUnits(b).anovaP_thresh    = anovaP_thresh(b);
    allDotUnits(b).blankMean        = blankMean(b);
    allDotUnits(b).blankSD          = blankSD(b);
    allDotUnits(b).meanPreferredRaw = meanPreferredRaw(b);
    allDotUnits(b).prefSpeedIdx     = prefSpeedIdx(b);
    allDotUnits(b).isResponsive     = isResponsive(b);
end

%%  t-test / Wilcoxon: preferred speed vs blank (protected by ANOVA) 
% Replaces the SD-multiplier heuristic with an ACTUAL statistical test
% comparing the preferred speed's single trials directly against the
% blank trials. ANOVA (already computed above) is kept as a "protected"
% gate: we only trust the direct preferred-vs-blank comparison if the
% ANOVA across ALL conditions was already significant. This
% guards against the selection-bias problem of picking whichever
% direction happens to look best and testing ONLY that one against
% blank with no correction -- picking the max out of nSpeeds conditions
% and then testing it in isolation is itself a form of multiple
% comparisons / "double dipping" if left unguarded.

ttestP    = nan(nBoutonsTotal, 1);
ranksumP  = nan(nBoutonsTotal, 1);
isResponsive_ttest   = false(nBoutonsTotal, 1);
isResponsive_ranksum = false(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    if isnan(anovaP_thresh(b)) || isnan(prefSpeedIdx(b))
        continue;
    end

    prefTrials  = allDotUnits(b).alltraces{prefSpeedIdx(b)}(:);
    prefTrials  = prefTrials(~isnan(prefTrials));
    blankTrials = allDotUnits(b).blankTrials(:);
    blankTrials = blankTrials(~isnan(blankTrials));

    if numel(prefTrials) < 2 || numel(blankTrials) < 2
        continue; % need at least 2 trials per group for a meaningful test
    end

    % Welch's t-test (does not assume equal variances between groups)
    [~, ttestP(b)] = ttest2(prefTrials, blankTrials, 'Vartype', 'unequal');

    % Wilcoxon rank-sum test (nonparametric, no normality assumption
    % more appropriate if dF/F trial distributions are skewed/heavy-tailed)
    ranksumP(b) = ranksum(prefTrials, blankTrials);

    %  only count as responsive if the  ANOVA
    % across all conditions was ALSO significant (guards against
    % selection bias from testing only the cherry-picked max speeds)
    isResponsive_ttest(b)   = (anovaP_thresh(b) < alpha) && (ttestP(b) < alpha) && ...
                              (meanPreferredRaw(b) > blankMean(b)); % direction check: must be ABOVE blank, not just different
    isResponsive_ranksum(b) = (anovaP_thresh(b) < alpha) && (ranksumP(b) < alpha) && ...
                              (meanPreferredRaw(b) > blankMean(b));
end

fprintf('\n%d / %d boutons responsive (ANOVA-protected Welch''s t-test, p<%.2f).\n', ...
    sum(isResponsive_ttest), nBoutonsTotal, alpha);
fprintf('%d / %d boutons responsive (ANOVA-protected Wilcoxon rank-sum, p<%.2f).\n', ...
    sum(isResponsive_ranksum), nBoutonsTotal, alpha);

for b = 1:nBoutonsTotal
    allDotUnits(b).ttestP               = ttestP(b);
    allDotUnits(b).ranksumP             = ranksumP(b);
    allDotUnits(b).isResponsive_ttest   = isResponsive_ttest(b);
    allDotUnits(b).isResponsive_ranksum = isResponsive_ranksum(b);
end

%% comparison: SD-heuristic vs t-test vs Wilcoxon vs cross-val R^2 
fprintf('\n Method comparison (n=%d boutons) \n', nBoutonsTotal);
fprintf('%-35s %6d\n', 'SD-heuristic (original)',        sum(isResponsive));
fprintf('%-35s %6d\n', 'ANOVA-protected Welch t-test',    sum(isResponsive_ttest));
fprintf('%-35s %6d\n', 'ANOVA-protected Wilcoxon',        sum(isResponsive_ranksum));

%% ANOVA + cross-validated R^2 
anovaP = nan(nBoutonsTotal, 1);
cvR2   = nan(nBoutonsTotal, 1);
cvPval = nan(nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    alltraces = allDotUnits(b).alltraces;
    nSpeeds = numel(alltraces);

    % ANOVA across speeds
    y = []; grp = [];
    for s = 1:nSpeeds
        vals = alltraces{s}(:);
        vals = vals(~isnan(vals));
        y   = [y; vals];
        grp = [grp; repmat(s, numel(vals), 1)];
    end
    if numel(y) >= nSpeeds * 2
        anovaP(b) = anova1(y, grp, 'off');
    end

    % cross-validated R^2 (downsample to equal trial count per speed first)
    trialCounts = cellfun(@(x) sum(~isnan(x)), alltraces);
    if all(trialCounts > 0)
        minTrial = min(trialCounts);
        if minTrial >= r2opts.kval
            gcaDownsampled = cellfun(@(x) x(~isnan(x)), alltraces, 'UniformOutput', false);
            gcaDownsampled = cellfun(@(x) x(1:minTrial), gcaDownsampled, 'UniformOutput', false);
            [cvR2(b), cvPval(b)] = calc_kfold_R2(gcaDownsampled, r2opts.kval, r2opts.nPerms, ...
                r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);
        end
    end

    if mod(b, 100) == 0
        fprintf('Processed %d / %d boutons...\n', b, nBoutonsTotal);
    end
end

isTunedCV = cvPval < ALPHA;
fprintf('\n%d / %d boutons show significant cross-validated speed tuning (p<%.2f).\n', ...
    sum(isTunedCV, 'omitnan'), nBoutonsTotal, ALPHA);

for b = 1:nBoutonsTotal
    allDotUnits(b).anovaP    = anovaP(b);
    allDotUnits(b).cvR2      = cvR2(b);
    allDotUnits(b).cvPval    = cvPval(b);
    allDotUnits(b).isTunedCV = isTunedCV(b);
end

%% population summary 
figure('Position', [100 100 900 400]);
subplot(1,2,1);
validR2 = ~isnan(cvR2);
histogram(cvR2(validR2), 30);
xlabel('Cross-validated R^2'); ylabel('Number of boutons');
title(sprintf('Distribution of tuning R^2 (n=%d)', sum(validR2)));

subplot(1,2,2);
cdfplot(cvR2(validR2));
xlabel('Cross-validated R^2'); ylabel('Proportion of boutons');
title('CDF of tuning R^2 (all trials, no state split)');

sgtitle('DotFields speed tuning: pooled across sessions/mice, no behavior split');

%%  Gaussian tuning-curve fits (single state) 
options.binSpacing = 2.0;

gaussChar   = nan(nBoutonsTotal, 1);
gaussR2fit  = nan(nBoutonsTotal, 1);
prefSpeed   = nan(nBoutonsTotal, 1);
gaussParamsAll = repmat({nan(1,4)}, nBoutonsTotal, 1);

for b = 1:nBoutonsTotal
    alltraces = allDotUnits(b).alltraces; % nSpeeds x 1
    tuning    = allDotUnits(b).tuning;

    % remove NaNs before fitting -- lsqcurvefit (inside
    % fitGaussianTemplates_tuning) requires finite YDATA, and a bad
    % trial (all-NaN trace) or empty cell will otherwise error out.
    alltracesClean = cellfun(@(x) x(~isnan(x)), alltraces, 'UniformOutput', false);
    if any(cellfun(@isempty, alltracesClean))
        continue; % at least one speed has zero valid trials -- leave NaN placeholders for this bouton
    end
    
    % this also calls edd's function
    [gp, gc, gr2] = fitGaussianTemplates_tuning(alltracesClean, 0.5, false);
    gaussParamsAll{b} = gp;
    gaussChar(b)  = gc;
    gaussR2fit(b) = gr2;

    if gc == 4
        [~, prefSpeed(b)] = min(tuning);
    else
        [~, prefSpeed(b)] = max(tuning);
    end

    if mod(b, 100) == 0
        fprintf('Gaussian fits: %d / %d boutons...\n', b, nBoutonsTotal);
    end
end

for b = 1:nBoutonsTotal
    allDotUnits(b).gaussParams = gaussParamsAll{b};
    allDotUnits(b).gaussChar   = gaussChar(b);
    allDotUnits(b).gaussR2fit  = gaussR2fit(b);
    allDotUnits(b).prefSpeed   = prefSpeed(b);
end

%%  dynamic range and Fano factor
for b = 1:nBoutonsTotal
    allDotUnits(b).dynamicRange = range(allDotUnits(b).tuning) / options.binSpacing;
    allDotUnits(b).fanoFactor   = mean(cellfun(@(x) var(x)/mean(x), allDotUnits(b).alltraces), 1);
end

%% population summaries: gaussChar, prefSpeed, dynamic range, Fano factor (edd's)
r2_thresh  = 0.1;
r2p_thresh = 0.05;
validIdx = find(cvR2 > r2_thresh & cvPval < r2p_thresh);
fprintf('\n%d / %d boutons pass R^2 filter (R^2>%.2f, p<%.2f) for descriptive summaries.\n', ...
    numel(validIdx), nBoutonsTotal, r2_thresh, r2p_thresh);

figure('Position', [100 100 1000 700]);
subplot(2,2,1);
histogram(gaussChar(validIdx), 'BinMethod', 'integers');
xlabel('Gaussian tuning category (gaussChar)'); ylabel('Count');
title('Tuning shape category');

subplot(2,2,2);
histogram(prefSpeed(validIdx), 'BinMethod', 'integers');
xlabel('Preferred speed index'); ylabel('Count');
title('Preferred speed distribution');

subplot(2,2,3);
histogram([allDotUnits(validIdx).dynamicRange]);
xlabel('Dynamic range'); ylabel('Count');
title('Dynamic range distribution');

subplot(2,2,4);
histogram([allDotUnits(validIdx).fanoFactor]);
xlabel('Fano factor'); ylabel('Count');
title('Fano factor distribution');

sgtitle(sprintf('Descriptive tuning summaries (pooled, filtered, n=%d)', numel(validIdx)));

%% individual Gaussian fit examples (first 6 cross-validated units) 
if ~isempty(validIdx)
    maxPlots = min(6, numel(validIdx));
    xDense = linspace(1, numel(allDotUnits(1).tuning), 100);

    figure('Color', 'w', 'Name', 'Cross-Validated Gaussian Fit Examples (pooled)', 'Position', [100, 50, 1100, 800]);
    for iPlot = 1:maxPlots
        b = validIdx(iPlot);
        subplot(2, 3, iPlot); hold on;

        nSpeedsHere = numel(allDotUnits(b).tuning);
        y_data = allDotUnits(b).tuning / options.binSpacing;

        p = allDotUnits(b).gaussParams;
        y_fit = p(1) + p(2) * exp(-((xDense - p(3)).^2) / (2 * p(4)^2));
        y_fit = y_fit / options.binSpacing;

        plot(1:nSpeedsHere, y_data, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5.5);
        plot(xDense, y_fit, 'k-', 'LineWidth', 2);

        allVals = [y_data(:); y_fit(:)];
        ylim([min(allVals) - 0.03, max(allVals) + 0.03]);
        xlim([0.5, nSpeedsHere + 0.5]); xticks(1:nSpeedsHere);
        ylabel('\DeltaF/F'); xlabel('Visual Speed Index');
        title(sprintf('Bouton %d (%s)\nType %d, R^2=%.2f', ...
            b, allDotUnits(b).sessionLabel, allDotUnits(b).gaussChar, allDotUnits(b).cvR2), ...
            'FontSize', 9, 'Interpreter', 'none');
    end
    sgtitle('Gaussian Model Fits (pooled, cross-validated cohort, no state split)', 'FontSize', 13, 'FontWeight', 'bold');
end

