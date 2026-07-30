% DotFields_TuningCurveAnalysis_pooledStates_2PData.m
%
% POOLED-STATE companion to
% DotFields_TuningCurveAnalysis_compareStatesV2_2PData.m
%
%   - Same wheel-based trial classification as the state-split script
%     (ambiguous trials are still excluded), BUT stationary and running
%     trials are pooled together per speed/blank condition rather than
%     kept separate. This reproduces the "pooled across states" numbers
%     used to motivate the state-split analysis (e.g. 550/2313 boutons
%     responsive when pooled vs. 261/2313 stationary-only and 962/2313
%     running-only).
%   - Same ANOVA(speeds+blank)+2SD-threshold, ANOVA-protected t-test /
%     Wilcoxon, cross-validated R^2, Gaussian fitting, dynamic
%     range/Fano factor, computed ONCE per bouton on the pooled trial
%     set (no per-state split).
%
% Requires calc_kfold_R2.m and fitGaussianTemplates_tuning.m from
% NeuralDataAnalyses (github repo)
%
% Use this alongside (not instead of) the state-split script -- this one
% answers "what does tuning look like if we don't consider behavioral
% state," which is the justification step before the state-split
% analysis, not a replacement for it.

%%
mouseList = {'M25132', 'M25133', 'M26003'};

stimWindowMask_range = [0.5 2];  % window for extracting per-trial calcium response
stimFramesMask_range = [0 2.0];  % window for wheel-speed behavior classification

runSpeedThresh  = 3;    % mean speed above this =  "running"
statSpeedThresh = 0.5;  % mean speed below this = candidate "stationary"
propThresh      = 0.75;

ALPHA = 0.05;
NSD   = 1;

r2opts.kval       = 3;
r2opts.nPerms     = 10;
r2opts.randFlag   = 1;
r2opts.validMeans = 1;
r2opts.nShuffle   = 100;

options.binSpacing = 2.0; % matches 2s-on/2s-off stimulus timing

%%
filteredTable = filterMasterTable('MouseID', mouseList, 'HasStimulus', 'DotMotion_SpeedTuning', ...
    'Suite2PPreprocessing', 1, 'Exclude', 0);
allMice    = filteredTable.MouseID;
uniqueMice = unique(allMice, 'stable');

allDotUnits = struct('alltraces', {}, 'blankTrials', {}, 'tuning', {}, 'sessionLabel', {}, ...
    'mouseID', {}, 'sessionName', {}, 'roiIdx', {});
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

        %%  build trial list + classify running/stationary for EVERY trial (incl. blank)
        % NOTE: classification is still computed (and ambiguous trials
        % still excluded) so the pooled analysis uses the SAME trial
        % inclusion criteria as the state-split script -- the only
        % difference is stat+run trials are pooled together below,
        % rather than kept in separate columns.
        nGroups = numel(response.wheelData);
        trialsSpeed2D = struct('VelX1', {}, 'numDots1', {}, 'runFlag', {}, ...
                                'origGroup', {}, 'origTrialInGroup', {});
        trialCounter = 1;

        for g = 1:nGroups
            grpWheel  = response.wheelData(g);
            grpBonsai = bonsaiData.trialGroups(g);

            speedMatrix = grpWheel.alignedResponses;
            tVec        = grpWheel.timeVector;
            stimFramesMask = (tVec >= stimFramesMask_range(1) & tVec <= stimFramesMask_range(2));

            for ti = 1:size(speedMatrix, 2)
                singleTrialTrace = speedMatrix(:, ti);
                if all(isnan(singleTrialTrace)), continue; end

                meanSpeed      = nanmean(singleTrialTrace(stimFramesMask));
                propRunning    = sum(singleTrialTrace(stimFramesMask) > statSpeedThresh) / sum(stimFramesMask);
                propStationary = sum(singleTrialTrace(stimFramesMask) < runSpeedThresh)  / sum(stimFramesMask);

                runFlag = NaN;
                if propRunning >= propThresh && meanSpeed > runSpeedThresh
                    runFlag = 1;
                elseif propStationary >= propThresh && meanSpeed < statSpeedThresh
                    runFlag = 0;
                end
                if isnan(runFlag), continue; end % ambiguous/ITI-crossing trial -- excluded, not guessed

                trialsSpeed2D(trialCounter).VelX1 = grpWheel.stimValue;
                if grpWheel.stimValue == 1
                    trialsSpeed2D(trialCounter).numDots1 = 0; % blank
                else
                    trialsSpeed2D(trialCounter).numDots1 = 573; % all other speeds
                end
                trialsSpeed2D(trialCounter).runFlag          = runFlag; % kept for bookkeeping only -- NOT used to split below
                trialsSpeed2D(trialCounter).origGroup        = g;
                trialsSpeed2D(trialCounter).origTrialInGroup = ti;

                trialCounter = trialCounter + 1;
            end
        end

        if trialCounter == 1
            warning('    No valid classified trials for %s -- skipping.', thisSessionName);
            continue;
        end

        tsd = trialsSpeed2D;
        temp_tsd  = tsd([tsd.numDots1] == 573); % real (non-blank) trials
        blank_tsd = tsd([tsd.numDots1] == 0);   % genuine BLANK trials

        if isempty(temp_tsd)
            warning('    No non-blank trials for %s -- skipping.', thisSessionName);
            continue;
        end
        if isempty(blank_tsd)
            warning('    No blank trials for %s -- blank-dependent stats will be NaN.', thisSessionName);
        end

        uniqueVelocities = unique(abs([temp_tsd.VelX1]));
        nSpeeds = numel(uniqueVelocities); % 0    16    32    64   128   256 (speed 0 is also included here)

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

        %%  per-bouton, POOLED responses (speeds AND blank; stat+run combined)
        for thisROI = 1:nBoutons
            alltraces   = cell(nSpeeds, 1); % single column -- pooled across state
            blankTrials = cell(1, 1);

            for thisSpeed = 1:nSpeeds
                % NOTE: no runFlag filter here -- this is the only
                % structural difference from the state-split script.
                matchingTrials = find(abs([temp_tsd.VelX1]) == uniqueVelocities(thisSpeed));
                traceAccumulator = nan(1, numel(matchingTrials));
                for mt = 1:numel(matchingTrials)
                    origGroup = temp_tsd(matchingTrials(mt)).origGroup;
                    origTi    = temp_tsd(matchingTrials(mt)).origTrialInGroup;
                    fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
                    traceAccumulator(mt) = nanmean(fullTrace(stimWindowMask));
                end
                alltraces{thisSpeed, 1} = traceAccumulator;
            end

            blankAccum = nan(1, numel(blank_tsd));
            for mt = 1:numel(blank_tsd)
                origGroup = blank_tsd(mt).origGroup;
                origTi    = blank_tsd(mt).origTrialInGroup;
                fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
                blankAccum(mt) = nanmean(fullTrace(stimWindowMask));
            end
            blankTrials{1} = blankAccum;

            allDotUnits(end+1) = struct( ...
                'alltraces',    {alltraces}, ...
                'blankTrials',  {blankTrials}, ...
                'tuning',       cellfun(@nanmean, alltraces), ...
                'sessionLabel', sprintf('%s_%s', thisMouse, thisSessionName), ...
                'mouseID',      thisMouse, ...
                'sessionName',  thisSessionName, ...
                'roiIdx',       thisROI);
        end
        fprintf('    Added %d boutons (running total: %d).\n', nBoutons, numel(allDotUnits));
    end
end

nBoutonsTotal = numel(allDotUnits);
fprintf('\nDone pooling. %d boutons total across %d mice.\n', nBoutonsTotal, numel(uniqueMice));

%%  ANOVA(speeds+blank) + 2SD-above-blank threshold, POOLED
anovaP_thresh    = nan(nBoutonsTotal, 1);
blankMean        = nan(nBoutonsTotal, 1);
blankSD          = nan(nBoutonsTotal, 1);
meanPreferredRaw = nan(nBoutonsTotal, 1);
prefSpeedIdx     = nan(nBoutonsTotal, 1);
isResponsive     = false(nBoutonsTotal, 1);

for thisBouton = 1:nBoutonsTotal
    alltraces   = allDotUnits(thisBouton).alltraces(:, 1);
    blankTrials = allDotUnits(thisBouton).blankTrials{1}(:);
    blankTrials = blankTrials(~isnan(blankTrials));
    nSpeeds = numel(alltraces);

    y = [];
    grp = [];
    for thisSpeed = 1:nSpeeds
        vals = alltraces{thisSpeed}(:);
        vals = vals(~isnan(vals));
        y   = [y; vals];
        grp = [grp; repmat(thisSpeed, numel(vals), 1)];
    end
    y   = [y;   blankTrials];
    grp = [grp; repmat(nSpeeds + 1, numel(blankTrials), 1)];

    if numel(y) < (nSpeeds + 1) * 2 || isempty(blankTrials)
        continue;
    end
    anovaP_thresh(thisBouton) = anova1(y, grp, 'off');
    blankMean(thisBouton) = mean(blankTrials);
    blankSD(thisBouton)   = std(blankTrials);

    tuningThisState = cellfun(@(x) mean(x, 'omitnan'), alltraces);
    [~, prefSpeedIdx(thisBouton)] = max(tuningThisState);
    meanPreferredRaw(thisBouton) = mean(alltraces{prefSpeedIdx(thisBouton)}, 'omitnan');

    isResponsive(thisBouton) = (anovaP_thresh(thisBouton) < ALPHA) && ...
                      (meanPreferredRaw(thisBouton) > blankMean(thisBouton) + NSD * blankSD(thisBouton));
end

fprintf('\n[Pooled] %d / %d boutons responsive (ANOVA[speeds+blank]+threshold, p<%.2f, >%dSD above blank).\n', ...
    sum(isResponsive), nBoutonsTotal, ALPHA, NSD);

for thisBouton = 1:nBoutonsTotal
    allDotUnits(thisBouton).anovaP_thresh_pooled    = anovaP_thresh(thisBouton);
    allDotUnits(thisBouton).blankMean_pooled        = blankMean(thisBouton);
    allDotUnits(thisBouton).blankSD_pooled          = blankSD(thisBouton);
    allDotUnits(thisBouton).meanPreferredRaw_pooled = meanPreferredRaw(thisBouton);
    allDotUnits(thisBouton).prefSpeedIdx_pooled     = prefSpeedIdx(thisBouton);
    allDotUnits(thisBouton).isResponsive_pooled     = isResponsive(thisBouton);
end

%%  ANOVA-protected t-test / Wilcoxon vs blank, POOLED
ttestP    = nan(nBoutonsTotal, 1);
ranksumP  = nan(nBoutonsTotal, 1);
isResp_ttest   = false(nBoutonsTotal, 1);
isResp_ranksum = false(nBoutonsTotal, 1);

for thisBouton = 1:nBoutonsTotal
    anovaP_thresh_b = allDotUnits(thisBouton).anovaP_thresh_pooled;
    prefSpeedIdx_b  = allDotUnits(thisBouton).prefSpeedIdx_pooled;
    if isnan(anovaP_thresh_b) || isnan(prefSpeedIdx_b)
        continue;
    end

    prefTrials  = allDotUnits(thisBouton).alltraces{prefSpeedIdx_b, 1}(:);
    prefTrials  = prefTrials(~isnan(prefTrials));
    blankTrials = allDotUnits(thisBouton).blankTrials{1}(:);
    blankTrials = blankTrials(~isnan(blankTrials));

    if numel(prefTrials) < 2 || numel(blankTrials) < 2
        continue;
    end

    [~, ttestP(thisBouton)] = ttest2(prefTrials, blankTrials, 'Vartype', 'unequal');
    ranksumP(thisBouton) = ranksum(prefTrials, blankTrials);

    blankMean_b = allDotUnits(thisBouton).blankMean_pooled;
    meanPreferredRaw_b = allDotUnits(thisBouton).meanPreferredRaw_pooled;

    isResp_ttest(thisBouton)   = (anovaP_thresh_b < ALPHA) && (ttestP(thisBouton) < ALPHA) && (meanPreferredRaw_b > blankMean_b);
    isResp_ranksum(thisBouton) = (anovaP_thresh_b < ALPHA) && (ranksumP(thisBouton) < ALPHA) && (meanPreferredRaw_b > blankMean_b);
end

fprintf('[Pooled] %d / %d boutons responsive (ANOVA-protected Welch t-test, p<%.2f).\n', ...
    sum(isResp_ttest), nBoutonsTotal, ALPHA);
fprintf('[Pooled] %d / %d boutons responsive (ANOVA-protected Wilcoxon, p<%.2f).\n', ...
    sum(isResp_ranksum), nBoutonsTotal, ALPHA);

for thisBouton = 1:nBoutonsTotal
    allDotUnits(thisBouton).ttestP_pooled             = ttestP(thisBouton);
    allDotUnits(thisBouton).ranksumP_pooled           = ranksumP(thisBouton);
    allDotUnits(thisBouton).isResponsive_ttest_pooled   = isResp_ttest(thisBouton);
    allDotUnits(thisBouton).isResponsive_ranksum_pooled = isResp_ranksum(thisBouton);
end

%% method comparison summary, POOLED
fprintf('\n=== Method comparison [Pooled] (n=%d boutons) ===\n', nBoutonsTotal);
fprintf('%-35s %6d\n', 'SD-heuristic',                 sum([allDotUnits.isResponsive_pooled]));
fprintf('%-35s %6d\n', 'ANOVA-protected Welch t-test',  sum([allDotUnits.isResponsive_ttest_pooled]));
fprintf('%-35s %6d\n', 'ANOVA-protected Wilcoxon',      sum([allDotUnits.isResponsive_ranksum_pooled]));

%%  cross-validated R^2, POOLED
statR2_pooled = nan(nBoutonsTotal, 1); statR2_pooled_pval = nan(nBoutonsTotal, 1);

for thisBouton = 1:nBoutonsTotal
    alltraces = allDotUnits(thisBouton).alltraces; % nSpeeds x 1

    minTrial = min(cellfun(@(x) sum(~isnan(x)), alltraces));
    if minTrial < r2opts.kval
        continue;
    end
    alltracesDownsampled = cellfun(@(x) x(find(~isnan(x), minTrial)), alltraces, 'UniformOutput', false);

    gca_pooled = alltracesDownsampled(:, 1)';
    [statR2_pooled(thisBouton), statR2_pooled_pval(thisBouton)] = calc_kfold_R2(gca_pooled, r2opts.kval, r2opts.nPerms, ...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle);

    if mod(thisBouton, 100) == 0
        fprintf('Cross-val R^2: %d / %d boutons...\n', thisBouton, nBoutonsTotal);
    end
end

isTunedPooled = statR2_pooled_pval < ALPHA;
fprintf('\n[Pooled] %d / %d boutons show significant cross-validated tuning (p<%.2f).\n', ...
    sum(isTunedPooled, 'omitnan'), nBoutonsTotal, ALPHA);

for thisBouton = 1:nBoutonsTotal
    allDotUnits(thisBouton).R2_pooled      = statR2_pooled(thisBouton);
    allDotUnits(thisBouton).R2_pooled_pval = statR2_pooled_pval(thisBouton);
    allDotUnits(thisBouton).isTuned_pooled = isTunedPooled(thisBouton);
end

%%  Gaussian tuning-curve fits, POOLED
gaussParams_pooled_all = repmat({nan(1,4)}, nBoutonsTotal, 1);
gaussChar_pooled_all   = nan(nBoutonsTotal, 1);
gaussR2_pooled_all     = nan(nBoutonsTotal, 1);
prefSpeed_pooled_all   = nan(nBoutonsTotal, 1);

for thisBouton = 1:nBoutonsTotal
    alltraces = allDotUnits(thisBouton).alltraces;
    tuning    = allDotUnits(thisBouton).tuning;

    alltracesClean = cellfun(@(x) x(~isnan(x)), alltraces, 'UniformOutput', false);
    if any(cellfun(@isempty, alltracesClean(:)))
        continue;
    end

    [gp_pooled, gc_pooled, gr2_pooled] = fitGaussianTemplates_tuning(alltracesClean(:,1), 0.5, false);
    if gc_pooled == 4, [~, ps_pooled] = min(tuning(:,1)); else, [~, ps_pooled] = max(tuning(:,1)); end

    gaussParams_pooled_all{thisBouton} = gp_pooled;
    gaussChar_pooled_all(thisBouton)   = gc_pooled;
    gaussR2_pooled_all(thisBouton)     = gr2_pooled;
    prefSpeed_pooled_all(thisBouton)   = ps_pooled;

    if mod(thisBouton, 100) == 0
        fprintf('Gaussian fits: %d / %d boutons...\n', thisBouton, nBoutonsTotal);
    end
end

for thisBouton = 1:nBoutonsTotal
    allDotUnits(thisBouton).gaussParams_pooled = gaussParams_pooled_all{thisBouton};
    allDotUnits(thisBouton).gaussChar_pooled   = gaussChar_pooled_all(thisBouton);
    allDotUnits(thisBouton).gaussR2_pooled     = gaussR2_pooled_all(thisBouton);
    allDotUnits(thisBouton).prefSpeed_pooled   = prefSpeed_pooled_all(thisBouton);
end

%% dynamic range and Fano factor, POOLED
for thisBouton = 1:nBoutonsTotal
    alltraces = allDotUnits(thisBouton).alltraces;
    tuning    = allDotUnits(thisBouton).tuning;

    allDotUnits(thisBouton).dynamicRange_pooled = range(tuning(:,1)) / options.binSpacing;
    allDotUnits(thisBouton).fanoFactor_pooled   = mean(cellfun(@(x) var(x)/mean(x), alltraces(:,1)), 1);
end

%%  R2 filter + descriptive summaries, POOLED
r2_thresh  = 0.1;
r2p_thresh = 0.05;
validIdx_pooled = find(cat(1,allDotUnits.R2_pooled) > r2_thresh & cat(1,allDotUnits.R2_pooled_pval) < r2p_thresh);
fprintf('\n%d boutons pass the R^2 filter (pooled, R^2>%.2f, p<%.2f).\n', ...
    numel(validIdx_pooled), r2_thresh, r2p_thresh);

%% population R^2 comparison (single CDF, pooled)
figure('Color', 'w', 'Name', 'Population Tuning Quality (pooled, no state split)', 'Position', [200, 200, 500, 400]);
hold on;
allR2Pooled = statR2_pooled(~isnan(statR2_pooled));
hPooled = cdfplot(allR2Pooled);
set(hPooled, 'Color', [0.3 0.3 0.3], 'LineWidth', 2.5, 'DisplayName', 'Pooled (no state split)');
title(sprintf('Cross-validated R^2 distribution (pooled, n=%d boutons)', nBoutonsTotal), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Cross-validated R^2 score', 'FontSize', 11);
ylabel('Proportion of ROIs', 'FontSize', 11);
legend('Location', 'southeast');
hold off;

%% classification distribution (pooled) -- same style as Fig 3.5f
classNames = {'High-pass','Band-pass','Low-pass','Trough','Unclassified'};
idxValid = validIdx_pooled;
classCounts = nan(1,5);
for c = 1:4
    classCounts(c) = sum(gaussChar_pooled_all(idxValid) == c);
end
classCounts(5) = numel(idxValid) - sum(classCounts(1:4));

figure('Color','w','Position',[200 200 450 400]);
bar(classCounts);
set(gca, 'XTickLabel', classNames);
ylabel('# of boutons (pooled, R^2-filtered)');
title(sprintf('Classification distribution, pooled (n=%d)', numel(idxValid)));

%% example single-bouton pooled tuning curves (analogous to Fig 3.5b-e / running examples)
if ~isempty(validIdx_pooled)
    shownIdx = validIdx_pooled(randperm(numel(validIdx_pooled)));
    maxPlots = min(6, numel(shownIdx));
    xDense = linspace(1, 6, 100);
    figure('Color', 'w', 'Name', 'Pooled Gaussian Fit Examples (no state split)', 'Position', [100, 50, 1100, 500]);
    for iPlot = 1:maxPlots
        thisBouton = shownIdx(iPlot);
        subplot(2,3,iPlot); hold on;
        y_data = allDotUnits(thisBouton).tuning(:,1) / options.binSpacing;
        p_pooled = allDotUnits(thisBouton).gaussParams_pooled;
        y_fit = (p_pooled(1) + p_pooled(2)*exp(-((xDense-p_pooled(3)).^2)/(2*p_pooled(4)^2))) / options.binSpacing;

        plot(1:6, y_data, 'ko', 'MarkerFaceColor','k','MarkerSize',5.5);
        plot(xDense, y_fit, 'k-', 'LineWidth', 2);

        allVals = [y_data; y_fit'];
        ylim([min(allVals)-0.03, max(allVals)+0.03]); xlim([0.5,6.5]); xticks(1:6);
        ylabel('\DeltaF/F'); xlabel('Visual Speed Index');
        title(sprintf('Bouton %d (%s)\nPooled Type %d', thisBouton, allDotUnits(thisBouton).sessionLabel, ...
            allDotUnits(thisBouton).gaussChar_pooled), 'FontSize', 9, 'Interpreter', 'none');
    end
    sgtitle('Gaussian Model Fits (pooled across states, no split)', 'FontSize', 13, 'FontWeight', 'bold');
end

%% Summary printout for figure-legend / text use
fprintf('\n=== SUMMARY (paste into text) ===\n');
fprintf('Pooled: %d / %d boutons (%.1f%%) responsive (ANOVA[speeds+blank]+threshold, p<%.2f, >%dSD above blank).\n', ...
    sum(isResponsive), nBoutonsTotal, 100*sum(isResponsive)/nBoutonsTotal, ALPHA, NSD);