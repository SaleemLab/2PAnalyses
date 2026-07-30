% DotFields_TuningCurveAnalysis_compareStatesV2_2PData.m
%
%   - Edd's wheel-based running/stationary trial classification
%     (DotFields_TuningCurveAnalysis_compareStates_2PData.m)
%   - Multi-session/mouse pooling (same pattern as DirTuning scripts)
%   -  blank-trial retention (not discarded) as the literal blank
%     reference, same as DotFields_Pooled_AllTrials.m 
%   - ANOVA(speeds+blank)+2SD-threshold, ANOVA-protected t-test/Wilcoxon,
%     cross-validated R^2, Gaussian fitting, dynamic range/Fano factor --
%     ALL computed SEPARATELY per behavioral state (stationary/running),
%     with blank trials ALSO split by state so each state is compared
%     against its own matched-state blank, not a pooled-across-states one.

%   - Things to remember:
%      - black trials have been assigned stimValue = 1 
%      - as an additional check the number of dots is used (should be 0
%        relative to 576 dots for all other speed catergories) 
%
% Requires calc_kfold_R2.m and fitGaussianTemplates_tuning.m from
% NeuralDataAnalyses (github repo)
% 

%% to do
% check the stimulus durtation across these 6 sessions... were they always
% 2s on 2s off?  yrs

%%
mouseList = {'M26003'};
% mouseList = {'M25126', 'M26004', 'M26005', 'M25131'};

stimWindowMask_range = [0.1 3];  % window for extracting per-trial calcium response
stimFramesMask_range = [-0.2 2.8];  % window for wheel-speed behavior classification

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

options.binSpacing = 2.0; % matches 2s-on/2s-off stimulus timing (used for fano-factor)

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
                trialsSpeed2D(trialCounter).runFlag          = runFlag;
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
        blank_tsd = tsd([tsd.numDots1] == 0);   % genuine BLANK trials, ALSO state-classified

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

        %%  per-bouton, split responses by state (speeds AND blank) 
        for thisROI = 1:nBoutons
            alltraces   = cell(nSpeeds, 2); % column 1 = stationary, column 2 = running
            blankTrials = cell(1, 2); % same coloumn structure as above 

            for thisSpeed = 1:nSpeeds
                for istate = 1:2
                    matchingTrials = find(abs([temp_tsd.VelX1]) == uniqueVelocities(thisSpeed) & ...
                                          [temp_tsd.runFlag] == (istate - 1));
                    traceAccumulator = nan(1, numel(matchingTrials));
                    for mt = 1:numel(matchingTrials)
                        origGroup = temp_tsd(matchingTrials(mt)).origGroup;
                        origTi    = temp_tsd(matchingTrials(mt)).origTrialInGroup;
                        fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
                        traceAccumulator(mt) = nanmean(fullTrace(stimWindowMask));
                    end
                    alltraces{thisSpeed, istate} = traceAccumulator;
                end
            end

            for istate = 1:2
                matchingBlank = find([blank_tsd.runFlag] == (istate - 1));
                blankAccum = nan(1, numel(matchingBlank));
                for mt = 1:numel(matchingBlank)
                    origGroup = blank_tsd(matchingBlank(mt)).origGroup;
                    origTi    = blank_tsd(matchingBlank(mt)).origTrialInGroup;
                    fullTrace = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
                    blankAccum(mt) = nanmean(fullTrace(stimWindowMask));
                end
                blankTrials{istate} = blankAccum;
            end

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

stateNames  = {'stat', 'run'};
stateLabels = {'Stationary', 'Locomotion'};

%%  ANOVA(speeds+blank) + 2SD-above-blank threshold, PER STATE 
for si = 1:2
    thisState = stateNames{si};
    anovaP_thresh    = nan(nBoutonsTotal, 1);
    blankMean        = nan(nBoutonsTotal, 1);
    blankSD          = nan(nBoutonsTotal, 1);
    meanPreferredRaw = nan(nBoutonsTotal, 1);
    prefSpeedIdx     = nan(nBoutonsTotal, 1);
    isResponsive     = false(nBoutonsTotal, 1);

    for thisBouton = 1:nBoutonsTotal
        alltraces   = allDotUnits(thisBouton).alltraces(:, si); % speeds 
        blankTrials = allDotUnits(thisBouton).blankTrials{si}(:); % blank only 
        blankTrials = blankTrials(~isnan(blankTrials));
        nSpeeds = numel(alltraces);
        % every single trial's response value, from every speed condition and the blank condition
        y = [];  
        % same length as y telling which condition each entry in y came from
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

    fprintf('\n[%s] %d / %d boutons responsive (ANOVA[speeds+blank]+threshold, p<%.2f, >%dSD above blank).\n', ...
        stateLabels{si}, sum(isResponsive), nBoutonsTotal, ALPHA, NSD);

    for thisBouton = 1:nBoutonsTotal
        allDotUnits(thisBouton).(['anovaP_thresh_' thisState])    = anovaP_thresh(thisBouton);
        allDotUnits(thisBouton).(['blankMean_' thisState])        = blankMean(thisBouton);
        allDotUnits(thisBouton).(['blankSD_' thisState])          = blankSD(thisBouton);
        allDotUnits(thisBouton).(['meanPreferredRaw_' thisState]) = meanPreferredRaw(thisBouton);
        allDotUnits(thisBouton).(['prefSpeedIdx_' thisState])     = prefSpeedIdx(thisBouton);
        allDotUnits(thisBouton).(['isResponsive_' thisState])     = isResponsive(thisBouton); 
    end
end


%% check trial counts by states
%DotFields_CheckTrialCountsByState
%% are boutons responsive during locmotion also responsive during stationary
DotFields_ResponsiveOverlapCheck
%%  ANOVA-protected t-test / Wilcoxon vs blank, PER STATE 
for si = 1:2
    thisState = stateNames{si};
    ttestP    = nan(nBoutonsTotal, 1);
    ranksumP  = nan(nBoutonsTotal, 1);
    isResp_ttest   = false(nBoutonsTotal, 1);
    isResp_ranksum = false(nBoutonsTotal, 1);

    for thisBouton = 1:nBoutonsTotal
        anovaP_thresh_b = allDotUnits(thisBouton).(['anovaP_thresh_' thisState]);
        prefSpeedIdx_b  = allDotUnits(thisBouton).(['prefSpeedIdx_' thisState]);
        if isnan(anovaP_thresh_b) || isnan(prefSpeedIdx_b)
            continue;
        end

        prefTrials  = allDotUnits(thisBouton).alltraces{prefSpeedIdx_b, si}(:);
        prefTrials  = prefTrials(~isnan(prefTrials));
        blankTrials = allDotUnits(thisBouton).blankTrials{si}(:);
        blankTrials = blankTrials(~isnan(blankTrials));

        if numel(prefTrials) < 2 || numel(blankTrials) < 2
            continue;
        end

        [~, ttestP(thisBouton)] = ttest2(prefTrials, blankTrials, 'Vartype', 'unequal');
        ranksumP(thisBouton) = ranksum(prefTrials, blankTrials);

        blankMean_b = allDotUnits(thisBouton).(['blankMean_' thisState]);
        meanPreferredRaw_b = allDotUnits(thisBouton).(['meanPreferredRaw_' thisState]);

        isResp_ttest(thisBouton)   = (anovaP_thresh_b < ALPHA) && (ttestP(thisBouton) < ALPHA) && (meanPreferredRaw_b > blankMean_b);
        isResp_ranksum(thisBouton) = (anovaP_thresh_b < ALPHA) && (ranksumP(thisBouton) < ALPHA) && (meanPreferredRaw_b > blankMean_b);
    end

    fprintf('[%s] %d / %d boutons responsive (ANOVA-protected Welch t-test, p<%.2f).\n', ...
        stateLabels{si}, sum(isResp_ttest), nBoutonsTotal, ALPHA);
    fprintf('[%s] %d / %d boutons responsive (ANOVA-protected Wilcoxon, p<%.2f).\n', ...
        stateLabels{si}, sum(isResp_ranksum), nBoutonsTotal, ALPHA);

    for thisBouton = 1:nBoutonsTotal
        allDotUnits(thisBouton).(['ttestP_' thisState])             = ttestP(thisBouton);
        allDotUnits(thisBouton).(['ranksumP_' thisState])           = ranksumP(thisBouton);
        allDotUnits(thisBouton).(['isResponsive_ttest_' thisState])   = isResp_ttest(thisBouton);
        allDotUnits(thisBouton).(['isResponsive_ranksum_' thisState]) = isResp_ranksum(thisBouton);
    end
end

%% method comparison summary, PER STATE
for si = 1:2
    thisState = stateNames{si};
    fprintf('\n=== Method comparison [%s] (n=%d boutons) ===\n', stateLabels{si}, nBoutonsTotal);
    fprintf('%-35s %6d\n', 'SD-heuristic',                 sum([allDotUnits.(['isResponsive_' thisState])]));
    fprintf('%-35s %6d\n', 'ANOVA-protected Welch t-test',  sum([allDotUnits.(['isResponsive_ttest_' thisState])]));
    fprintf('%-35s %6d\n', 'ANOVA-protected Wilcoxon',      sum([allDotUnits.(['isResponsive_ranksum_' thisState])]));
end


%%
reliability_stat = nan(nBoutonsTotal, 1);
reliability_run  = nan(nBoutonsTotal, 1);

for thisBouton = 1:nBoutonsTotal
    reliability_stat(thisBouton) = computeVisualReliabilityIndex_DotFields_fromScalars(thisBouton, allDotUnits, 1);
    reliability_run(thisBouton)  = computeVisualReliabilityIndex_DotFields_fromScalars(thisBouton, allDotUnits, 2);
end

figure('Color','w');
subplot(1,2,1); histogram(reliability_stat, 30); xlabel('Reliability'); title('Dot Fields - Stationary');
subplot(1,2,2); histogram(reliability_run, 30); xlabel('Reliability'); title('Dot Fields - Locomotion');

%% cross-validated R^2, PER STATE -- Edd's joint per-session downsample + a fixed global floor
% Edd's original (DotFields_TuningCurveAnalysis_compareStates_2PData.m,
% run once per session) computes ONE minTrial per session:
%   minTrial = min(min(cellfun(@(x) size(x,2), units(1).allSpikes)));
% This is a single value taken jointly across ALL speeds AND BOTH states
% (not independent per state), using ROI 1 as a stand-in for the whole
% session -- valid because matchingTrials is built purely from session
% trial metadata (temp_tsd), identically for every ROI, so every ROI in a
% session has exactly the same trial COUNT per condition (individual
% values can still be NaN, but the array length can't differ). He then
% truncates every ROI's raw trial columns with x(:,1:minTrial) -- no NaN
% filtering, just the first minTrial trials as recorded.
%
% Our per-session minTrial ranged from 2 to 12 trials across  6
% sessions (checked via the per-session/per-condition diagnostic), which
% is too variable to trust R^2 magnitudes as comparable across sessions.
% We therefore ALSO apply a fixed global floor (nTrialsTarget = 9) on top
% of Edd's per-session logic:
%   - sessions whose own minTrial < nTrialsTarget are excluded entirely
%     (every bouton in that session gets NaN for both statR2 and runR2)
%   - sessions whose minTrial >= nTrialsTarget are truncated down to
%     EXACTLY nTrialsTarget (not their own higher available count)
% This keeps the stat-vs-run comparison fair (matched trial count per
% bouton, as before) AND keeps R^2 comparable across the sessions that
% remain (matched trial count across sessions too). With nTrialsTarget=9
% this retains 3 of 6 sessions (M25132_20260306, M26003_20260324,
% M26003_20260326); nTrialsTarget=10 would only retain 2.

nTrialsTarget = 4; % at least n trials included in each session 

sessionLabels_all = {allDotUnits.sessionLabel}';
uniqueSessionsAll  = unique(sessionLabels_all, 'stable');
sessionMinTrial    = containers.Map('KeyType', 'char', 'ValueType', 'double');

for s = 1:numel(uniqueSessionsAll)
    thisSession     = uniqueSessionsAll{s};
    firstBoutonIdx  = find(strcmp(sessionLabels_all, thisSession), 1);
    alltraces_ref   = allDotUnits(firstBoutonIdx).alltraces; % nSpeeds x 2, raw (unfiltered) trial counts
    sessionMinTrial(thisSession) = min(min(cellfun(@numel, alltraces_ref)));
end

fprintf('\nSessions retained at nTrialsTarget=%d:\n', nTrialsTarget);
for s = 1:numel(uniqueSessionsAll)
    thisSession = uniqueSessionsAll{s};
    thisMin = sessionMinTrial(thisSession);
    fprintf('  %s: minTrial=%d -> %s\n', thisSession, thisMin, ...
        string(thisMin >= nTrialsTarget));
end

statR2 = nan(nBoutonsTotal, 1); statR2_pval = nan(nBoutonsTotal, 1);
runR2  = nan(nBoutonsTotal, 1); runR2_pval  = nan(nBoutonsTotal, 1);
seed = 1;
for thisBouton = 1:nBoutonsTotal
    alltraces = allDotUnits(thisBouton).alltraces; % nSpeeds x 2
    minTrial  = sessionMinTrial(allDotUnits(thisBouton).sessionLabel);

    if minTrial < nTrialsTarget
        continue; % session excluded entirely -- doesn't meet the fixed global floor
    end

    alltracesDS = cellfun(@(x) x(:, 1:nTrialsTarget), alltraces, 'UniformOutput', false);

    gca_stat = alltracesDS(:, 1)';
    [statR2(thisBouton), statR2_pval(thisBouton)] = calc_kfold_R2(gca_stat, r2opts.kval, r2opts.nPerms, ...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle, seed);

    gca_run = alltracesDS(:, 2)';
    [runR2(thisBouton), runR2_pval(thisBouton)] = calc_kfold_R2(gca_run, r2opts.kval, r2opts.nPerms, ...
        r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle, seed);

    if mod(thisBouton, 100) == 0
        fprintf('Cross-val R^2: %d / %d boutons...\n', thisBouton, nBoutonsTotal);
    end
end

isTunedStat = statR2_pval < ALPHA;
isTunedRun  = runR2_pval  < ALPHA;

nStatValid = sum(~isnan(statR2_pval));
nRunValid  = sum(~isnan(runR2_pval));

fprintf('\n[Stationary] %d / %d boutons show significant cross-validated tuning (p<%.2f).\n', ...
    sum(isTunedStat, 'omitnan'), nBoutonsTotal, ALPHA);
fprintf('[Locomotion] %d / %d boutons show significant cross-validated tuning (p<%.2f).\n', ...
    sum(isTunedRun, 'omitnan'), nBoutonsTotal, ALPHA);

% Real proportions: out of boutons that had enough trials to get an R^2
% in that state at all, what fraction were tuned. With the joint
% per-session downsample (matching Edd's method), nStatValid and
% nRunValid will be equal -- this is still worth checking so you report
% percentages relative to the actual included pool, not nBoutonsTotal.
fprintf('[Stationary] %d / %d boutons WITH VALID R^2 were tuned (%.1f%%).\n', ...
    sum(isTunedStat, 'omitnan'), nStatValid, 100*sum(isTunedStat,'omitnan')/nStatValid);
fprintf('[Locomotion] %d / %d boutons WITH VALID R^2 were tuned (%.1f%%).\n', ...
    sum(isTunedRun, 'omitnan'), nRunValid, 100*sum(isTunedRun,'omitnan')/nRunValid);



for thisBouton = 1:nBoutonsTotal
    allDotUnits(thisBouton).statR2      = statR2(thisBouton);
    allDotUnits(thisBouton).statR2_pval = statR2_pval(thisBouton);
    allDotUnits(thisBouton).runR2       = runR2(thisBouton);
    allDotUnits(thisBouton).runR2_pval  = runR2_pval(thisBouton);
    allDotUnits(thisBouton).isTunedStat = isTunedStat(thisBouton);
    allDotUnits(thisBouton).isTunedRun  = isTunedRun(thisBouton);
end



%% cross-validated R^2, PER STATE (downsample WITHIN each state only)
% statR2 = nan(nBoutonsTotal, 1); statR2_pval = nan(nBoutonsTotal, 1);
% runR2  = nan(nBoutonsTotal, 1); runR2_pval  = nan(nBoutonsTotal, 1);
% seed = 1;
% for thisBouton = 1:nBoutonsTotal
%     alltraces = allDotUnits(thisBouton).alltraces; % nSpeeds x 2
% 
%     %  stationary: downsample using ONLY stationary's own 6 speed cells 
%     minTrial_stat = min(cellfun(@(x) sum(~isnan(x)), alltraces(:,1)));
%     if minTrial_stat >= r2opts.kval
%         gca_stat = cellfun(@(x) x(find(~isnan(x), minTrial_stat)), alltraces(:,1), 'UniformOutput', false)';
%         [statR2(thisBouton), statR2_pval(thisBouton)] = calc_kfold_R2(gca_stat, r2opts.kval, r2opts.nPerms, ...
%             r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle, seed);
%     end
% 
%     % locomotion: downsample using ONLY locomotion's own 6 speed cells
%     minTrial_run = min(cellfun(@(x) sum(~isnan(x)), alltraces(:,2)));
%     if minTrial_run >= r2opts.kval
%         gca_run = cellfun(@(x) x(find(~isnan(x), minTrial_run)), alltraces(:,2), 'UniformOutput', false)';
%         [runR2(thisBouton), runR2_pval(thisBouton)] = calc_kfold_R2(gca_run, r2opts.kval, r2opts.nPerms, ...
%             r2opts.randFlag, r2opts.validMeans, r2opts.nShuffle, seed);
%     end
% 
%     if mod(thisBouton, 100) == 0
%         fprintf('Cross-val R^2: %d / %d boutons...\n', thisBouton, nBoutonsTotal);
%     end
% end
% 
% 
% 
% isTunedStat = statR2_pval < ALPHA;
% isTunedRun  = runR2_pval < ALPHA;
% fprintf('\n[Stationary] %d / %d boutons show significant cross-validated tuning (p<%.2f).\n', ...
%     sum(isTunedStat, 'omitnan'), nBoutonsTotal, ALPHA);
% fprintf('[Locomotion] %d / %d boutons show significant cross-validated tuning (p<%.2f).\n', ...
%     sum(isTunedRun, 'omitnan'), nBoutonsTotal, ALPHA);
% 
% for thisBouton = 1:nBoutonsTotal
%     allDotUnits(thisBouton).statR2      = statR2(thisBouton);
%     allDotUnits(thisBouton).statR2_pval = statR2_pval(thisBouton);
%     allDotUnits(thisBouton).runR2       = runR2(thisBouton);
%     allDotUnits(thisBouton).runR2_pval  = runR2_pval(thisBouton);
%     allDotUnits(thisBouton).isTunedStat = isTunedStat(thisBouton);
%     allDotUnits(thisBouton).isTunedRun  = isTunedRun(thisBouton);
% end
% 
% r2_thresh  = 0.0;
% r2p_thresh = 0.05;
% 
% validIdx_stat = find(cat(1,allDotUnits.statR2) > r2_thresh & cat(1,allDotUnits.statR2_pval) < r2p_thresh);
% validIdx_run  = find(cat(1,allDotUnits.runR2)  > r2_thresh & cat(1,allDotUnits.runR2_pval)  < r2p_thresh);
% validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));
% 
% nStatOnly = numel(setdiff(validIdx_stat, validIdx_run));
% nRunOnly  = numel(setdiff(validIdx_run, validIdx_stat));
% nBoth     = numel(validIdx_both);
% nStatTotal = numel(validIdx_stat);
% nRunTotal  = numel(validIdx_run);
% nBoutonsTotal = numel(allDotUnits);
% nNeither  = nBoutonsTotal - numel(union(validIdx_stat, validIdx_run));
% 
% fprintf('\n--- R^2 filter breakdown (R^2>%.2f, p<%.2f) ---\n', r2_thresh, r2p_thresh);
% fprintf('%-30s %6d\n', 'Total boutons',            nBoutonsTotal);
% fprintf('%-30s %6d\n', 'Pass stationary only',      nStatOnly);
% fprintf('%-30s %6d\n', 'Pass locomotion only',      nRunOnly);
% fprintf('%-30s %6d\n', 'Pass BOTH (paired cohort)', nBoth);
% fprintf('%-30s %6d\n', 'Pass neither',              nNeither);
% fprintf('%-30s %6d\n', 'Total passing stationary',  nStatTotal);
% fprintf('%-30s %6d\n', 'Total passing locomotion',  nRunTotal);

%% check how many stationry vs running trials there rae
statTrials = arrayfun(@(x) min(cellfun(@(y) sum(~isnan(y)), x.alltraces(:,1))), allDotUnits);
runTrials  = arrayfun(@(x) min(cellfun(@(y) sum(~isnan(y)), x.alltraces(:,2))), allDotUnits);

figure('Color','w');
subplot(1,2,1); histogram(statTrials); title('Min trials/speed - Stationary'); xlabel('N trials');
subplot(1,2,2); histogram(runTrials);  title('Min trials/speed - Locomotion'); xlabel('N trials');
fprintf('Median min-trials: stationary=%.1f, locomotion=%.1f\n', median(statTrials), median(runTrials));

%%  population R^2 comparison (CDF, stat vs run) 
figure('Color', 'w', 'Name', 'Population Tuning Quality (pooled)', 'Position', [200, 200, 500, 400]);
hold on;
allStatR2 = statR2(~isnan(statR2));
allRunR2  = runR2(~isnan(runR2));

nStat = numel(allStatR2);
nRun  = numel(allRunR2);

hStat = cdfplot(allStatR2);
set(hStat, 'Color', 'k', 'LineWidth', 2.5, 'DisplayName', sprintf('Stationary (n=%d)', nStat));
hRun = cdfplot(allRunR2);
set(hRun, 'Color', 'r', 'LineWidth', 2.5, 'DisplayName', sprintf('Locomotion (n=%d)', nRun));

title(sprintf('Cross-validated R^2 distribution (pooled, n=%d boutons)', nBoutonsTotal), 'FontSize', 12, 'FontWeight', 'bold');
subtitle('Numbers in brackets indicate number of contributing tuning curves', 'FontSize', 9, 'FontAngle', 'italic');
xlabel('Cross-validated R^2 score', 'FontSize', 11);
ylabel('Proportion of ROIs', 'FontSize', 11);
legend('Location', 'southeast');
hold off;

%%  Gaussian tuning-curve fits, PER STATE 
gaussParams_stat_all = repmat({nan(1,4)}, nBoutonsTotal, 1);
gaussChar_stat_all   = nan(nBoutonsTotal, 1);
gaussR2_stat_all     = nan(nBoutonsTotal, 1);
prefSpeed_stat_all   = nan(nBoutonsTotal, 1);

gaussParams_run_all  = repmat({nan(1,4)}, nBoutonsTotal, 1);
gaussChar_run_all    = nan(nBoutonsTotal, 1);
gaussR2_run_all      = nan(nBoutonsTotal, 1);
prefSpeed_run_all    = nan(nBoutonsTotal, 1);

for thisBouton = 1:nBoutonsTotal
    alltraces = allDotUnits(thisBouton).alltraces;
    tuning    = allDotUnits(thisBouton).tuning;

    alltracesClean = cellfun(@(x) x(~isnan(x)), alltraces, 'UniformOutput', false);
    if any(cellfun(@isempty, alltracesClean(:)))
        continue;
    end

    istate = 1;
    [gp_stat, gc_stat, gr2_stat] = fitGaussianTemplates_tuning(alltracesClean(:,istate), 0.5, false);
    if gc_stat == 4, [~, ps_stat] = min(tuning(:,istate)); else, [~, ps_stat] = max(tuning(:,istate)); end

    istate = 2;
    [gp_run, gc_run, gr2_run] = fitGaussianTemplates_tuning(alltracesClean(:,istate), 0.5, false);
    if gc_run == 4, [~, ps_run] = min(tuning(:,istate)); else, [~, ps_run] = max(tuning(:,istate)); end

    gaussParams_stat_all{thisBouton} = gp_stat; gaussChar_stat_all(thisBouton) = gc_stat; gaussR2_stat_all(thisBouton) = gr2_stat; prefSpeed_stat_all(thisBouton) = ps_stat;
    gaussParams_run_all{thisBouton}  = gp_run;  gaussChar_run_all(thisBouton)  = gc_run;  gaussR2_run_all(thisBouton)  = gr2_run;  prefSpeed_run_all(thisBouton)  = ps_run;

    if mod(thisBouton, 100) == 0
        fprintf('Gaussian fits: %d / %d boutons...\n', thisBouton, nBoutonsTotal);
    end
end

for thisBouton = 1:nBoutonsTotal
    allDotUnits(thisBouton).gaussParams_stat = gaussParams_stat_all{thisBouton};
    allDotUnits(thisBouton).gaussChar_stat   = gaussChar_stat_all(thisBouton);
    allDotUnits(thisBouton).gaussR2_stat     = gaussR2_stat_all(thisBouton);
    allDotUnits(thisBouton).prefSpeed_stat   = prefSpeed_stat_all(thisBouton);

    allDotUnits(thisBouton).gaussParams_run  = gaussParams_run_all{thisBouton};
    allDotUnits(thisBouton).gaussChar_run    = gaussChar_run_all(thisBouton);
    allDotUnits(thisBouton).gaussR2_run      = gaussR2_run_all(thisBouton);
    allDotUnits(thisBouton).prefSpeed_run    = prefSpeed_run_all(thisBouton);
end

%% dynamic range and Fano factor, PER STATE
for thisBouton = 1:nBoutonsTotal
    alltraces = allDotUnits(thisBouton).alltraces;
    tuning    = allDotUnits(thisBouton).tuning;

    istate = 1;
    allDotUnits(thisBouton).dynamicRange_stat = range(tuning(:,istate)) / options.binSpacing;
    allDotUnits(thisBouton).fanoFactor_stat   = mean(cellfun(@(x) var(x)/mean(x), alltraces(:,istate)), 1);

    istate = 2;
    allDotUnits(thisBouton).dynamicRange_run = range(tuning(:,istate)) / options.binSpacing;
    allDotUnits(thisBouton).fanoFactor_run   = mean(cellfun(@(x) var(x)/mean(x), alltraces(:,istate)), 1);
end

%%  R2 filter for both stat and loco + descriptive summaries 

r2_thresh  = 0.1;
r2p_thresh = 0.05;
validIdx_stat = find(cat(1,allDotUnits.statR2) > r2_thresh & cat(1,allDotUnits.statR2_pval) < r2p_thresh);
validIdx_run  = find(cat(1,allDotUnits.runR2)  > r2_thresh & cat(1,allDotUnits.runR2_pval)  < r2p_thresh);
validIdx_both = validIdx_stat(ismember(validIdx_stat, validIdx_run));
fprintf('\n%d boutons pass the dual R^2 filter (both stationary AND running, R^2>%.2f, p<%.2f).\n', ...
    numel(validIdx_both), r2_thresh, r2p_thresh);

if ~isempty(validIdx_both)
    allPrefs = [cat(1,allDotUnits(validIdx_both).prefSpeed_stat), cat(1,allDotUnits(validIdx_both).prefSpeed_run)];
    vals = histcounts2(allPrefs(:,1), allPrefs(:,2), 1:7, 1:7);
    vals = vals./sum(vals,2);
    figure
    imagesc(vals')
    hold on
    plot([0.5 6.5],[0.5 6.5],'r')
    axis xy
    colorbar
    xlabel('Preferred stimulus (stationary)')
    ylabel('Preferred stimulus (locomotion)')
    title(sprintf('Preferred speed agreement (n=%d cross-validated units)', numel(validIdx_both)));

    figure
    plot([allDotUnits(validIdx_both).dynamicRange_stat], [allDotUnits(validIdx_both).dynamicRange_run], 'k.')
    xlabel('Dynamic range (stationary)'); ylabel('Dynamic range (locomotion)');
    title(sprintf('Dynamic range (n=%d)', numel(validIdx_both)));

    figure
    plot([allDotUnits(validIdx_both).fanoFactor_stat], [allDotUnits(validIdx_both).fanoFactor_run], 'k.')
    xlabel('Fano Factor (stationary)'); ylabel('Fano Factor (locomotion)');
    title(sprintf('Fano factor (n=%d)', numel(validIdx_both)));
end

idx = find(cat(1,allDotUnits.runR2)>r2_thresh & cat(1,allDotUnits.runR2_pval)<r2p_thresh...
    & cat(1,allDotUnits.statR2)>r2_thresh & cat(1,allDotUnits.statR2_pval)<r2p_thresh...
    & cat(1,allDotUnits.gaussChar_stat)==3 & cat(1,allDotUnits.gaussChar_run)==3);
if ~isempty(idx)
    allStatParams = cat(1,allDotUnits(idx).gaussParams_stat);
    allRunParams  = cat(1,allDotUnits(idx).gaussParams_run);
    paramNames = {'Baseline', 'Response amplitude', '\mu (preferred speed location)', '\sigma (tuning width)'};
    paramColors = {[0,0.4470,0.7410],[0.85,0.325,0.098],[0.929,0.694,0.125],[0.494,0.184,0.556]};

    figure('Color', 'w', 'Position', [100, 100, 900, 750]);
    for iparam = 1:4
        subplot(2,2,iparam); hold on;
        x_stat = allStatParams(:,iparam); y_run = allRunParams(:,iparam);
        if iparam == 1 || iparam == 2
            x_stat = x_stat / options.binSpacing; y_run = y_run / options.binSpacing;
        end
        scatter(x_stat, y_run, 45, 'MarkerFaceColor', paramColors{iparam}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
        minVal = min([x_stat; y_run]); maxVal = max([x_stat; y_run]);
        if minVal == maxVal, minVal = minVal-0.1; maxVal = maxVal+0.1; end
        plot([minVal,maxVal],[minVal,maxVal],'k--','LineWidth',1.2);
        p_val = signrank(x_stat, y_run);
        xlim([minVal,maxVal]); ylim([minVal,maxVal]);
        xlabel('Stationary Value'); ylabel('Locomotion Value');
        title(sprintf('%s (signrank p=%.3g)', paramNames{iparam}, p_val));
    end
    sgtitle(sprintf('Gaussian params: bandpass cells only (n=%d)', numel(idx)));
end
validIdx_both = validIdx_both(randperm(length(validIdx_both)));
if ~isempty(validIdx_both)
    maxPlots = min(6, numel(validIdx_both));
    xDense = linspace(1, 6, 100);
    figure('Color', 'w', 'Name', 'Cross-Validated Gaussian Fit Overlays', 'Position', [100, 50, 1100, 800]);
    for iPlot = 1:maxPlots
        
        thisBouton = validIdx_both(iPlot);
        subplot(2,3,iPlot); hold on;
        y_data_stat = allDotUnits(thisBouton).tuning(:,1) / options.binSpacing;
        y_data_run  = allDotUnits(thisBouton).tuning(:,2) / options.binSpacing;
        p_stat = allDotUnits(thisBouton).gaussParams_stat;
        y_fit_stat = (p_stat(1) + p_stat(2)*exp(-((xDense-p_stat(3)).^2)/(2*p_stat(4)^2))) / options.binSpacing;
        p_run = allDotUnits(thisBouton).gaussParams_run;
        y_fit_run = (p_run(1) + p_run(2)*exp(-((xDense-p_run(3)).^2)/(2*p_run(4)^2))) / options.binSpacing;

        plot(1:6, y_data_stat, 'ko', 'MarkerFaceColor','k','MarkerSize',5.5);
        plot(xDense, y_fit_stat, 'k-', 'LineWidth', 2, 'DisplayName', 'Stationary');
        plot(1:6, y_data_run, 'ro', 'MarkerFaceColor','r','MarkerSize',5.5);
        plot(xDense, y_fit_run, 'r-', 'LineWidth', 2, 'DisplayName', 'Locomotion');

        allVals = [y_data_stat; y_data_run; y_fit_stat'; y_fit_run'];
        ylim([min(allVals)-0.03, max(allVals)+0.03]); xlim([0.5,6.5]); xticks(1:6);
        ylabel('\DeltaF/F'); xlabel('Visual Speed Index');
        title(sprintf('Bouton %d (%s)\nStat Type %d | Run Type %d', thisBouton, allDotUnits(thisBouton).sessionLabel, ...
            allDotUnits(thisBouton).gaussChar_stat, allDotUnits(thisBouton).gaussChar_run), 'FontSize', 9, 'Interpreter', 'none');
        if iPlot == 1, legend('Location','best'); end
    end
    sgtitle('Gaussian Model Fits Across States (cross-validated cohort)', 'FontSize', 13, 'FontWeight', 'bold');
end

%%
r2_thresh  = 0.1;
r2p_thresh = 0.05;
 
nBoutonsTotal = numel(allDotUnits);
 
stateFieldSuffix = {'stat', 'run'};
stateTitles      = {'Stationary', 'Running'};
 
for si = 1:2
    suffix = stateFieldSuffix{si};
 
    cvR2   = cat(1, allDotUnits.([suffix 'R2']));
    cvPval = cat(1, allDotUnits.([suffix 'R2_pval']));
 
    validIdx = find(cvR2 > r2_thresh & cvPval < r2p_thresh);
    fprintf('\n[%s] %d / %d boutons pass R^2 filter (R^2>%.2f, p<%.2f) for descriptive summaries.\n', ...
        stateTitles{si}, numel(validIdx), nBoutonsTotal, r2_thresh, r2p_thresh);
 
    gaussChar    = cat(1, allDotUnits.(['gaussChar_' suffix]));
    prefSpeed    = cat(1, allDotUnits.(['prefSpeed_' suffix]));
    dynamicRange = cat(1, allDotUnits.(['dynamicRange_' suffix]));
    fanoFactor   = cat(1, allDotUnits.(['fanoFactor_' suffix]));
 
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
    histogram(dynamicRange(validIdx));
    xlabel('Dynamic range'); ylabel('Count');
    title('Dynamic range distribution');
 
    subplot(2,2,4);
    histogram(fanoFactor(validIdx));
    xlabel('Fano factor'); ylabel('Count');
    title('Fano factor distribution');
 
    sgtitle(sprintf('Descriptive tuning summaries [%s] (R^2-filtered, n=%d)', ...
        stateTitles{si}, numel(validIdx)));
end

% tests
% DotFields_SelectRepresentativeBouton



% plotDotFieldsExampleBoutonStateSplit('M25132', '20260226',12)
% plotDotFieldsExampleBoutonStateSplit('M25132', '20260226',9)
% plotDotFieldsExampleBoutonStateSplit('M25133', '20260224',38)
% 
% plotDotFieldsExampleBoutonStateSplit('M26003', '20260326',20)


%DotFields_TrialCountMatchedComparison

