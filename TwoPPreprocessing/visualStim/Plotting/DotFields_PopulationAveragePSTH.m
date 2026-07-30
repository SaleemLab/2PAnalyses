% DotFields_PopulationAveragePSTH.m
%
% Builds a population-averaged PSTH (full time-course, not just a
% windowed mean) for SD-heuristic responsive boutons, separately for
% stationary and running trials, using each bouton's OWN preferred speed
% (so a bouton's contribution to the average reflects its best response,
% not diluted by non-preferred speeds).
%
% Purpose: eyeball rise/peak/decay timing across the population (not just
% one example bouton) to choose a defensible response window for the
% ANOVA/threshold responsiveness test and tuning-curve fits.
%
% Requires the SAME mouseList / classification thresholds as
% DotFields_TuningCurveAnalysis_compareStatesV2_2PData.m -- kept
% identical here so the population this figure describes matches the
% population your responsiveness numbers come from.

%%
mouseList = {'M25132', 'M25133', 'M26003'};

stimFramesMask_range = [-0.2 3];   % window for wheel-speed behavior classification
runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;

ALPHA = 0.05;
NSD   = 1;

% NOTE: this is the window used ONLY to pick each bouton's preferred
% speed / determine responsiveness for inclusion here -- deliberately
% kept the same as your main script's stimWindowMask_range so the
% "responsive" bouton pool matches. The whole point of this script is to
% then look at the FULL PSTH shape, unconstrained by that window, to
% sanity-check whether the window itself is well chosen.
prefSpeedWindow_range = [0.5 2];

%%
filteredTable = filterMasterTable('MouseID', mouseList, 'HasStimulus', 'DotMotion_SpeedTuning', ...
    'Suite2PPreprocessing', 1, 'Exclude', 0);
allMice    = filteredTable.MouseID;
uniqueMice = unique(allMice, 'stable');

% accumulate: for each responsive bouton, its full PSTH trace AT ITS
% PREFERRED SPEED, per state
psthTraces_stat = {}; % each cell: [1 x nTimepoints], mean across that bouton's preferred-speed trials
psthTraces_run  = {};
timeVecRef = []; % populated once, assumed consistent across sessions

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

        %% classify every trial (incl. blank), same as main script
        nGroups = numel(response.wheelData);
        trialsSpeed2D = struct('VelX1', {}, 'numDots1', {}, 'runFlag', {}, ...
                                'origGroup', {}, 'origTrialInGroup', {});
        trialCounter = 1;
        for g = 1:nGroups
            grpWheel  = response.wheelData(g);
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
                if isnan(runFlag), continue; end
                trialsSpeed2D(trialCounter).VelX1 = grpWheel.stimValue;
                trialsSpeed2D(trialCounter).numDots1 = (grpWheel.stimValue ~= 1) * 573;
                trialsSpeed2D(trialCounter).runFlag          = runFlag;
                trialsSpeed2D(trialCounter).origGroup        = g;
                trialsSpeed2D(trialCounter).origTrialInGroup = ti;
                trialCounter = trialCounter + 1;
            end
        end
        if trialCounter == 1, continue; end

        tsd = trialsSpeed2D;
        temp_tsd  = tsd([tsd.numDots1] == 573);
        blank_tsd = tsd([tsd.numDots1] == 0);
        if isempty(temp_tsd), continue; end

        uniqueVelocities = unique(abs([temp_tsd.VelX1]));
        nSpeeds = numel(uniqueVelocities);
        nBoutons = size(response.psthData(1).alignedResponses, 1);
        timeVec = response.psthData(1).timeVector(:)';
        if isempty(timeVecRef), timeVecRef = timeVec; end
        prefWinMask = timeVec >= prefSpeedWindow_range(1) & timeVec <= prefSpeedWindow_range(2);

        %% per bouton: find preferred speed + responsiveness per state, then grab full PSTH at that speed
        for thisROI = 1:nBoutons
            for istate = 0:1
                % gather windowed means per speed (to find preferred speed + run ANOVA/threshold)
                meansPerSpeed = nan(nSpeeds, 1);
                fullTracesPerSpeed = cell(nSpeeds, 1);
                yAll = []; grpAll = [];
                for s = 1:nSpeeds
                    matchIdx = find(abs([temp_tsd.VelX1]) == uniqueVelocities(s) & [temp_tsd.runFlag] == istate);
                    traceMat = nan(numel(matchIdx), numel(timeVec));
                    for mt = 1:numel(matchIdx)
                        og = temp_tsd(matchIdx(mt)).origGroup;
                        ot = temp_tsd(matchIdx(mt)).origTrialInGroup;
                        traceMat(mt,:) = squeeze(response.psthData(og).alignedResponses(thisROI, :, ot));
                    end
                    fullTracesPerSpeed{s} = traceMat;
                    windowedVals = mean(traceMat(:, prefWinMask), 2, 'omitnan');
                    meansPerSpeed(s) = mean(windowedVals, 'omitnan');
                    yAll = [yAll; windowedVals(~isnan(windowedVals))]; %#ok<AGROW>
                    grpAll = [grpAll; repmat(s, sum(~isnan(windowedVals)), 1)]; %#ok<AGROW>
                end

                blankIdx = find([blank_tsd.runFlag] == istate);
                blankVals = nan(numel(blankIdx), 1);
                for mt = 1:numel(blankIdx)
                    og = blank_tsd(blankIdx(mt)).origGroup;
                    ot = blank_tsd(blankIdx(mt)).origTrialInGroup;
                    fullTrace = squeeze(response.psthData(og).alignedResponses(thisROI, :, ot));
                    blankVals(mt) = mean(fullTrace(prefWinMask), 'omitnan');
                end
                blankVals = blankVals(~isnan(blankVals));
                if isempty(blankVals) || numel(yAll) < (nSpeeds+1)*2, continue; end

                yAll   = [yAll; blankVals]; %#ok<AGROW>
                grpAll = [grpAll; repmat(nSpeeds+1, numel(blankVals), 1)]; %#ok<AGROW>
                pAnova = anova1(yAll, grpAll, 'off');

                [maxVal, prefIdx] = max(meansPerSpeed);
                isResp = (pAnova < ALPHA) && (maxVal > mean(blankVals) + NSD*std(blankVals));
                if ~isResp, continue; end

                prefTraceMat = fullTracesPerSpeed{prefIdx};
                if isempty(prefTraceMat), continue; end
                meanTrace = mean(prefTraceMat, 1, 'omitnan');

                if istate == 0
                    psthTraces_stat{end+1} = meanTrace; %#ok<SAGROW>
                else
                    psthTraces_run{end+1} = meanTrace; %#ok<SAGROW>
                end
            end
        end
        fprintf('    Running totals -- stat: %d | run: %d responsive-bouton PSTHs\n', ...
            numel(psthTraces_stat), numel(psthTraces_run));
    end
end

%% assemble into matrices + plot population average +/- SEM
matStat = cat(1, psthTraces_stat{:});
matRun  = cat(1, psthTraces_run{:});

meanStat = mean(matStat, 1, 'omitnan'); semStat = std(matStat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(matStat),1));
meanRun  = mean(matRun, 1, 'omitnan');  semRun  = std(matRun, 0, 1, 'omitnan')  ./ sqrt(sum(~isnan(matRun),1));

figure('Color','w','Position',[200 200 800 450]);
popPSTHFig = gcf;
hold on;
fill([timeVecRef, fliplr(timeVecRef)], [meanStat+semStat, fliplr(meanStat-semStat)], 'k', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
plot(timeVecRef, meanStat, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Stationary (n=%d)', size(matStat,1)));
fill([timeVecRef, fliplr(timeVecRef)], [meanRun+semRun, fliplr(meanRun-semRun)], 'r', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');
plot(timeVecRef, meanRun, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Running (n=%d)', size(matRun,1)));
xline(0, 'k--', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('\DeltaF/F (mean across responsive boutons, at each bouton''s preferred speed)');
title('Population-averaged PSTH (SD-heuristic responsive boutons, preferred speed)');
legend('Location','best');
box off;

%% quantify rise/peak/decay to help pick a window
[peakValStat, peakIdxStat] = max(meanStat);
[peakValRun,  peakIdxRun]  = max(meanRun);
fprintf('\nStationary: peak %.3f at t=%.2fs\n', peakValStat, timeVecRef(peakIdxStat));
fprintf('Running:    peak %.3f at t=%.2fs\n', peakValRun, timeVecRef(peakIdxRun));

% time to decay back to 50% of peak-above-baseline, post-peak
baselineStat = mean(meanStat(timeVecRef < 0), 'omitnan');
baselineRun  = mean(meanRun(timeVecRef < 0), 'omitnan');
halfStat = baselineStat + 0.5*(peakValStat - baselineStat);
halfRun  = baselineRun  + 0.5*(peakValRun  - baselineRun);

postPeakIdxStat = peakIdxStat:numel(timeVecRef);
decayIdxStat = postPeakIdxStat(find(meanStat(postPeakIdxStat) < halfStat, 1));
postPeakIdxRun = peakIdxRun:numel(timeVecRef);
decayIdxRun = postPeakIdxRun(find(meanRun(postPeakIdxRun) < halfRun, 1));

if ~isempty(decayIdxStat)
    fprintf('Stationary: decays to half-peak at t=%.2fs\n', timeVecRef(decayIdxStat));
end
if ~isempty(decayIdxRun)
    fprintf('Running:    decays to half-peak at t=%.2fs\n', timeVecRef(decayIdxRun));
end

%% data-driven onset latency: first post-stimulus time point where the
% population mean trace exceeds baseline + 2*SD(baseline), and STAYS
% above that threshold for at least minSustainSec (avoids flagging a
% single noisy sample as "onset").
minSustainSec = 0.1; % require the crossing to hold for this long
dt = mean(diff(timeVecRef));
minSustainSamples = max(1, round(minSustainSec / dt));

baselineMaskStat = timeVecRef < 0;
baselineSDStat = std(meanStat(baselineMaskStat), 'omitnan');
threshStat = baselineStat + 2*baselineSDStat;

baselineMaskRun = timeVecRef < 0;
baselineSDRun = std(meanRun(baselineMaskRun), 'omitnan');
threshRun = baselineRun + 2*baselineSDRun;

postStimIdx = find(timeVecRef >= 0);

aboveStat = meanStat(postStimIdx) > threshStat;
onsetIdxStat = [];
for k = 1:(numel(aboveStat) - minSustainSamples + 1)
    if all(aboveStat(k:k+minSustainSamples-1))
        onsetIdxStat = postStimIdx(k);
        break;
    end
end

aboveRun = meanRun(postStimIdx) > threshRun;
onsetIdxRun = [];
for k = 1:(numel(aboveRun) - minSustainSamples + 1)
    if all(aboveRun(k:k+minSustainSamples-1))
        onsetIdxRun = postStimIdx(k);
        break;
    end
end

fprintf('\n--- Data-driven onset latency (first sustained crossing of baseline + 2SD) ---\n');
if ~isempty(onsetIdxStat)
    fprintf('Stationary: onset at t=%.2fs (baseline=%.3f, threshold=%.3f)\n', ...
        timeVecRef(onsetIdxStat), baselineStat, threshStat);
else
    fprintf('Stationary: no sustained threshold crossing found -- inspect trace manually.\n');
end
if ~isempty(onsetIdxRun)
    fprintf('Running:    onset at t=%.2fs (baseline=%.3f, threshold=%.3f)\n', ...
        timeVecRef(onsetIdxRun), baselineRun, threshRun);
else
    fprintf('Running:    no sustained threshold crossing found -- inspect trace manually.\n');
end

if ~isempty(onsetIdxStat) && ~isempty(onsetIdxRun)
    earliestOnset = min(timeVecRef(onsetIdxStat), timeVecRef(onsetIdxRun));
    fprintf('\nSuggested window start (earlier of the two onsets): t=%.2fs\n', earliestOnset);
    fprintf('Suggested window end (past both half-decay points): t=%.2fs\n', ...
        max(timeVecRef(decayIdxStat), timeVecRef(decayIdxRun)));
end

%% mark onset/peak/half-decay on the PSTH plot for visual confirmation
figure(popPSTHFig);
hold on;
if ~isempty(onsetIdxStat), xline(timeVecRef(onsetIdxStat), 'k:', 'Stat onset', 'HandleVisibility','off'); end
if ~isempty(onsetIdxRun),  xline(timeVecRef(onsetIdxRun),  'r:', 'Run onset',  'HandleVisibility','off'); end
