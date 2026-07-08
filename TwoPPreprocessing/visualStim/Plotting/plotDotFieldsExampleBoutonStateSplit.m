function plotDotFieldsExampleBoutonStateSplit(thisMouse, thisSessionName, thisROI)

%%
stimFramesMask_range = [0 2.0];
runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;
respWin_range   = [0 2.5]; % window for tuning-curve/response-magnitude calculation

%% load raw session data 
infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
if ~isfile(infoPath), error('sessionFileInfo not found for %s %s', thisMouse, thisSessionName); end
loadedInfo      = load(infoPath, 'sessionFileInfo');
sessionFileInfo = loadedInfo.sessionFileInfo;
stimNames       = {sessionFileInfo.stimFiles.name};
dotIdx = find(contains(stimNames, 'DotMotion_SpeedTuning'), 1);
if isempty(dotIdx), error('No DotMotion_SpeedTuning file found for %s %s', thisMouse, thisSessionName); end

load(sessionFileInfo.stimFiles(dotIdx).Response, 'response');
load(sessionFileInfo.stimFiles(dotIdx).BonsaiData, 'bonsaiData');

%% rebuild trial classification 
nGroups = numel(response.wheelData);
trialsSpeed2D = struct('VelX1', {}, 'numDots1', {}, 'runFlag', {}, 'origGroup', {}, 'origTrialInGroup', {});
trialCounter = 1;

for g = 1:nGroups
    grpWheel  = response.wheelData(g);
    speedMatrix = grpWheel.alignedResponses;
    tVecWheel   = grpWheel.timeVector;
    stimFramesMask = (tVecWheel >= stimFramesMask_range(1) & tVecWheel <= stimFramesMask_range(2));

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
        trialsSpeed2D(trialCounter).origTrialInGroup  = ti;
        trialCounter = trialCounter + 1;
    end
end

tsd = trialsSpeed2D;
temp_tsd = tsd([tsd.numDots1] == 573); % non-blank only, for this figure

uniqueVelocities = unique(abs([temp_tsd.VelX1]));
nSpeeds = numel(uniqueVelocities);

timeVec = response.psthData(1).timeVector(:)';
respIdx = timeVec >= respWin_range(1) & timeVec <= respWin_range(2);

stateNames  = {'Stationary', 'Running'};
stateColors = {'k', 'r'};

%%  extract full timecourses + tuning values, per speed x state 
fullTraces = cell(nSpeeds, 2); % {speed, state} -> [nTrials x nTimepoints]
tuningVals = cell(nSpeeds, 2); % {speed, state} -> [nTrials x 1], mean within respWin

for s = 1:nSpeeds
    for istate = 1:2
        matchingTrials = find(abs([temp_tsd.VelX1]) == uniqueVelocities(s) & [temp_tsd.runFlag] == (istate - 1));
        traceMat = nan(numel(matchingTrials), numel(timeVec));
        for mt = 1:numel(matchingTrials)
            origGroup = temp_tsd(matchingTrials(mt)).origGroup;
            origTi    = temp_tsd(matchingTrials(mt)).origTrialInGroup;
            traceMat(mt, :) = squeeze(response.psthData(origGroup).alignedResponses(thisROI, :, origTi));
        end
        fullTraces{s, istate} = traceMat;
        tuningVals{s, istate} = mean(traceMat(:, respIdx), 2, 'omitnan');
    end
end

%%  Fig
figure('Color', 'w', 'Position', [50 50 1400 900]);
sgtitle(sprintf('%s | %s | ROI %d', thisMouse, thisSessionName, thisROI), 'Interpreter', 'none', 'FontWeight', 'bold');
speedCmap = parula(nSpeeds);

%%  Row 1: PSTH, one subplot per state 
for istate = 1:2
    subplot(3, 2, istate); hold on;
    for s = 1:nSpeeds
        traceMat = fullTraces{s, istate};
        if isempty(traceMat), continue; end
        mTrace = mean(traceMat, 1, 'omitnan');
        semTrace = std(traceMat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(traceMat), 1));
        fill([timeVec, fliplr(timeVec)], [mTrace+semTrace, fliplr(mTrace-semTrace)], ...
            speedCmap(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(timeVec, mTrace, 'Color', speedCmap(s,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%.0f', uniqueVelocities(s)));
    end
    xline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('Time (s)'); ylabel('\DeltaF/F');
    title(sprintf('PSTH: %s', stateNames{istate}));
    if istate == 1, legend('Location', 'best', 'FontSize', 7); end
end

%%  Row 2: tuning curve, both states overlaid 
subplot(3, 2, [3 4]); hold on;
for istate = 1:2
    meanVals = nan(nSpeeds, 1); semVals = nan(nSpeeds, 1);
    for s = 1:nSpeeds
        vals = tuningVals{s, istate};
        meanVals(s) = mean(vals, 'omitnan');
        semVals(s)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end
    errorbar(1:nSpeeds, meanVals, semVals, 'o-', 'Color', stateColors{istate}, ...
        'LineWidth', 1.5, 'MarkerFaceColor', stateColors{istate}, 'DisplayName', stateNames{istate});
end
set(gca, 'XTick', 1:nSpeeds, 'XTickLabel', arrayfun(@(v) sprintf('%.0f', v), uniqueVelocities, 'UniformOutput', false));
xlabel('Visual speed (deg/s)'); ylabel('Avg \DeltaF/F');
title('Speed tuning curve'); legend('Location', 'best');

%%  Row 3: per-speed heatmaps, stationary block on top, running block below 
for s = 1:nSpeeds
    subplot(3, nSpeeds, 2*nSpeeds + s);

    statMat = fullTraces{s, 1};
    runMat  = fullTraces{s, 2};
    combinedMat = [statMat; runMat];
    nStatTrials = size(statMat, 1);

    imagesc(timeVec, 1:size(combinedMat,1), combinedMat);
    hold on;
    xline(0, 'w:', 'LineWidth', 1);
    if nStatTrials > 0 && nStatTrials < size(combinedMat,1)
        yline(nStatTrials + 0.5, 'w-', 'LineWidth', 2); % divider between stat/run blocks
    end
    xlabel('Time (s)');
    if s == 1, ylabel('Trials (stat top, run bottom)'); end
    title(sprintf('Sp: %.0f', uniqueVelocities(s)), 'FontSize', 9);
    colormap(gca, 'parula');
end

end