function plotDotFieldsExampleBouton(allDotUnits, mouseID, sessionName, roiIdx, respWin_range)
% plotDotFieldsExampleBouton(allDotUnits, mouseID, sessionName, roiIdx, respWin_range)
%
% Panel A (top-left):   PSTH, POOLED ACROSS BEHAVIORAL STATE (not split),
%                        one line per speed, color-coded, shaded SEM.
% Panel B (top-right):  Speed tuning curve, SPLIT BY STATE (black =
%                        Stationary, red = Running), overlaid.
% Panel C (bottom row):  Per-speed trial rasters, manually positioned
%                        axes. Trials sorted stat-block-top /
%                        run-block-bottom, divider line, shared
%                        quantile-based color scale + one shared colorbar.
%
% USAGE: pass mouseID, sessionName, and roiIdx directly -- e.g. matching
% what you see printed on an existing figure title ("M26003 | 20260326 |
% ROI 20"):
%
%   plotDotFieldsExampleBouton(allDotUnits, 'M26003', '20260326', 20)
%
% The (mouseID, sessionName, roiIdx) -> allDotUnits index lookup is a
% pure in-memory search over fields already stored in allDotUnits --
% NO session/disk loading is needed for this step. Disk is only touched
% afterwards, to pull the full per-trial traces for plotting (and only
% once per session, thanks to internal caching -- so plotting several
% ROIs from the same session stays fast after the first call).
%
% respWin_range optional, defaults to [0.1 3].

persistent cachedSessionLabel cachedResponse cachedTemp_tsd cachedUniqueVelocities cachedNSpeeds cachedTimeVec

if nargin < 5 || isempty(respWin_range), respWin_range = [0.1 3]; end

%% look up the allDotUnits index from mouseID/sessionName/roiIdx -- no disk access
matchIdx = find(strcmp({allDotUnits.mouseID}, mouseID) & ...
                strcmp({allDotUnits.sessionName}, sessionName) & ...
                [allDotUnits.roiIdx] == roiIdx);

if isempty(matchIdx)
    error('No bouton found in allDotUnits matching mouseID=%s, sessionName=%s, roiIdx=%d.', ...
        mouseID, sessionName, roiIdx);
elseif numel(matchIdx) > 1
    warning('Multiple matches found (%d) -- using the first one.', numel(matchIdx));
    matchIdx = matchIdx(1);
end
boutonIdx = matchIdx;
unit = allDotUnits(boutonIdx);
thisMouse       = unit.mouseID;
thisSessionName = unit.sessionName;
thisROI         = unit.roiIdx;
sessionLabel    = unit.sessionLabel;

%% classification thresholds (same as main pipeline)
stimFramesMask_range = [0 2.0];
runSpeedThresh  = 3;
statSpeedThresh = 0.5;
propThresh      = 0.75;


xWindow = [-0.1 3];
%% load or reuse cached session data
if isequal(cachedSessionLabel, sessionLabel)
    response         = cachedResponse;
    temp_tsd         = cachedTemp_tsd;
    uniqueVelocities = cachedUniqueVelocities;
    nSpeeds          = cachedNSpeeds;
    timeVec          = cachedTimeVec;
else
    fprintf('Loading session %s from disk (not cached)...\n', sessionLabel);
    infoPath = findSessionFileInfoFilePath(thisMouse, thisSessionName);
    if ~isfile(infoPath), error('sessionFileInfo not found for %s %s', thisMouse, thisSessionName); end
    loadedInfo      = load(infoPath, 'sessionFileInfo');
    sessionFileInfo = loadedInfo.sessionFileInfo;
    stimNames       = {sessionFileInfo.stimFiles.name};
    dotIdx = find(contains(stimNames, 'DotMotion_SpeedTuning'), 1);
    if isempty(dotIdx), error('No DotMotion_SpeedTuning file found for %s %s', thisMouse, thisSessionName); end

    load(sessionFileInfo.stimFiles(dotIdx).Response, 'response');
    load(sessionFileInfo.stimFiles(dotIdx).BonsaiData, 'bonsaiData'); %#ok<NASGU>

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
    temp_tsd = tsd([tsd.numDots1] == 573);
    uniqueVelocities = unique(abs([temp_tsd.VelX1]));
    nSpeeds = numel(uniqueVelocities);
    timeVec = response.psthData(1).timeVector(:)';

    cachedSessionLabel     = sessionLabel;
    cachedResponse         = response;
    cachedTemp_tsd         = temp_tsd;
    cachedUniqueVelocities = uniqueVelocities;
    cachedNSpeeds          = nSpeeds;
    cachedTimeVec          = timeVec;
end

respIdx = timeVec >= respWin_range(1) & timeVec <= respWin_range(2);
speedCmap = parula(nSpeeds);
% hexColors = {'#d53e4f', '#fc8d59', '#fdae61', '#e6f598', '#99d594', '#3288bd'};
% speedCmap = cell2mat(cellfun(@(h) hex2rgb(h), hexColors', 'UniformOutput', false));
 
% if nSpeeds ~= 6, trim or recycle as needed, e.g.:
speedCmap = speedCmap(1:nSpeeds, :);
 
stateColors = {'k', 'r'};


%% extract full traces (per speed x state) AND pooled-across-state traces (per speed)
fullTraces   = cell(nSpeeds, 2);
tuningVals   = cell(nSpeeds, 2);
pooledTraces = cell(nSpeeds, 1);

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
    pooledTraces{s} = [fullTraces{s,1}; fullTraces{s,2}];
end

%% Fig
hFig = figure('Color', 'w', 'Position', [80 50 1400 850]);

statLabel = ''; runLabel = '';
if isfield(unit,'isResponsive_stat'), statLabel = sprintf('Resp(stat)=%d', unit.isResponsive_stat); end
if isfield(unit,'isResponsive_run'),  runLabel  = sprintf('Resp(run)=%d', unit.isResponsive_run); end
sgtitle(sprintf('%s | ROI %d | %s %s', strrep(sessionLabel,'_',' '), thisROI, statLabel, runLabel), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

%% Panel A (top-left): PSTH, POOLED ACROSS STATE, one line per speed
subplot(2, 2, 1); hold on;
pHandles = gobjects(nSpeeds, 1);
yl = ylim;
patch([0 2 2 0], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
uistack(findobj(gca,'Type','patch'), 'bottom');
for s = 1:nSpeeds
    traceMat = pooledTraces{s};
    if isempty(traceMat), continue; end
    mTrace = mean(traceMat, 1, 'omitnan');
    sTrace = std(traceMat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(traceMat),1));
    fill([timeVec(:); flipud(timeVec(:))], [mTrace(:)-sTrace(:); flipud(mTrace(:)+sTrace(:))], ...
        speedCmap(s,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    pHandles(s) = plot(timeVec, mTrace, 'Color', speedCmap(s,:), 'LineWidth', 3);
end
xlim(xWindow)
xlabel('Time (s)'); ylabel('\DeltaF/F');
lgd = legend(pHandles, arrayfun(@(v) sprintf('%.0f', v), uniqueVelocities, 'UniformOutput', false), ...
    'Location', 'bestoutside');
title(lgd, 'Speeds');
title('PSTH (pooled across behavioral state)');
defaultAxesProperties(gca, 1);
    
%% Panel B (top-right): tuning curve, split by state
subplot(2, 2, 2); hold on;
stateNames = {'Stationary', 'Running'};
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
title('Speed tuning curve (split by state)'); legend('Location', 'best');
defaultAxesProperties(gca, 1);
    

%% Panel C (bottom row): per-speed rasters, manually positioned axes, sorted stat-top/run-bottom
startX = 0.05;
widthPerRaster = 0.85 / nSpeeds;
rasterHeight = 0.28;
bottomY = 0.08;


% xWindow = [-0.5, 3.5];


% 'parula'      : MATLAB default, perceptually reasonable, good general choice
% flipud(gray)  : classic grayscale, dark=high activity -- matches Fig 3.5a style
% 'turbo'       : higher contrast than parula, good for subtle differences
% 'bone'        : cooler grayscale variant, softer than plain gray
heatmapColormap = flipud(gray); %flipud(gray);  % <-- change to 'parula', 'turbo', 'bone', etc. to try alternatives

% mask restricting which timepoints count toward the normalization range
% (so activity outside your displayed window doesn't affect the scale)
xWindowMask = timeVec >= xWindow(1) & timeVec <= xWindow(2);

allDataAcrossSpeeds = [];
for s = 1:nSpeeds
    thisData = pooledTraces{s}(:, xWindowMask);
    allDataAcrossSpeeds = [allDataAcrossSpeeds; thisData(:)];
end
allDataAcrossSpeeds = allDataAcrossSpeeds(~isnan(allDataAcrossSpeeds));

% --- NORMALIZATION: min-max rescale this bouton's own 1st-99th
% percentile range (within the displayed window) to [0, 1]
dataMin = quantile(allDataAcrossSpeeds, 0.01);
dataMax = quantile(allDataAcrossSpeeds, 0.99);
if dataMax <= dataMin, dataMax = dataMin + 1e-6; end % guard against degenerate range

globalCLim = [0, 1];
climMid = 0.5;

for s = 1:nSpeeds
    axPos = [startX + (s-1)*widthPerRaster, bottomY, widthPerRaster*0.85, rasterHeight];
    axR = axes(hFig, 'Position', axPos);

    % apply the same min-max normalization to this speed's data before plotting
    combinedMat = (pooledTraces{s} - dataMin) ./ (dataMax - dataMin);
    nStatTrials = size(fullTraces{s,1}, 1);

    if ~isempty(combinedMat)
        imagesc(axR, timeVec, 1:size(combinedMat,1), combinedMat);
        colormap(axR, heatmapColormap);
        hold(axR, 'on');
        xline(axR, 0, 'w', 'LineWidth', 2);
        if nStatTrials > 0 && nStatTrials < size(combinedMat,1)
            yline(axR, nStatTrials + 0.5, 'w-', 'LineWidth', 2);
        end
        set(axR, 'CLim', globalCLim);
        xlim(axR, [-1.5 3]);   %

        set(axR, 'FontSize', 12);
        title(axR, sprintf('Sp: %.0f', uniqueVelocities(s)), 'FontSize', 12);
        if s == 1
            ylabel(axR, 'Trials (stat top, run bottom)', 'FontSize', 12);
        else
            set(axR, 'YTickLabel', []);
        end
        xlabel(axR, 'Time (s)');

        if s == nSpeeds
            cbPos = [startX + s*widthPerRaster, bottomY, 0.012, rasterHeight];
            cb = colorbar(axR, 'Position', cbPos, 'Location', 'eastoutside');
            cb.Ticks = [globalCLim(1), climMid, globalCLim(2)];
            cb.TickLabels = {sprintf('%.2f', globalCLim(1)), sprintf('%.2f', climMid), sprintf('%.2f', globalCLim(2))};
            cb.TickDirection = 'out';
            cb.Box = 'off';
            cb.FontName = 'Arial';
            cb.FontSize = 10;
            cb.Label.String = 'Normalized activity (a.u.)';
            cb.Label.FontName = 'Arial';
            cb.Label.FontSize = 12;
            cb.Label.Rotation = 90;
            cb.Label.VerticalAlignment = 'bottom';
        end
    end
end

baseFileName = sprintf('%s_%s_exampleROI', mouseID, sessionName);
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section3_Fig4_3\exmaple_traces_trials';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = fullfile(outputDir, baseFileName);
saveFigureFormats(hFig, fullSavePath);

end
