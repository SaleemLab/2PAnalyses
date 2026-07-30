function hFig=plotExampleDirTuningBouton(allDirTuning, boutonIdx, hFig, smoothWindow)
% PLOTEXAMPLEDIRTUNINGBOUTON  Reproduce a reference-figure-style summary
% for one bouton: polar tuning plot + text (left) and mini per-direction
% PSTH traces (individual trials + mean) with directional arrows (bottom
% row), annotated with preferred direction, OSI, and DSI. Landscape
% layout.
%
%   plotExampleDirTuningBouton(allDirTuning, boutonIdx)
%   plotExampleDirTuningBouton(allDirTuning, boutonIdx, hFig)
%   plotExampleDirTuningBouton(allDirTuning, boutonIdx, hFig, smoothWindow)
%
% If hFig is provided, draws into that existing figure (cleared first)
% instead of creating a new one -- used for batch PDF export across many
% boutons without popping open a new visible window each time.
%
% smoothWindow (optional, default 1 = no smoothing): window size in
% FRAMES for a moving-average smooth (movmean) applied to each trial's
% trace before plotting/averaging. DISPLAY ONLY -- this never touches
% meanDirResponse, trialRawResp, OSI/DSI, or any other stored statistic;
% it's purely cosmetic smoothing to make the response shape easier to
% see by eye when trial counts are low and noisy (e.g. the running-only
% cohort). Keep it modest (e.g. 3-5 frames) so you're smoothing out
% frame-to-frame sampling noise, not smearing over real response
% kinetics -- check your frame rate (1/diff(timeVec(1:2))) to judge
% what's reasonable.

if nargin < 3 || isempty(hFig)
    hFig = figure('Color', 'w', 'Position', [100 100 1000 550]);
else
    clf(hFig);
end
if nargin < 4 || isempty(smoothWindow)
    smoothWindow = 1; % no smoothing by default
end

s = allDirTuning(boutonIdx);

if ~isfield(s, 'OSI_simple') || ~isfield(s, 'DSI_simple')
    error('Bouton %d is missing OSI/DSI fields -- run computeDirTuningOSI and computeDirTuningDSI first.', boutonIdx);
end
if ~isfield(s, 'fullTraceSub') || ~isfield(s, 'timeVec')
    error(['Bouton %d is missing fullTraceSub/timeVec -- re-run the pooling script ' ...
        '(updated version that stores the full trace) before calling this function.'], boutonIdx);
end

thetaDeg = s.stimValues(:)';
[thetaSorted, sortIdx] = sort(mod(round(thetaDeg), 360));
nDir = numel(thetaSorted);
R = s.meanDirResponse(sortIdx);
R = R(:)';
R_plot = max(R, 0);

% arrow glyphs for standard math-convention angles (0 = rightward/east,
% increasing counterclockwise). Adjust mapping here if your rig uses a
% different convention for 0 degrees / rotation direction.
arrowMap = containers.Map( ...
    {0, 45, 90, 135, 180, 225, 270, 315}, ...
    {char(8594), char(8599), char(8593), char(8598), char(8592), char(8601), char(8595), char(8600)});

% arrow label for each sorted direction, used both under the mini panels
% and as the polar-plot tick labels
arrowLabels = cell(1, nDir);
for dIdx = 1:nDir
    if isKey(arrowMap, thetaSorted(dIdx))
        arrowLabels{dIdx} = arrowMap(thetaSorted(dIdx));
    else
        arrowLabels{dIdx} = '?';
    end
end

figure(hFig);

%%  left: polar tuning plot
axPolar = polaraxes('Position', [0.04 0.30 0.24 0.62]); % [left bottom width height], normalized
thetaRadClosed = deg2rad([thetaSorted, thetaSorted(1)]);
Rclosed = [R_plot, R_plot(1)];
polarplot(axPolar, thetaRadClosed, Rclosed, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
axPolar.ThetaZeroLocation = 'right';
axPolar.ThetaDir = 'counterclockwise';
% tick labels at the actual stimulus angles (from stimValues), in
% degrees -- arrows are used only for the mini-panel labels below
axPolar.ThetaTick = thetaSorted;
axPolar.ThetaTickLabel = arrayfun(@(x) sprintf('%d%s', x, char(176)), thetaSorted, 'UniformOutput', false);
title(axPolar, sprintf('Bouton %d', boutonIdx));

%%  text annotations (right of polar plot)
annotStr = sprintf(['Preferred direction: %.0f%s\n' ...
                     'Orientation sel. index: %.3f\n' ...
                     'Direction sel. index : %.3f'], ...
    s.prefDirectionDeg_simple, char(176), s.OSI_simple, s.DSI_simple); % simple-ratio preferred direction (vector method is commented out in computeDirTuningDSI.m)
annotation(hFig, 'textbox', [0.30 0.55 0.28 0.35], 'String', annotStr, ...
    'EdgeColor', 'none', 'FontSize', 9, 'VerticalAlignment', 'top');

%% bottom: mini PSTH row (fixed display window, spans full width)
timeVec = s.timeVec(:)';
stimOnDuration = s.stimOnDuration;

% fixed display window for all mini panels -- previously the full
% pre+post trace was shown, which could run well past the interesting
% part of the response; clip display (not the underlying data) to
% [-0.5, 4] s.
displayWindow = [-0.5, 4];
displayMask   = timeVec >= displayWindow(1) & timeVec <= displayWindow(2);
timeVecDisp   = timeVec(displayMask);

% compute yLims ONCE, before the per-direction panel loop, from the raw
% (smoothed-for-display) individual trial values across ALL directions
% -- since the mini panels show the raw trial overlay (not a mean+/-SEM
% band), yLims needs to span the actual trial-to-trial range, not a
% narrower mean+/-SEM range, or trials get clipped at the panel edges.
allBounds = [];
for dIdx = 1:nDir
    tm = s.fullTraceSub{sortIdx(dIdx)};
    tm = tm(:, displayMask);
    if isempty(tm), continue; end
    if smoothWindow > 1
        tm = movmean(tm, smoothWindow, 2, 'omitnan'); % display-only smoothing, per trial
    end
    allBounds = [allBounds, tm(:)']; %#ok<AGROW> % raw trial values, not mean+/-SEM
end
yLims = [min(allBounds, [], 'omitnan'), max(allBounds, [], 'omitnan')];
if any(isnan(yLims)) || diff(yLims) == 0
    yLims = [-0.1 0.1];
end

panelLeft   = 0.04;
panelWidth  = 0.90 / nDir;
panelBottom = 0.12;
panelHeight = 0.28;

% stimulus-on indicator: a short bar near the bottom of each panel
% instead of a full-height shaded box, which was visually squashing the
% trace against the top/bottom of the panel
stimBarHeight = yLims(1) + 0.08 * diff(yLims);

for d = 1:nDir
    axP = axes('Position', [panelLeft + (d-1)*panelWidth, panelBottom, panelWidth*0.9, panelHeight]);
    hold(axP, 'on');

    % short stimulus-on indicator bar, from t=0 to t=stimOnDuration
    if ~isnan(stimOnDuration)
        patch(axP, [0 stimOnDuration stimOnDuration 0], ...
            [yLims(1) yLims(1) stimBarHeight stimBarHeight], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
    end

    traceMat = s.fullTraceSub{sortIdx(d)}; % [nTrials x nFrames]
    traceMat = traceMat(:, displayMask);
    if smoothWindow > 1
        traceMat = movmean(traceMat, smoothWindow, 2, 'omitnan'); % display-only smoothing, per trial -- same window used for yLims above
    end

    % individual-trial overlay + mean
    if ~isempty(traceMat)
        hTrials = plot(axP, timeVecDisp, traceMat', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
        for k = 1:numel(hTrials)
            hTrials(k).Color(4) = 0.3; % alpha: 0 = fully transparent, 1 = fully opaque
        end
        plot(axP, timeVecDisp, mean(traceMat, 1, 'omitnan'), 'k-', 'LineWidth', 1.5);
    end

    xline(axP, 0, 'k:', 'LineWidth', 0.75); % stimulus onset marker

    ylim(axP, yLims);
    xlim(axP, displayWindow);

    if d == 1
        % show a real y-axis on the first panel instead of relying only
        % on the separate scale-bar axis
        axP.XColor = 'none';   % hide x-axis line/ticks, keep y-axis
        axP.YColor = 'k';
        axP.Box    = 'off';
        axP.TickDir = 'out';
        axP.FontSize = 7;
        ylabel(axP, '\DeltaF/F', 'FontSize', 8);
    else
        axis(axP, 'off');
    end

    % arrow label below each mini-panel
    text(axP, 0.5, -0.18, arrowLabels{d}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
end

% scale bar to the right of the last panel -- height picked
% automatically as a "nice" round number close to a target fraction of
% the actual y-range, since yLims is data-derived and a hardcoded height
% could be wildly mismatched to the real scale of a given bouton's data
targetFrac = 0.4; % scale bar spans roughly this fraction of the panel's y-range
rawScaleBarHeight = targetFrac * diff(yLims);
niceSteps = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10];
[~, niceIdx] = min(abs(niceSteps - rawScaleBarHeight));
scaleBarHeight = niceSteps(niceIdx);

axScale = axes('Position', [panelLeft + nDir*panelWidth, panelBottom, 0.05, panelHeight], 'Visible', 'off');
xlim(axScale, [0 1]); ylim(axScale, yLims);
line(axScale, [0.2 0.2], [yLims(1) yLims(1)+scaleBarHeight], 'Color', 'k', 'LineWidth', 1.5); % anchored to the bottom of the panel, not a fixed y=0
text(axScale, 0.3, yLims(1) + scaleBarHeight/2, sprintf('%.3g \\DeltaF/F', scaleBarHeight), ...
    'FontSize', 8, 'Units', 'data');

[~, prefIdxSimple] = max(s.meanDirResponse);
prefDirSimpleDeg = s.stimValues(prefIdxSimple);

sgtitle(sprintf('Bouton %d | pref dir=%.0f%s | OSI=%.3f | DSI=%.3f', ...
    boutonIdx, prefDirSimpleDeg, char(176), s.OSI_simple, s.DSI_simple));

end