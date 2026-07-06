function plotExampleDirTuningBouton(allDirTuning, boutonIdx, hFig)
% PLOTEXAMPLEDIRTUNINGBOUTON  Reproduce a reference-figure-style summary
% for one bouton: polar tuning plot + text (left) and mini per-direction
% PSTH traces with directional arrows (bottom row), annotated with
% preferred direction, OSI, and DSI. Landscape layout.
%
%   plotExampleDirTuningBouton(allDirTuning, boutonIdx)
%   plotExampleDirTuningBouton(allDirTuning, boutonIdx, hFig)
%
% If hFig is provided, draws into that existing figure (cleared first)
% instead of creating a new one -- used for batch PDF export across many
% boutons without popping open a new visible window each time.
%
% REQUIRES: allDirTuning(boutonIdx) must already have OSI/DSI/
% prefDirectionDeg computed (run computeDirTuningOSI and
% computeDirTuningDSI first), and must have been pooled with the
% updated pooling script that stores fullTraceSub/timeVec/stimOnDuration
% (full pre+post trace per trial, not just the response-window slice).

if nargin < 3 || isempty(hFig)
    hFig = figure('Color', 'w', 'Position', [100 100 1000 550]);
else
    clf(hFig);
end

s = allDirTuning(boutonIdx);

if ~isfield(s, 'OSI') || ~isfield(s, 'DSI')
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
R = R(:)'; % force row vector -- meanDirResponse may be stored as a column, causing a shape mismatch below otherwise

% arrow glyphs for standard math-convention angles (0 = rightward/east,
% increasing counterclockwise). Adjust mapping here if your rig uses a
% different convention for 0 degrees / rotation direction.
arrowMap = containers.Map( ...
    {0, 45, 90, 135, 180, 225, 270, 315}, ...
    {char(8594), char(8599), char(8593), char(8598), char(8592), char(8601), char(8595), char(8600)});

figure(hFig);

%%  left: polar tuning plot 
axPolar = polaraxes('Position', [0.04 0.30 0.24 0.62]); % [left bottom width height], normalized
thetaRadClosed = deg2rad([thetaSorted, thetaSorted(1)]);
Rclosed = [R, R(1)];
polarplot(axPolar, thetaRadClosed, Rclosed, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
axPolar.ThetaZeroLocation = 'right';
axPolar.ThetaDir = 'counterclockwise';
title(axPolar, sprintf('Bouton %d', boutonIdx));

%%  text annotations (right of polar plot)
annotStr = sprintf(['Preferred direction: %.0f%s\n' ...
                     'Orientation sel. index: %.3f\n' ...
                     'Direction sel. index : %.3f'], ...
    s.prefDirectionDeg, char(176), s.OSI_simple, s.DSI_simple);
annotation(hFig, 'textbox', [0.30 0.55 0.28 0.35], 'String', annotStr, ...
    'EdgeColor', 'none', 'FontSize', 9, 'VerticalAlignment', 'top');

%% bottom: mini PSTH row (full timecourse, spans full width)
timeVec = s.timeVec(:)';
stimOnDuration = s.stimOnDuration;

allTraces = [];
for d = 1:nDir
    allTraces = [allTraces; s.fullTraceSub{sortIdx(d)}(:)]; 
end
yLims = [min(allTraces, [], 'omitnan'), max(allTraces, [], 'omitnan')];
if any(isnan(yLims)) || diff(yLims) == 0
    yLims = [-0.1 0.1];
end

panelLeft   = 0.04;
panelWidth  = 0.90 / nDir;
panelBottom = 0.12;
panelHeight = 0.28;

for d = 1:nDir
    axP = axes('Position', [panelLeft + (d-1)*panelWidth, panelBottom, panelWidth*0.9, panelHeight]);
    hold(axP, 'on');

    % shaded stimulus-on box, from t=0 to t=stimOnDuration
    if ~isnan(stimOnDuration)
        patch(axP, [0 stimOnDuration stimOnDuration 0], ...
            [yLims(1) yLims(1) yLims(2) yLims(2)], ...
            [0.9 0.9 0.9], 'EdgeColor', 'none');
    end

    traceMat = s.fullTraceSub{sortIdx(d)}; % [nTrials x nFrames]
    if ~isempty(traceMat)
        plot(axP, timeVec, traceMat', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
        plot(axP, timeVec, mean(traceMat, 1, 'omitnan'), 'k-', 'LineWidth', 1.5);
    end

    xline(axP, 0, 'k:', 'LineWidth', 0.75); % stimulus onset marker

    ylim(axP, yLims);
    xlim(axP, [timeVec(1), timeVec(end)]);
    axis(axP, 'off');

    % arrow label below each mini-panel
    thisAngle = thetaSorted(d);
    if isKey(arrowMap, thisAngle)
        arrowChar = arrowMap(thisAngle);
    else
        arrowChar = '?';
    end
    text(axP, 0.5, -0.18, arrowChar, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
end

% scale bar to the right of the last panel
scaleBarHeight = 0.5; %
axScale = axes('Position', [panelLeft + nDir*panelWidth, panelBottom, 0.05, panelHeight], 'Visible', 'off');
xlim(axScale, [0 1]); ylim(axScale, yLims);
line(axScale, [0.2 0.2], [0 scaleBarHeight], 'Color', 'k', 'LineWidth', 1.5);
text(axScale, 0.3, scaleBarHeight/2, sprintf('%.2g \\DeltaF/F', scaleBarHeight), ...
    'FontSize', 8, 'Units', 'data');

[~, prefIdxSimple] = max(s.meanDirResponse);
prefDirSimpleDeg = s.stimValues(prefIdxSimple);

sgtitle(sprintf('Bouton %d | pref dir=%.0f%s | OSI=%.3f | DSI=%.3f', ...
    boutonIdx, prefDirSimpleDeg, char(176), s.OSI_simple, s.DSI_simple));

end
