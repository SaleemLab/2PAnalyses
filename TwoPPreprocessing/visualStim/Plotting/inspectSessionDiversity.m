%% Check whether within-session position DIVERSITY is real
% Picks one representative bouton from EACH distinct preferred Az/El
% position within a session, plots them side by side in the SAME order as
% their assigned position (sorted by azimuth then elevation) so you can
% directly check: does the RF's actual bright spot visibly move to match
% its assigned preferred position, or do they all look similar/noisy
% regardless of what position got assigned?

% reuse targetSession from inspectClusteredBoutons.m, or set explicitly:
% targetSession = 'M25132_20260306';

sessMask = strcmp(respSessionLabels, targetSession);
sessBoutonGlobalIdx = candidateRespIdxList(sessMask);
sessAz = prefAz(sessMask);
sessEl = prefEl(sessMask);

[uniquePos, ~, idx] = unique([sessAz(:), sessEl(:)], 'rows');
nUniquePos = size(uniquePos, 1);

fprintf('Session: %s has %d responsive boutons across %d distinct positions.\n', ...
    targetSession, numel(sessBoutonGlobalIdx), nUniquePos);

% one representative bouton per position (just the first one found there)
representativeBoutons = nan(nUniquePos, 1);
for p = 1:nUniquePos
    idxAtPos = find(idx == p, 1, 'first');
    representativeBoutons(p) = sessBoutonGlobalIdx(idxAtPos);
end

% sort by position (azimuth, then elevation) so the plotted order visually
% sweeps across the grid, making any real shift easy to see
[sortedPos, sortOrder] = sortrows(uniquePos, [1 2]);
representativeBoutons = representativeBoutons(sortOrder);

fprintf('\nOne representative bouton per position (sorted by Az, then El):\n');
for p = 1:nUniquePos
    fprintf('  Az=%.0f, El=%.0f -> Bouton %d\n', sortedPos(p,1), sortedPos(p,2), representativeBoutons(p));
end

%% Plot all representatives side by side, in position order
nB = nUniquePos;
figDivReal = figure('Color', 'w', 'Position', [50 50 nB*350 450], ...
    'Name', sprintf('%s -- one bouton per distinct position', targetSession));

colW = 1/nB;
for i = 1:nB
    b = allRFMapping(representativeBoutons(i));
    xBase = (i-1)*colW;
    axPanel = axes('Position', [xBase+0.03 0.18 colW-0.05 0.65]);
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    title(axPanel, sprintf('Bouton %d\n(assigned Az=%.0f,El=%.0f)', ...
        representativeBoutons(i), sortedPos(i,1), sortedPos(i,2)), ...
        'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
end
