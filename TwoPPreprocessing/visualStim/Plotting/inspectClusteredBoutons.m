%% Inspect boutons sharing the same preferred position, within one session
% Picks the session with the most responsive boutons, finds its most
% "crowded" Az/El position, and plots every bouton at that position side
% by side so you can visually check whether they look like genuinely
% independent RFs (different noise, slightly different shapes) or
% suspiciously identical (a sign of duplicate ROIs / bleed-through /
% shared motion artifact rather than real independent tuning).

targetSession = '';  % 'M25132_20260306'

if isempty(targetSession)
    sessCountsAll = countcats(categorical(respSessionLabels));
    sessNamesAll  = categories(categorical(respSessionLabels));
    [~, maxIdx] = max(sessCountsAll);
    targetSession = sessNamesAll{maxIdx};
end

sessMask = strcmp(respSessionLabels, targetSession);
sessBoutonGlobalIdx = candidateRespIdxList(sessMask);  % indices into allRFMapping
sessAz = prefAz(sessMask);
sessEl = prefEl(sessMask);

[uniquePos, ~, idx] = unique([sessAz(:), sessEl(:)], 'rows');
counts = accumarray(idx, 1);
[~, mostCrowdedRank] = max(counts);
targetPos = uniquePos(mostCrowdedRank, :);

boutonsAtPos = sessBoutonGlobalIdx(idx == mostCrowdedRank);

fprintf('Session: %s\n', targetSession);
fprintf('Most crowded position: Az=%.0f, El=%.0f (%d boutons)\n', ...
    targetPos(1), targetPos(2), numel(boutonsAtPos));
fprintf('Bouton IDs: %s\n', mat2str(boutonsAtPos));

%% Plot all boutons at this position, side by side
nB = numel(boutonsAtPos);
figInspect = figure('Color', 'w', 'Position', [50 50 nB*400 450], ...
    'Name', sprintf('%s -- boutons at Az=%.0f,El=%.0f', targetSession, targetPos(1), targetPos(2)));

colW = 1/nB;
for i = 1:nB
    b = allRFMapping(boutonsAtPos(i));
    xBase = (i-1)*colW;
    axPanel = axes('Position', [xBase+0.03 0.15 colW-0.05 0.70]);
    plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
    title(axPanel, sprintf('Bouton %d', boutonsAtPos(i)), ...
        'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
end

%% Quick check: are these boutons' traces suspiciously correlated?
% (a sign of shared contamination rather than independent real tuning)
fprintf('\nPairwise correlation of mean grid responses (flattened), boutons at this position:\n');
for i = 1:nB
    for j = i+1:nB
        gi = allRFMapping(boutonsAtPos(i)).meanGridResponse;
        gj = allRFMapping(boutonsAtPos(j)).meanGridResponse;
        r = corr(gi(:), gj(:), 'rows', 'complete');
        fprintf('  Bouton %d vs %d: r = %.2f\n', boutonsAtPos(i), boutonsAtPos(j), r);
    end
end
fprintf('(Very high r, e.g. >0.95, across many pairs would suggest shared\n');
fprintf(' contamination rather than independently real tuning at this position.)\n');
