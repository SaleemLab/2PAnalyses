%% Plot ALL responsive boutons in tiled grid pages (for manual inspection/note-taking)
% Requires: allRFMapping, candidateRespIdxList, uAz, uEl_plot, timeVector,
%           plotRFHeatmapWithTraces.m already on path (as built in your main script).
% Each panel is labeled with its bouton ID (index into allRFMapping), so you
% can note down interesting bouton numbers directly while looking through.

nPerFig = 12;   % boutons per figure page
nCols   = 4;    % grid columns per page (nPerFig should be divisible by nCols)
nRows   = ceil(nPerFig / nCols);

boutonsToPlot = candidateRespIdxList;   % change to respIdxList / gaussianRespIdxList if desired
nTotal = numel(boutonsToPlot);
nPages = ceil(nTotal / nPerFig);

fprintf('Plotting %d responsive boutons across %d page(s), %d per page.\n', ...
    nTotal, nPages, nPerFig);

for pg = 1:nPages
    idxStart = (pg-1)*nPerFig + 1;
    idxEnd   = min(pg*nPerFig, nTotal);
    thisPageIdx = boutonsToPlot(idxStart:idxEnd);
    nThisPage = numel(thisPageIdx);

    figPage = figure('Color', 'w', 'Position', [50 50 nCols*350 nRows*320], ...
        'Name', sprintf('Responsive boutons, page %d/%d', pg, nPages));

    for i = 1:nThisPage
        boutonID = thisPageIdx(i);
        b = allRFMapping(boutonID);

        r = ceil(i / nCols);
        c = mod(i-1, nCols) + 1;

        left   = (c-1)/nCols + 0.02;
        bottom = 1 - r/nRows + 0.06;
        width  = 1/nCols - 0.04;
        height = 1/nRows - 0.12;

        axPanel = axes('Position', [left bottom width height]);
        plotRFHeatmapWithTraces(axPanel, b, uAz, uEl_plot, timeVector, 'Smooth', true);
        title(axPanel, sprintf('Bouton %d', boutonID), ...
            'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    end
end
