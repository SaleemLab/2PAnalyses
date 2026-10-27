%% Panels e/f style: trusted Gaussian centres (filled) with excluded peak positions (outline)
load('Z:\ibn-vision\USERS\Sonali\Analysis\RFMapping\RF_rsp_boutons.mat');

%%

G = gaussFitResults_rsp;

x0 = -[G.x0];   % positive azimuth convention
y0 =  [G.y0];

relAz  = x0(trustedMask_rsp);   relEl  = y0(trustedMask_rsp);
exclAz = x0(~trustedMask_rsp);  exclEl = y0(~trustedMask_rsp);
fprintf('Reliable fits: %d | Excluded fits: %d\n', numel(relAz), numel(exclAz));

% Edges of the stimulated area
stimAz = [-20 80];
stimEl = [-25 40];


azEdges = -30:5:90;
elEdges = -30:5:50;
% Bins: one per grid position inside the stimulated area,
% plus one bin beyond each edge (up to the fitting limit)
% azEdges = [-30 -20 0 20 40 60 80 90];
% elEdges = [-30 -25 0 20 40 50];

colRelAz = [0.15 0.35 0.35];
colRelEl = [0.8 0.4 0.2];
colExcl  = [0.4 0.4 0.4];

figExt = figure('Color', 'w', 'Position', [100, 100, 650, 280]);

% --- Azimuth ---
ax1 = subplot(1,2,1); hold(ax1, 'on');
histogram(ax1, 'BinEdges', azEdges, 'BinCounts', 100*histcounts(relAz, azEdges)/numel(relAz), ...
    'FaceColor', colRelAz, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('Reliable fit (n=%d)', numel(relAz)));
histogram(ax1, 'BinEdges', azEdges, 'BinCounts', 100*histcounts(exclAz, azEdges)/numel(exclAz), ...
    'DisplayStyle', 'stairs', 'EdgeColor', colExcl, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Excluded fit (n=%d)', numel(exclAz)));
xline(ax1, stimAz(1), '--k', 'HandleVisibility', 'off');
xline(ax1, stimAz(2), '--k', 'HandleVisibility', 'off');
set(ax1, 'XDir', 'reverse', 'XTick', azEdges, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
xlim(ax1, [azEdges(1) azEdges(end)]);
xlabel(ax1, 'Fitted centre azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax1, '% of boutons', 'FontName', 'Arial', 'FontSize', 9);
legend(ax1, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
defaultAxesProperties(ax1, true);

% --- Elevation ---
ax2 = subplot(1,2,2); hold(ax2, 'on');
histogram(ax2, 'BinEdges', elEdges, 'BinCounts', 100*histcounts(relEl, elEdges)/numel(relEl), ...
    'FaceColor', colRelEl, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('Reliable fit (n=%d)', numel(relEl)));
histogram(ax2, 'BinEdges', elEdges, 'BinCounts', 100*histcounts(exclEl, elEdges)/numel(exclEl), ...
    'DisplayStyle', 'stairs', 'EdgeColor', colExcl, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Excluded fit (n=%d)', numel(exclEl)));
xline(ax2, stimEl(1), '--k', 'HandleVisibility', 'off');
xline(ax2, stimEl(2), '--k', 'HandleVisibility', 'off');
set(ax2, 'XTick', elEdges, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
xlim(ax2, [elEdges(1) elEdges(end)]);
xlabel(ax2, 'Fitted centre elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax2, '% of boutons', 'FontName', 'Arial', 'FontSize', 9);
legend(ax2, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
defaultAxesProperties(ax2, true);

set(ax1, 'XDir', 'reverse', 'XTick', -30:10:90, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
set(ax2, 'XTick', -30:10:50, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);


% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\pref_az_el';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figExt, fullfile(outputDir, 'fittedCentre_reliable_vs_excluded_alignedBins'));

%% Counts for the Results
beyondAz  = exclAz > stimAz(2) | exclAz < stimAz(1);
beyondEl  = exclEl > stimEl(2) | exclEl < stimEl(1);
beyondAny = beyondAz | beyondEl;

fprintf('\nExcluded fits with centre beyond the stimulated area:\n');
fprintf('  Azimuth:   %d/%d (%.0f%%)  [peripheral %d, nasal %d]\n', sum(beyondAz), numel(exclAz), 100*mean(beyondAz), ...
    sum(exclAz > stimAz(2)), sum(exclAz < stimAz(1)));
fprintf('  Elevation: %d/%d (%.0f%%)  [top %d, bottom %d]\n', sum(beyondEl), numel(exclEl), 100*mean(beyondEl), ...
    sum(exclEl > stimEl(2)), sum(exclEl < stimEl(1)));
fprintf('  Either:    %d/%d (%.0f%%)\n', sum(beyondAny), numel(exclAz), 100*mean(beyondAny));

% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\pref_az_el';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figExt, fullfile(outputDir, 'fittedCentre_reliable_vs_excluded_extended'));


x0 = -[G.x0];   % positive azimuth convention
y0 =  [G.y0];

nBeyondPeriph = sum(x0 > 80);    % beyond the peripheral edge
nBeyondNasal  = sum(x0 < -20);   % beyond the nasal edge
nBeyondTop    = sum(y0 > 40);    % above the top of the screen
nBeyondBottom = sum(y0 < -25);   % below the bottom of the screen

fprintf('Fitted centres beyond the stimulated area (of %d fits):\n', numel(G));
fprintf('  Peripheral (> 80 deg azimuth): %d\n', nBeyondPeriph);
fprintf('  Nasal (< -20 deg azimuth):     %d\n', nBeyondNasal);
fprintf('  Top (> 40 deg elevation):      %d\n', nBeyondTop);
fprintf('  Bottom (< -25 deg elevation):  %d\n', nBeyondBottom);



%% 
%% Exclusion counts for the RF fits (final criteria: R2 >= 0.1)
load('Z:\ibn-vision\USERS\Sonali\Analysis\RFMapping\RF_rsp_boutons.mat');
G   = gaussFitResults_rsp;
RFM = allRFMapping_rsp;

respROIs = find([RFM.isResponsive]);
nResp    = numel(respROIs);

% Criterion 1: poor fit (R2 < 0.1, so no stored fit)
poorFitROIs = setdiff(respROIs, [G.iROI]);

% Criteria 2-4, from the stored fits
isRobust = false(1, numel(G));
for gi = 1:numel(G)
    vb = G(gi).bootR2(~isnan(G(gi).bootR2));
    if ~isempty(vb), isRobust(gi) = prctile(vb, 5) >= 0.1; end
end
failBoot   = ~isRobust;
failWidth  = [G.isDegenerateSigma];
failCentre = [G.isBeyondRange];

% Sequential: first criterion failed
seqBoot   = failBoot;
seqWidth  = ~failBoot & failWidth;
seqCentre = ~failBoot & ~failWidth & failCentre;
reliable  = ~failBoot & ~failWidth & ~failCentre;

fprintf('\n=== Sequential (first criterion failed), of %d responsive boutons ===\n', nResp);
fprintf('  1. Poor fit (R2 < 0.1):                %3d (%4.1f%%)\n', numel(poorFitROIs), 100*numel(poorFitROIs)/nResp);
fprintf('  2. Unreliable across trials:           %3d (%4.1f%%)\n', sum(seqBoot),   100*sum(seqBoot)/nResp);
fprintf('  3. Not resolvable (width):             %3d (%4.1f%%)\n', sum(seqWidth),  100*sum(seqWidth)/nResp);
fprintf('  4. Centre outside stimulated area:     %3d (%4.1f%%)\n', sum(seqCentre), 100*sum(seqCentre)/nResp);
fprintf('  Reliable fit:                          %3d (%4.1f%%)\n', sum(reliable),  100*sum(reliable)/nResp);
fprintf('  Total:                                 %3d\n', numel(poorFitROIs) + sum(seqBoot) + sum(seqWidth) + sum(seqCentre) + sum(reliable));

fprintf('\n=== Overlapping (fails this criterion, regardless of others), of %d fits ===\n', numel(G));
fprintf('  Unreliable across trials:   %3d\n', sum(failBoot));
fprintf('  Not resolvable (width):     %3d\n', sum(failWidth));
fprintf('  Centre outside:             %3d\n', sum(failCentre));
fprintf('  Fails ONLY unreliable:      %3d\n', sum( failBoot & ~failWidth & ~failCentre));
fprintf('  Fails ONLY width:           %3d\n', sum(~failBoot &  failWidth & ~failCentre));
fprintf('  Fails ONLY centre:          %3d\n', sum(~failBoot & ~failWidth &  failCentre));

% Bouton IDs for picking examples (fails only one criterion, best R2 first)
R2all = [G.R2];
showTop = @(mask, name) fprintf('\n%s (top 10 by R2): %s\n', name, ...
    mat2str(sortByR2([G(mask).iROI], R2all(mask), 10)));
sortByR2 = @(ids, r2, n) ids(subsref(sortrows([r2(:) (1:numel(r2))'], -1), struct('type', '()', 'subs', {{1:min(n, numel(r2)), 2}})));

showTop( failBoot & ~failWidth & ~failCentre, 'Only unreliable');
showTop(~failBoot &  failWidth & ~failCentre, 'Only not resolvable');
showTop(~failBoot & ~failWidth &  failCentre, 'Only centre outside');
showTop(reliable, 'Reliable');
fprintf('\nPoor fit (no stored fit), first 10: %s\n', mat2str(poorFitROIs(1:min(10, numel(poorFitROIs)))));


%% Five examples per exclusion category, one figure per category
load('Z:\ibn-vision\USERS\Sonali\Analysis\RFMapping\RF_rsp_boutons.mat');
G   = gaussFitResults_rsp;
RFM = allRFMapping_rsp;
uAzP = uAz_rsp;  uElP = uEl_rsp;
nPerCat = 5;

gaussFit2D = @(p, xy) p(1) * exp( -( (xy(:,1)-p(2)).^2 ./ (2*p(4)^2) + (xy(:,2)-p(3)).^2 ./ (2*p(5)^2) ) );
[AzGridP, ElGridP] = meshgrid(uAzP, uElP);
xyList  = [AzGridP(:), ElGridP(:)];
respIdx = timeVector >= 0.1 & timeVector <= 3;
azStepP = mean(diff(sort(uAzP)));  elStepP = mean(diff(sort(uElP)));

% --- Categories (fails only one criterion) ---
isRobust = false(1, numel(G));
for gi = 1:numel(G)
    vb = G(gi).bootR2(~isnan(G(gi).bootR2));
    if ~isempty(vb), isRobust(gi) = prctile(vb, 5) >= 0.1; end
end
failBoot = ~isRobust;  failWidth = [G.isDegenerateSigma];  failCentre = [G.isBeyondRange];
R2all = [G.R2];

catMasks = { failBoot & ~failWidth & ~failCentre, ...
            ~failBoot &  failWidth & ~failCentre, ...
            ~failBoot & ~failWidth &  failCentre };
catNames = {'Unreliable across trials', 'Not resolvable', 'Centre outside stimulated area'};
catFiles = {'unreliable', 'notResolvable', 'centreOutside'};

% Pick the 5 best-fitting examples in each category
catGI = cell(1, 3);
for c = 1:3
    idx = find(catMasks{c});
    [~, o] = sort(R2all(idx), 'descend');
    catGI{c} = idx(o(1:min(nPerCat, numel(idx))));
end

% --- Poor fits: refit 5 at random (no stored fit, no bootstrap) ---
respROIs    = find([RFM.isResponsive]);
poorFitROIs = setdiff(respROIs, [G.iROI]);
lb = [0,   min(uAzP)-azStepP, min(uElP)-elStepP, azStepP/4, elStepP/4];
ub = [Inf, max(uAzP)+azStepP, max(uElP)+elStepP, 3*range(uAzP), 3*range(uElP)];
rng(1);
poorPick = poorFitROIs(randperm(numel(poorFitROIs), min(nPerCat, numel(poorFitROIs))));
Gpoor = struct('iROI', {}, 'A', {}, 'x0', {}, 'y0', {}, 'sigmaX', {}, 'sigmaY', {}, 'R2', {}, 'bootR2', {});
for iROI = poorPick
    tm = RFM(iROI).baselineSubtracted;
    rv = cellfun(@(t) mean(mean(double(t(:, respIdx)), 2, 'omitnan'), 'omitnan'), tm(:));
    [~, mi] = max(rv);
    p0 = [max(max(rv), 1e-3), xyList(mi,1), xyList(mi,2), azStepP, elStepP];
    [pF, rn] = lsqcurvefit(gaussFit2D, p0, xyList, rv, lb, ub, optimoptions('lsqcurvefit', 'Display', 'off'));
    Gpoor(end+1) = struct('iROI', iROI, 'A', pF(1), 'x0', pF(2), 'y0', pF(3), 'sigmaX', pF(4), 'sigmaY', pF(5), ...
        'R2', 1 - rn/sum((rv - mean(rv)).^2), 'bootR2', nan(1000,1)); %#ok<SAGROW>
end

% --- All four categories: fits to plot and their names ---
plotSets  = {Gpoor, G(catGI{1}), G(catGI{2}), G(catGI{3})};
plotNames = [{'Poor fit'}, catNames];
plotFiles = [{'poorFit'}, catFiles];

% --- Shared plotting settings ---
limAz = [min(uAzP) - azStepP/2, max(uAzP) + azStepP/2];
limEl = [min(uElP) - elStepP/2, max(uElP) + elStepP/2];
[xq, yq] = meshgrid(linspace(limAz(1), limAz(2), 100), linspace(limEl(1), limEl(2), 100));
flipLabels = @(ax) set(ax, 'XTickLabel', arrayfun(@(v) sprintf('%d', -v), get(ax, 'XTick'), 'UniformOutput', false));

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\excluded_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

for c = 1:numel(plotSets)
    S = plotSets{c};
    nRows = numel(S);
    if nRows == 0, fprintf('No examples for %s\n', plotNames{c}); continue; end
    fig = figure('Color', 'w', 'Position', [100 50 900 200*nRows + 100]);

    leftMargin = 0.08; rightMargin = 0.05; topMargin = 0.06; bottomMargin = 0.06;
    hGap = 0.08; vGap = 0.05;
    colWidth  = (1 - leftMargin - rightMargin - hGap) / 2;
    rowHeight = (1 - topMargin - bottomMargin - (nRows-1)*vGap) / nRows;
    insetFrac = 0.85;
    padW = colWidth * (1 - insetFrac) / 2;  padH = rowHeight * (1 - insetFrac) / 2;

    for r = 1:nRows
        res   = S(r);
        bData = RFM(res.iROI);
        yPos  = 1 - topMargin - r*rowHeight - (r-1)*vGap;

        tm = bData.baselineSubtracted;
        respVec = nan(numel(tm), 1);
        for p = 1:numel(tm)
            if ~isempty(tm{p})
                respVec(p) = mean(mean(double(tm{p}(:, respIdx)), 2, 'omitnan'), 'omitnan');
            end
        end
        boutonData = bData;
        boutonData.meanGridResponse = reshape(respVec, size(AzGridP));
        boutonData.peakAmplitude    = max(respVec) * 2;

        % Left: grid with traces
        axL = axes('Position', [leftMargin+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
        plotRFHeatmapWithTraces(axL, boutonData, uAzP, uElP, timeVector, ...
            'Colormap', 'bone', 'Smooth', false, 'AxisMode', 'image', 'LabelUnits', 'deg', ...
            'PlotRespWin', true, 'RespWin', [-1 4]);
        xlim(axL, limAz); ylim(axL, limEl); set(axL, 'YDir', 'normal');
        set(axL, 'XTick', sort(uAzP), 'YTick', sort(uElP), 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
        flipLabels(axL);
        xlabel(axL, 'Azimuth (\circ)'); ylabel(axL, 'Elevation (\circ)');
        title(axL, sprintf('Bouton %d: %s', res.iROI, plotNames{c}), 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');

        % Right: fitted Gaussian
        pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
        fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));
        axR = axes('Position', [leftMargin+colWidth+hGap+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
        imagesc(axR, xq(1,:), yq(:,1), fittedGrid);
        axis(axR, 'image'); colormap(axR, 'bone'); colorbar(axR);
        xlim(axR, limAz); ylim(axR, limEl); set(axR, 'YDir', 'normal'); hold(axR, 'on');
        plot(axR, res.x0, res.y0, 'r+', 'MarkerSize', 6, 'LineWidth', 2);
        vb = res.bootR2(~isnan(res.bootR2));
        if isempty(vb)
            ttl = sprintf('R^2 = %.2f, \\sigma_{az} = %.0f\\circ, \\sigma_{el} = %.0f\\circ', res.R2, res.sigmaX, res.sigmaY);
        else
            ttl = sprintf('R^2 = %.2f, boot 5th pct = %.2f, \\sigma_{az} = %.0f\\circ, \\sigma_{el} = %.0f\\circ', ...
                res.R2, prctile(vb, 5), res.sigmaX, res.sigmaY);
        end
        title(axR, ttl, 'FontName', 'Arial', 'FontSize', 9);
        set(axR, 'XTick', sort(uAzP), 'YTick', sort(uElP), 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 9);
        flipLabels(axR);
        xlabel(axR, 'Azimuth (\circ)'); ylabel(axR, 'Elevation (\circ)');
    end

    saveFigureFormats(fig, fullfile(outputDir, sprintf('RF_excluded_%s', plotFiles{c})));
    fprintf('Saved %s (%d examples)\n', plotNames{c}, nRows);
end


%% 
RFM  = allRFMapping_rsp;
sess = sessionLabels_rsp;
isResp = [RFM.isResponsive];
isRel  = false(1, numel(RFM));  isRel([gaussFitResults_rsp(trustedMask_rsp).iROI]) = true;

uS = unique(sess);
pctResp = nan(numel(uS), 1);  pctRel = nan(numel(uS), 1);
for s = 1:numel(uS)
    m = strcmp(sess, uS{s});
    pctResp(s) = 100 * mean(isResp(m));
    pctRel(s)  = 100 * mean(isRel(m));
end
fprintf('Responsive:   %.1f +/- %.1f%% (mean +/- SD, n = %d FOVs), range %.1f-%.1f%%\n', ...
    mean(pctResp), std(pctResp), numel(uS), min(pctResp), max(pctResp));
fprintf('Reliable fit: %.1f +/- %.1f%% (mean +/- SD, n = %d FOVs), range %.1f-%.1f%%\n', ...
    mean(pctRel), std(pctRel), numel(uS), min(pctRel), max(pctRel));


%% Exclusion reasons for RF fits (final criteria: R2 >= 0.1)
load('Z:\ibn-vision\USERS\Sonali\Analysis\RFMapping\RF_rsp_boutons.mat');
G   = gaussFitResults_rsp;
RFM = allRFMapping_rsp;

respROIs    = find([RFM.isResponsive]);
nResp       = numel(respROIs);
poorFitROIs = setdiff(respROIs, [G.iROI]);

isRobust = false(1, numel(G));
for gi = 1:numel(G)
    vb = G(gi).bootR2(~isnan(G(gi).bootR2));
    if ~isempty(vb), isRobust(gi) = prctile(vb, 5) >= 0.1; end
end
failBoot   = ~isRobust;
failWidth  = [G.isDegenerateSigma];
failCentre = [G.isBeyondRange];

% Sequential counts (first criterion failed)
seqCounts = [numel(poorFitROIs), ...
             sum(failBoot), ...
             sum(~failBoot & failWidth), ...
             sum(~failBoot & ~failWidth & failCentre), ...
             sum(~failBoot & ~failWidth & ~failCentre)];
seqLabels = {'Poor fit', 'Unreliable', 'Not resolvable', 'Centre outside', 'Reliable'};

% Overlapping counts (of fitted boutons) and "only this criterion"
anyCounts  = [sum(failBoot), sum(failWidth), sum(failCentre)];
onlyCounts = [sum( failBoot & ~failWidth & ~failCentre), ...
              sum(~failBoot &  failWidth & ~failCentre), ...
              sum(~failBoot & ~failWidth &  failCentre)];
ovLabels = {'Unreliable', 'Not resolvable', 'Centre outside'};

fprintf('Sequential (of %d responsive): %s\n', nResp, mat2str(seqCounts));
fprintf('Total: %d\n', sum(seqCounts));

colExcl = [0.6 0.6 0.6];
colRel  = [0.15 0.35 0.35];

figExcl = figure('Color', 'w', 'Position', [100 100 750 300]);

% --- Left: sequential ---
ax1 = subplot(1,2,1);
b = bar(ax1, 100*seqCounts/nResp, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.6);
b.CData = [repmat(colExcl, 4, 1); colRel];
for k = 1:numel(seqCounts)
    text(ax1, k, 100*seqCounts(k)/nResp + 1.5, sprintf('%d', seqCounts(k)), ...
        'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 8);
end
set(ax1, 'XTick', 1:numel(seqLabels), 'XTickLabel', seqLabels, 'XTickLabelRotation', 30, ...
    'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
ylabel(ax1, '% of visually responsive boutons', 'FontName', 'Arial', 'FontSize', 9);
title(ax1, sprintf('First criterion failed (n = %d)', nResp), 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
defaultAxesProperties(ax1, true);

% --- Right: overlapping ---
ax2 = subplot(1,2,2); hold(ax2, 'on');
bar(ax2, 100*anyCounts/numel(G),  'FaceColor', colExcl,     'EdgeColor', 'none', 'BarWidth', 0.6);
bar(ax2, 100*onlyCounts/numel(G), 'FaceColor', colExcl*0.5, 'EdgeColor', 'none', 'BarWidth', 0.6);
for k = 1:3
    text(ax2, k, 100*anyCounts(k)/numel(G) + 1.5, sprintf('%d', anyCounts(k)), ...
        'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 8);
end
set(ax2, 'XTick', 1:3, 'XTickLabel', ovLabels, 'XTickLabelRotation', 30, ...
    'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
ylabel(ax2, '% of fitted boutons', 'FontName', 'Arial', 'FontSize', 9);
title(ax2, sprintf('Any vs only this criterion (n = %d fits)', numel(G)), 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
legend(ax2, {'Failed this criterion', 'Failed only this criterion'}, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
defaultAxesProperties(ax2, true);

% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\excluded_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figExcl, fullfile(outputDir, 'RF_exclusion_reasons'));


%%
saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\data';
files = {'RSP boutons',  'RF_rsp_bouton.mat'; ...
         'VISp boutons', 'RF_v1_boutons.mat'; ...     
         'VISp somata',  'RF_v1_somas.mat'};          
for d = 1:size(files, 1)
    S = load(fullfile(saveDir, files{d, 2}));
    fn = fieldnames(S);

    % Find the variables by name pattern, whatever their suffix
    G    = S.(fn{startsWith(fn, 'gaussFitResults')});
    RFM  = S.(fn{startsWith(fn, 'allRFMapping')});
    sess = S.(fn{startsWith(fn, 'sessionLabels')});

    % Trusted set: same criteria for every dataset (R2 >= 0.1, bootstrap >= 0.1, width, centre)
    isRobust = false(1, numel(G));
    for gi = 1:numel(G)
        vb = G(gi).bootR2(~isnan(G(gi).bootR2));
        if ~isempty(vb), isRobust(gi) = prctile(vb, 5) >= 0.1; end
    end
    tm = ([G.R2] >= 0.1) & isRobust & ~[G.isDegenerateSigma] & ~[G.isBeyondRange];

    nTotal = numel(RFM);
    nResp  = sum([RFM.isResponsive]);
    nRel   = sum(tm);
    sigAz  = [G(tm).sigmaX];   sigEl = [G(tm).sigmaY];
    cAz    = -[G(tm).x0];      cEl   = [G(tm).y0];       % positive azimuth convention
    nFOV   = numel(unique(sess));

    fprintf('\n===== %s (%s, %d field(s) of view) =====\n', files{d,1}, files{d,2}, nFOV);
    fprintf('ROIs recorded:          %d\n', nTotal);
    fprintf('Visually responsive:    %d (%.1f%%)\n', nResp, 100*nResp/nTotal);
    fprintf('Reliable fit:           %d (%.1f%% of all, %.1f%% of responsive)\n', nRel, 100*nRel/nTotal, 100*nRel/nResp);
    fprintf('Width, azimuth:         %.1f +/- %.1f deg (mean +/- SD)\n', mean(sigAz), std(sigAz));
    fprintf('Width, elevation:       %.1f +/- %.1f deg\n', mean(sigEl), std(sigEl));
    fprintf('Centre, azimuth:        %.1f +/- %.1f deg (range %.1f to %.1f)\n', mean(cAz), std(cAz), min(cAz), max(cAz));
    fprintf('Centre, elevation:      %.1f +/- %.1f deg (range %.1f to %.1f)\n', mean(cEl), std(cEl), min(cEl), max(cEl));

    % Within-FOV scatter: distance of each centre from its FOV's mean (FOVs with >= 2 reliable ROIs)
    sessRel = sess([G(tm).iROI]);
    uS = unique(sessRel);
    dAz = nan(1, nRel);  dEl = nan(1, nRel);
    for s = 1:numel(uS)
        m = strcmp(sessRel, uS{s});
        if sum(m) >= 2
            dAz(m) = abs(cAz(m) - mean(cAz(m)));
            dEl(m) = abs(cEl(m) - mean(cEl(m)));
        end
    end
    fprintf('Within-FOV |dAz|:       %.1f +/- %.1f deg (n = %d ROIs)\n', mean(dAz, 'omitnan'), std(dAz, 'omitnan'), sum(~isnan(dAz)));
    fprintf('Within-FOV |dEl|:       %.1f +/- %.1f deg\n', mean(dEl, 'omitnan'), std(dEl, 'omitnan'));
end

%% 
G  = gaussFitResults_rsp;
tm = trustedMask_rsp;
cAz = -[G(tm).x0];  cEl = [G(tm).y0];
sessRel = sessionLabels_rsp([G(tm).iROI]);
uS = unique(sessRel);
dist2D = nan(1, sum(tm));
for s = 1:numel(uS)
    m = strcmp(sessRel, uS{s});
    if sum(m) >= 2
        dist2D(m) = sqrt((cAz(m) - mean(cAz(m))).^2 + (cEl(m) - mean(cEl(m))).^2);
    end
end
fprintf('Within-FOV distance from FOV mean: %.1f +/- %.1f deg (mean +/- SD, n = %d boutons)\n', ...
    mean(dist2D, 'omitnan'), std(dist2D, 'omitnan'), sum(~isnan(dist2D)));


%%
%% Part 1: VR cross-validated R^2 for each RSP bouton, aligned with allRFMapping_rsp
saveDir = 'Z:\ibn-vision\USERS\Sonali\Analysis\RFMapping';
load(fullfile(saveDir, 'RF_rsp_boutons.mat'));

crossValVarName = 'crossValExpVar';
signalName      = 'dFFNeuropilCorrected';

sessOrder = unique(sessionLabels_rsp, 'stable');   % same order as pooling
vrR2_rsp  = nan(numel(allRFMapping_rsp), 1);
vrP_rsp   = nan(numel(allRFMapping_rsp), 1);

for s = 1:numel(sessOrder)
    idxThis = find(strcmp(sessionLabels_rsp, sessOrder{s}));
    parts   = strsplit(sessOrder{s}, '_');           % e.g. 'M25132_20260219'
    mouseID = parts{1};  sessName = parts{2};

    try
        infoPath = findSessionFileInfoFilePath(mouseID, sessName);
        sfi = load(infoPath, 'sessionFileInfo');  sfi = sfi.sessionFileInfo;
        roiPath = sfi.otherSessFilePaths.sessionROIData;
        if ~isfile(roiPath)   % fallback, as in the pooling script
            cand = dir(fullfile(sfi.Directories.save_folder, sprintf('%s_%s_sessionROIData*.mat', sfi.animal_name, sfi.session_name)));
            roiPath = fullfile(cand(1).folder, cand(1).name);
        end
        cv = load(roiPath, crossValVarName);
        sig = cv.(crossValVarName).(signalName);
    catch ME
        warning('%s: could not load VR R^2 (%s)', sessOrder{s}, ME.message);
        continue;
    end

    r2 = sig.medianExpVar(:);  pv = sig.pValues(:);
    if numel(r2) ~= numel(idxThis)
        warning('%s: %d VR values vs %d RF boutons -- skipped', sessOrder{s}, numel(r2), numel(idxThis));
        continue;
    end
    vrR2_rsp(idxThis) = r2;
    vrP_rsp(idxThis)  = pv;
    fprintf('%s: %d boutons matched\n', sessOrder{s}, numel(idxThis));
end

fprintf('\nBoutons with VR R^2: %d / %d\n', sum(~isnan(vrR2_rsp)), numel(vrR2_rsp));
save(fullfile(saveDir, 'RF_rsp_bouton.mat'), 'vrR2_rsp', 'vrP_rsp', '-append');


%% Part 2: VR tuning of visually responsive vs non-responsive RSP boutons
isResp  = [allRFMapping_rsp.isResponsive]';
hasVR   = ~isnan(vrR2_rsp);
vrTuned = vrR2_rsp > 0.1 & vrP_rsp < 0.01;     % Chapter 4 criterion

r2Resp = vrR2_rsp(isResp & hasVR);
r2Non  = vrR2_rsp(~isResp & hasVR);

% 1. VR R^2: responsive vs non-responsive
pR2 = ranksum(r2Resp, r2Non);
fprintf('\nVR R^2, median [IQR]:\n');
fprintf('  Visually responsive (n=%d): %.3f [%.3f, %.3f]\n', numel(r2Resp), median(r2Resp), prctile(r2Resp,25), prctile(r2Resp,75));
fprintf('  Not responsive     (n=%d): %.3f [%.3f, %.3f]\n', numel(r2Non),  median(r2Non),  prctile(r2Non,25),  prctile(r2Non,75));
fprintf('  Wilcoxon rank-sum p = %.3g\n', pR2);

% 2. Fraction reliably tuned in VR
fracResp = mean(vrTuned(isResp & hasVR));
fracNon  = mean(vrTuned(~isResp & hasVR));
[~, chi2, pChi] = crosstab(isResp(hasVR), vrTuned(hasVR));
fprintf('\nReliably tuned in VR:\n');
fprintf('  Visually responsive: %.1f%%\n', 100*fracResp);
fprintf('  Not responsive:      %.1f%%\n', 100*fracNon);
fprintf('  Chi-square = %.2f, p = %.3g\n', chi2, pChi);

% 3. Per field of view (boutons in the same FOV are not independent)
sessOrder = unique(sessionLabels_rsp, 'stable');
perFOV = nan(numel(sessOrder), 2);   % [responsive, not responsive] fraction tuned in VR
for s = 1:numel(sessOrder)
    m = strcmp(sessionLabels_rsp(:), sessOrder{s}) & hasVR;
    if sum(m & isResp) >= 3            % need a few responsive boutons for a fraction
        perFOV(s,1) = mean(vrTuned(m & isResp));
        perFOV(s,2) = mean(vrTuned(m & ~isResp));
    end
end
okFOV = all(~isnan(perFOV), 2);
pFOV  = signrank(perFOV(okFOV,1), perFOV(okFOV,2));
fprintf('\nPer FOV (n=%d FOVs with >= 3 responsive boutons):\n', sum(okFOV));
fprintf('  Tuned in VR, responsive:     %.1f +/- %.1f%%\n', 100*mean(perFOV(okFOV,1)), 100*std(perFOV(okFOV,1)));
fprintf('  Tuned in VR, not responsive: %.1f +/- %.1f%%\n', 100*mean(perFOV(okFOV,2)), 100*std(perFOV(okFOV,2)));
fprintf('  Wilcoxon signed-rank p = %.3g\n', pFOV);

%% Figure
colResp = [0.85 0.33 0.10];  colNon = [0.5 0.5 0.5];
figVR = figure('Color', 'w', 'Position', [100 100 650 280]);

ax1 = subplot(1,2,1); hold(ax1, 'on');
[f1, x1] = ecdf(r2Non);   plot(ax1, x1, f1, 'Color', colNon,  'LineWidth', 1.5);
[f2, x2] = ecdf(r2Resp);  plot(ax1, x2, f2, 'Color', colResp, 'LineWidth', 1.5);
xline(ax1, 0.1, '--k', 'HandleVisibility', 'off');
xlabel(ax1, 'Cross-validated R^2 in VR', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax1, 'Cumulative fraction', 'FontName', 'Arial', 'FontSize', 9);
legend(ax1, {sprintf('Not visually responsive (n=%d)', numel(r2Non)), sprintf('Visually responsive (n=%d)', numel(r2Resp))}, ...
    'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
set(ax1, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
defaultAxesProperties(ax1, true);

ax2 = subplot(1,2,2); hold(ax2, 'on');
plot(ax2, [1 2], 100*perFOV(okFOV,:)', '-', 'Color', [0.8 0.8 0.8]);
scatter(ax2, ones(sum(okFOV),1),   100*perFOV(okFOV,1), 25, colResp, 'filled');
scatter(ax2, 2*ones(sum(okFOV),1), 100*perFOV(okFOV,2), 25, colNon,  'filled');
set(ax2, 'XTick', [1 2], 'XTickLabel', {'Visually responsive', 'Not responsive'}, 'XLim', [0.5 2.5], ...
    'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
ylabel(ax2, '% reliably tuned in VR', 'FontName', 'Arial', 'FontSize', 9);
defaultAxesProperties(ax2, true);

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\vr_comparison';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figVR, fullfile(outputDir, 'VR_R2_responsive_vs_not'));


%%
%% Step 1: running classification of RF mapping trials, per session
saveDir = 'Z:\ibn-vision\USERS\Sonali\Analysis\RFMapping';
load(fullfile(saveDir, 'RF_rsp_boutons.mat'));

runSpeedThresh = 3;  statSpeedThresh = 0.5;  propThresh = 0.75;
classWin = [-0.2 2.8];

sessOrder = unique(sessionLabels_rsp, 'stable');
sessionStats = struct('sessionLabel', {}, 'nGridRun', {}, 'nGridStat', {}, 'nGridAmbig', {}, ...
                      'nBlankRun', {}, 'nBlankStat', {}, 'nBlankAmbig', {});

for s = 1:numel(sessOrder)
    parts = strsplit(sessOrder{s}, '_');
    try
        sfi = load(findSessionFileInfoFilePath(parts{1}, parts{2}), 'sessionFileInfo');
        sfi = sfi.sessionFileInfo;
        rfIdx = find(contains({sfi.stimFiles.name}, 'RFMapping'), 1);
        R = load(sfi.stimFiles(rfIdx).Response, 'response');  R = R.response;
    catch ME
        warning('%s: could not load RF response file (%s)', sessOrder{s}, ME.message);
        continue;
    end
    if ~isfield(R, 'wheelData')
        warning('%s: no wheelData in RF response file', sessOrder{s});
        continue;
    end

    counts = zeros(2, 3);   % rows: grid, blank | cols: run, stat, ambiguous
    for g = 1:numel(R.wheelData)
        sv = R.wheelData(g).stimValue;
        isBlank = sv(1) == 200;                                   % blank = [200 0]
        tVec = R.wheelData(g).timeVector(:);
        win  = tVec >= classWin(1) & tVec <= classWin(2);
        speed = R.wheelData(g).alignedResponses;                  % samples x trials
        for ti = 1:size(speed, 2)
            tr = speed(win, ti);
            if all(isnan(tr)), continue; end
            mS = mean(tr, 'omitnan');
            if mean(tr > statSpeedThresh) >= propThresh && mS > runSpeedThresh
                col = 1;
            elseif mean(tr < runSpeedThresh) >= propThresh && mS < statSpeedThresh
                col = 2;
            else
                col = 3;
            end
            counts(isBlank + 1, col) = counts(isBlank + 1, col) + 1;
        end
    end

    sessionStats(end+1) = struct('sessionLabel', sessOrder{s}, ...
        'nGridRun', counts(1,1), 'nGridStat', counts(1,2), 'nGridAmbig', counts(1,3), ...
        'nBlankRun', counts(2,1), 'nBlankStat', counts(2,2), 'nBlankAmbig', counts(2,3)); %#ok<SAGROW>
    fprintf('%s: run %d | stationary %d | ambiguous %d\n', sessOrder{s}, ...
        counts(1,1)+counts(2,1), counts(1,2)+counts(2,2), counts(1,3)+counts(2,3));
end
save(fullfile(saveDir, 'RF_rsp_bouton.mat'), 'sessionStats', '-append');

%%
% Within M25132 only
isM25132 = startsWith(labelsPreserved, 'M25132');
[rhoM, pM] = corr(xData(isM25132), yData(isM25132), 'type', 'Spearman');
fprintf('Within M25132 (n=%d sessions): Spearman rho = %.2f, p = %.3f\n', sum(isM25132), rhoM, pM);

% Responsiveness in the sessions with the most running
highRun = xData >= 80;
fprintf('Sessions with >= 80%% running (n=%d): responsive %.1f%% to %.1f%%\n', ...
    sum(highRun), min(yData(highRun)), max(yData(highRun)));


%%
%% Receptive fields of boutons within each field of view
G  = gaussFitResults_rsp;
tm = trustedMask_rsp;
cAz = -[G(tm).x0]';  cEl = [G(tm).y0]';
sAz =  [G(tm).sigmaX]';  sEl = [G(tm).sigmaY]';
fov = sessionLabels_rsp([G(tm).iROI]);  fov = fov(:);

uF = unique(fov, 'stable');
nPer = cellfun(@(f) sum(strcmp(fov, f)), uF);
uF = uF(nPer >= 2);                          % FOVs with at least 2 reliable boutons

nCol = 4;  nRow = ceil(numel(uF) / nCol);
figPerFOV = figure('Color', 'w', 'Position', [50 50 1100 250*nRow]);
t = linspace(0, 2*pi, 60);
fprintf('\nWithin-FOV spread (mean absolute deviation from FOV mean):\n');
for f = 1:numel(uF)
    m = strcmp(fov, uF{f});
    ax = subplot(nRow, nCol, f); hold(ax, 'on');
    for k = find(m)'
        patch(ax, cAz(k) + sAz(k)*cos(t), cEl(k) + sEl(k)*sin(t), [0.2 0.35 0.35], ...
            'FaceAlpha', 0.12, 'EdgeColor', [0.2 0.35 0.35], 'EdgeAlpha', 0.5);
    end
    plot(ax, cAz(m), cEl(m), 'k+', 'MarkerSize', 5);
    rectangle(ax, 'Position', [-20 -25 100 65], 'LineStyle', '--', 'EdgeColor', [0.6 0.6 0.6]);  % screen limits
    set(ax, 'XDir', 'reverse', 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 7);
    xlim(ax, [-30 90]);  ylim(ax, [-30 50]);  axis(ax, 'equal');
    title(ax, sprintf('%s (n=%d)', strrep(uF{f}, '_', ' '), sum(m)), 'FontSize', 8, 'FontWeight', 'normal');
    fprintf('  %-18s n=%2d | az %.1f deg | el %.1f deg\n', uF{f}, sum(m), ...
        mean(abs(cAz(m) - mean(cAz(m)))), mean(abs(cEl(m) - mean(cEl(m)))));
end
xlabel(ax, 'Azimuth (\circ)'); ylabel(ax, 'Elevation (\circ)');

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\withinFOV';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figPerFOV, fullfile(outputDir, 'RFs_per_FOV'));

%% Pairwise RF centre distances: same FOV vs different FOV
G  = gaussFitResults_rsp;
tm = trustedMask_rsp;
cAz = -[G(tm).x0]';  cEl = [G(tm).y0]';
fov = sessionLabels_rsp([G(tm).iROI]);  fov = fov(:);
[~, ~, fovID] = unique(fov);

n = numel(cAz);
[I, J] = find(triu(true(n), 1));                         % all unique pairs
dPair  = hypot(cAz(I) - cAz(J), cEl(I) - cEl(J));        % distance between centres (deg)
sameFOV = fovID(I) == fovID(J);

dSame = dPair(sameFOV);  dDiff = dPair(~sameFOV);
fprintf('Pairwise centre distance: same FOV %.1f deg (median, %d pairs) | different FOV %.1f deg (%d pairs)\n', ...
    median(dSame), numel(dSame), median(dDiff), numel(dDiff));

% Permutation test on the median same-FOV distance (shuffle FOV labels)
nPerm = 10000;  rng(1);
nullMed = nan(nPerm, 1);
for k = 1:nPerm
    idS = fovID(randperm(n));
    nullMed(k) = median(dPair(idS(I) == idS(J)));
end
pPerm = (sum(nullMed <= median(dSame)) + 1) / (nPerm + 1);
fprintf('Permutation test: p = %.4f\n', pPerm);

figPair = figure('Color', 'w', 'Position', [100 100 330 300]);
ax = axes; hold(ax, 'on');
[f1, x1] = ecdf(dDiff);  plot(ax, x1, f1, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
[f2, x2] = ecdf(dSame);  plot(ax, x2, f2, 'Color', [0.15 0.35 0.35], 'LineWidth', 1.5);
xlabel(ax, 'Distance between RF centres (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax, 'Cumulative fraction of pairs', 'FontName', 'Arial', 'FontSize', 9);
legend(ax, {sprintf('Different FOV (%d pairs)', numel(dDiff)), sprintf('Same FOV (%d pairs)', numel(dSame))}, ...
    'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
title(ax, sprintf('Permutation p = %.3g', pPerm), 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8);
defaultAxesProperties(ax, true);

outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\withinFOV';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figPair, fullfile(outputDir, 'pairwiseRFdistance_sameVsDiffFOV'));