%% Excluded examples (RSP boutons), same style as the example figure
G   = gaussFitResults_rsp;
RFM = allRFMapping_rsp;
uAzP = uAz_rsp;  uElP = uEl_rsp;

% Definitions needed for plotting (in case they're not in the workspace)
gaussFit2D = @(p, xy) p(1) * exp( -( (xy(:,1)-p(2)).^2 ./ (2*p(4)^2) + (xy(:,2)-p(3)).^2 ./ (2*p(5)^2) ) );
[AzGrid, ~] = meshgrid(uAzP, uElP);
if ~exist('timeVector', 'var'), error('timeVector not in workspace: run the pooling script, or save timeVector with the data.'); end
respIdx = timeVector >= 0.1 & timeVector <= 3;

% Exclusion reason for each fit (final criteria: R2 >= 0.1, bootstrap >= 0.1, resolvable, centre within stimulated area)
isRobust = false(1, numel(G));
for gi = 1:numel(G)
    vb = G(gi).bootR2(~isnan(G(gi).bootR2));
    if ~isempty(vb), isRobust(gi) = prctile(vb, 5) >= 0.1; end
end
isDegen  = [G.isDegenerateSigma];
isBeyond = [G.isBeyondRange];
R2all    = [G.R2];

reasonIdx = { find(~isRobust), ...
              find(isRobust & isDegen), ...
              find(isRobust & ~isDegen & isBeyond) };
reasonNames = {'Unreliable across trials', 'Not resolvable', 'Centre outside stimulated area'};
for c = 1:3
    [~, o] = sort(R2all(reasonIdx{c}), 'descend');   % clearest examples first
    reasonIdx{c} = reasonIdx{c}(o);
    fprintf('%s: %d boutons\n', reasonNames{c}, numel(reasonIdx{c}));
end

% Take examples from each reason in turn, up to 10
nWanted = 10;
pickGI = []; pickReason = {};
k = 1;
while numel(pickGI) < nWanted && any(cellfun(@numel, reasonIdx) >= k)
    for c = 1:3
        if numel(reasonIdx{c}) >= k && numel(pickGI) < nWanted
            pickGI(end+1) = reasonIdx{c}(k); %#ok<SAGROW>
            pickReason{end+1} = reasonNames{c}; %#ok<SAGROW>
        end
    end
    k = k + 1;
end

% Shared axis limits: tested grid +/- half a step
azStepP = mean(diff(sort(uAzP)));  elStepP = mean(diff(sort(uElP)));
limAz = [min(uAzP) - azStepP/2, max(uAzP) + azStepP/2];
limEl = [min(uElP) - elStepP/2, max(uElP) + elStepP/2];
[xq, yq] = meshgrid(linspace(limAz(1), limAz(2), 100), linspace(limEl(1), limEl(2), 100));

rowsPerFig = 5;
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\excluded_rfs';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

for f = 1:ceil(numel(pickGI) / rowsPerFig)
    theseRows = (f-1)*rowsPerFig + 1 : min(f*rowsPerFig, numel(pickGI));
    nRows = numel(theseRows);
    fig = figure('Color', 'w', 'Position', [100 50 900 200*nRows + 100]);

    leftMargin = 0.08; rightMargin = 0.05; topMargin = 0.06; bottomMargin = 0.06;
    hGap = 0.08; vGap = 0.05;
    colWidth  = (1 - leftMargin - rightMargin - hGap) / 2;
    rowHeight = (1 - topMargin - bottomMargin - (nRows-1)*vGap) / nRows;
    insetFrac = 0.85;
    padW = colWidth * (1 - insetFrac) / 2;  padH = rowHeight * (1 - insetFrac) / 2;

    for r = 1:nRows
        gi   = pickGI(theseRows(r));
        res  = G(gi);
        bData = RFM(res.iROI);
        yPos = 1 - topMargin - r*rowHeight - (r-1)*vGap;

        % Mean response at each grid position
        tm = bData.baselineSubtracted;
        respVec = nan(numel(tm), 1);
        for p = 1:numel(tm)
            if ~isempty(tm{p})
                respVec(p) = mean(mean(double(tm{p}(:, respIdx)), 2, 'omitnan'), 'omitnan');
            end
        end
        respGrid = reshape(respVec, size(AzGrid));

        boutonData = bData;
        boutonData.meanGridResponse = respGrid;
        boutonData.peakAmplitude    = max(respGrid(:)) * 2;

        % Left: grid with traces
        axL = axes('Position', [leftMargin+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
        plotRFHeatmapWithTraces(axL, boutonData, uAzP, uElP, timeVector, ...
            'Colormap', 'bone', 'Smooth', false, 'AxisMode', 'image', 'LabelUnits', 'deg', ...
            'PlotRespWin', true, 'RespWin', [-1 4]);
        xlim(axL, limAz); ylim(axL, limEl); set(axL, 'YDir', 'normal');
        set(axL, 'XTick', sort(uAzP), 'YTick', sort(uElP), 'Box', 'off', 'TickDir', 'out', ...
            'FontName', 'Arial', 'FontSize', 9);
        title(axL, sprintf('Bouton %d: %s', res.iROI, pickReason{theseRows(r)}), ...
            'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');

        % Right: fitted Gaussian
        pFit = [res.A, res.x0, res.y0, res.sigmaX, res.sigmaY];
        fittedGrid = reshape(gaussFit2D(pFit, [xq(:), yq(:)]), size(xq));
        axR = axes('Position', [leftMargin+colWidth+hGap+padW, yPos+padH, colWidth*insetFrac, rowHeight*insetFrac]);
        imagesc(axR, xq(1,:), yq(:,1), fittedGrid);
        axis(axR, 'image'); colormap(axR, 'bone'); colorbar(axR);
        xlim(axR, limAz); ylim(axR, limEl); set(axR, 'YDir', 'normal'); hold(axR, 'on');
        plot(axR, res.x0, res.y0, 'r+', 'MarkerSize', 6, 'LineWidth', 2);
        vb = res.bootR2(~isnan(res.bootR2));
        title(axR, sprintf('R^2 = %.2f, boot 5th pct = %.2f, \\sigma_{az} = %.0f\\circ, \\sigma_{el} = %.0f\\circ', ...
            res.R2, prctile(vb, 5), res.sigmaX, res.sigmaY), 'FontName', 'Arial', 'FontSize', 9);
        xlabel(axR, 'Azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
        ylabel(axR, 'Elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
        set(axR, 'XTick', sort(uAzP), 'YTick', sort(uElP), 'Box', 'off', 'TickDir', 'out', ...
            'FontName', 'Arial', 'FontSize', 9);
    end

    saveFigureFormats(fig, fullfile(outputDir, sprintf('RF_excluded_examples_%d', f)));
end




%% Peak azimuth and elevation (grid position of the largest response), all responsive boutons
RFM  = allRFMapping_rsp;
G    = gaussFitResults_rsp;
uAzP = sort(uAz_rsp(:))';
uElP = sort(uEl_rsp(:))';

isResp = [RFM.isResponsive];
isTrustedROI = false(1, numel(RFM));
isTrustedROI([G(trustedMask_rsp).iROI]) = true;

peakAz_resp  = [RFM(isResp).centerAz];
peakEl_resp  = [RFM(isResp).centerEl];
peakAz_trust = [RFM(isResp & isTrustedROI).centerAz];
peakEl_trust = [RFM(isResp & isTrustedROI).centerEl];

% One bin per grid position
azStepP = mean(diff(uAzP));  elStepP = mean(diff(uElP));
azEdgesP = [uAzP - azStepP/2, uAzP(end) + azStepP/2];
elEdgesP = [uElP - elStepP/2, uElP(end) + elStepP/2];

colResp  = [0.6 0.6 0.6];
colTrust = [0.85 0.33 0.10];

figPeak = figure('Position', [100, 100, 600, 280]);

% --- Azimuth ---
ax1 = subplot(1,2,1); hold(ax1, 'on');
histogram(ax1, 'BinEdges', azEdgesP, 'BinCounts', 100*histcounts(peakAz_resp, azEdgesP)/numel(peakAz_resp), ...
    'FaceColor', colResp, 'FaceAlpha', 0.6, 'EdgeColor', colResp*0.6, 'LineWidth', 1, ...
    'DisplayName', sprintf('Responsive (n=%d)', numel(peakAz_resp)));
histogram(ax1, 'BinEdges', azEdgesP, 'BinCounts', 100*histcounts(peakAz_trust, azEdgesP)/numel(peakAz_trust), ...
    'FaceColor', colTrust, 'FaceAlpha', 0.6, 'EdgeColor', colTrust*0.6, 'LineWidth', 1, ...
    'DisplayName', sprintf('Trusted (n=%d)', numel(peakAz_trust)));
xlabel(ax1, 'Peak azimuth (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax1, '% of boutons', 'FontName', 'Arial', 'FontSize', 9);
legend(ax1, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
set(ax1, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', uAzP);
defaultAxesProperties(ax1, true);

% --- Elevation ---
ax2 = subplot(1,2,2); hold(ax2, 'on');
histogram(ax2, 'BinEdges', elEdgesP, 'BinCounts', 100*histcounts(peakEl_resp, elEdgesP)/numel(peakEl_resp), ...
    'FaceColor', colResp, 'FaceAlpha', 0.6, 'EdgeColor', colResp*0.6, 'LineWidth', 1, ...
    'DisplayName', sprintf('Responsive (n=%d)', numel(peakEl_resp)));
histogram(ax2, 'BinEdges', elEdgesP, 'BinCounts', 100*histcounts(peakEl_trust, elEdgesP)/numel(peakEl_trust), ...
    'FaceColor', colTrust, 'FaceAlpha', 0.6, 'EdgeColor', colTrust*0.6, 'LineWidth', 1, ...
    'DisplayName', sprintf('Trusted (n=%d)', numel(peakEl_trust)));
xlabel(ax2, 'Peak elevation (\circ)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax2, '% of boutons', 'FontName', 'Arial', 'FontSize', 9);
set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 8, 'XTick', uElP);
defaultAxesProperties(ax2, true);

% --- Save ---
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter4-RSP-VisualStim\Section1_Fig4_1\pref_az_el_peaksFromData';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
saveFigureFormats(figPeak, fullfile(outputDir, 'peak_az_el_responsive_vs_trusted'));