function plotPeaksCorridorSplit(RSPData, VISpData, varargin)
    % Added an SMI threshold so this analysis can be restricted to
    % modulated cells only (e.g. to test whether cells with a strong
    % position preference show a first-vs-second-half bias, as in
    % Diamanti et al.'s contrast adaptation analysis), instead of pooling
    % every cell with a valid peak regardless of modulation strength.
    p = inputParser;
    addRequired(p, 'RSPData', @isstruct);
    addRequired(p, 'VISpData', @isstruct);
    addParameter(p, 'SMIThreshold', [], @(x) isempty(x) || isnumeric(x)); % e.g. 0.25; [] = no filtering (all cells)
    parse(p, RSPData, VISpData, varargin{:});
    smiThreshold = p.Results.SMIThreshold;

    trackStart = 30;
    trackEnd = 170;
    midPoint = 100; % Exact middle of the 30-170cm valid corridor

    fprintf('Extracting RSP peaks from SMI_Metrics fields...\n');
    rspPeaks = extractSavedStructPeaks(RSPData, trackStart, trackEnd, smiThreshold);

    fprintf('Extracting VISp peaks from SMI_Metrics fields...\n');
    vispPeaks = extractSavedStructPeaks(VISpData, trackStart, trackEnd, smiThreshold);

    %% Per-session extraction (shared by scatter plot, printed breakdown, and stats)
    rspSession  = computeSessionHalfStats(RSPData, trackStart, midPoint, trackEnd, smiThreshold);
    vispSession = computeSessionHalfStats(VISpData, trackStart, midPoint, trackEnd, smiThreshold);

    if isempty(rspSession.pctFirst) || isempty(vispSession.pctFirst)
        error('No sessions had cells passing the SMI threshold (%.3g) — try relaxing SMIThreshold.', smiThreshold);
    end

    rspMeanFirst   = mean(rspSession.pctFirst);
    rspMeanSecond  = mean(rspSession.pctSecond);
    vispMeanFirst  = mean(vispSession.pctFirst);
    vispMeanSecond = mean(vispSession.pctSecond);

    p_RSP_1st_vs_2nd  = signrank(rspSession.pctFirst,  rspSession.pctSecond);
    p_VISp_1st_vs_2nd = signrank(vispSession.pctFirst, vispSession.pctSecond);
    p_First_RSPvsVISp  = ranksum(rspSession.pctFirst,  vispSession.pctFirst);
    p_Second_RSPvsVISp = ranksum(rspSession.pctSecond, vispSession.pctSecond);

    if isempty(smiThreshold)
        thrLabel = 'all cells';
    else
        thrLabel = sprintf('SMI > %.3g', smiThreshold);
    end

    fprintf('\n=== Peak Position Half-Split: Session-Level Statistics (%s) ===\n', thrLabel);
    fprintf('RSP:  First %.1f +/- %.1f%% (SD)  Second %.1f +/- %.1f%% (SD)  [n=%d sessions]  signrank p = %.3e\n', ...
        rspMeanFirst, std(rspSession.pctFirst), rspMeanSecond, std(rspSession.pctSecond), ...
        length(rspSession.pctFirst), p_RSP_1st_vs_2nd);
    fprintf('VISp: First %.1f +/- %.1f%% (SD)  Second %.1f +/- %.1f%% (SD)  [n=%d sessions]  signrank p = %.3e\n', ...
        vispMeanFirst, std(vispSession.pctFirst), vispMeanSecond, std(vispSession.pctSecond), ...
        length(vispSession.pctFirst), p_VISp_1st_vs_2nd);
    fprintf('Cross-region: First-half RSP vs VISp ranksum p = %.3e | Second-half RSP vs VISp ranksum p = %.3e\n', ...
        p_First_RSPvsVISp, p_Second_RSPvsVISp);

    % CHANGED: pooled cell-level first-half bias among modulated cells,
    % matching Diamanti et al.'s "X/Y fired more strongly in the first
    % half" reporting style, with a binomial test against 50/50.
    if ~isempty(smiThreshold)
        rspFirstCount  = sum(rspPeaks <= midPoint);
        rspTotalCount  = length(rspPeaks);
        vispFirstCount = sum(vispPeaks <= midPoint);
        vispTotalCount = length(vispPeaks);

        pRSPBinom  = 2 * min(binocdf(rspFirstCount, rspTotalCount, 0.5), ...
                              1 - binocdf(rspFirstCount - 1, rspTotalCount, 0.5));
        pVISpBinom = 2 * min(binocdf(vispFirstCount, vispTotalCount, 0.5), ...
                              1 - binocdf(vispFirstCount - 1, vispTotalCount, 0.5));

        fprintf('\n=== First-Half Bias Among Modulated Cells (%s) ===\n', thrLabel);
        fprintf('RSP boutons:  %d/%d (%.1f%%) peaked in the first half [binomial test vs 50/50: p = %.3g]\n', ...
            rspFirstCount, rspTotalCount, (rspFirstCount/rspTotalCount)*100, pRSPBinom);
        fprintf('VISp somas:   %d/%d (%.1f%%) peaked in the first half [binomial test vs 50/50: p = %.3g]\n', ...
            vispFirstCount, vispTotalCount, (vispFirstCount/vispTotalCount)*100, pVISpBinom);
    end

    barData = [rspMeanFirst, vispMeanFirst; rspMeanSecond, vispMeanSecond];

    %% plotting
    fig = figure('Color', 'w', 'Position', [150 150 550 480]);
    hold on;

    b = bar(barData, 'EdgeColor', 'none', 'BarWidth', 0.8);
    b(1).FaceColor = 'k';           % RSP = Black
    b(2).FaceColor = [0.6 0.6 0.6]; % VISp = Gray

    jitterWidth = 0.12;
    scatterFaceAlpha = 0.6;

    rspX_first   = b(1).XEndPoints(1) + (rand(size(rspSession.pctFirst))  - 0.5) * jitterWidth;
    rspX_second  = b(1).XEndPoints(2) + (rand(size(rspSession.pctSecond)) - 0.5) * jitterWidth;
    vispX_first  = b(2).XEndPoints(1) + (rand(size(vispSession.pctFirst))  - 0.5) * jitterWidth;
    vispX_second = b(2).XEndPoints(2) + (rand(size(vispSession.pctSecond)) - 0.5) * jitterWidth;

    scatter(rspX_first,   rspSession.pctFirst,   20, 'k', 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    scatter(rspX_second,  rspSession.pctSecond,  20, 'k', 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    scatter(vispX_first,  vispSession.pctFirst,  20, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    scatter(vispX_second, vispSession.pctSecond, 20, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', scatterFaceAlpha, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');

    rspSDFirst   = std(rspSession.pctFirst);
    rspSDSecond  = std(rspSession.pctSecond);
    vispSDFirst  = std(vispSession.pctFirst);
    vispSDSecond = std(vispSession.pctSecond);

    errorbar([b(1).XEndPoints(1), b(1).XEndPoints(2)], ...
             [rspMeanFirst, rspMeanSecond], ...
             [rspSDFirst, rspSDSecond], ...
             'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.2, 'CapSize', 8, 'HandleVisibility', 'off');
    errorbar([b(2).XEndPoints(1), b(2).XEndPoints(2)], ...
             [vispMeanFirst, vispMeanSecond], ...
             [vispSDFirst, vispSDSecond], ...
             'LineStyle', 'none', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'CapSize', 8, 'HandleVisibility', 'off');

    topOfBars = max([rspMeanFirst+rspSDFirst, rspMeanSecond+rspSDSecond, ...
                      vispMeanFirst+vispSDFirst, vispMeanSecond+vispSDSecond]);
    stepHeight = 4;

    y1 = topOfBars + stepHeight;
    y2 = topOfBars + stepHeight;
    y3 = topOfBars + stepHeight * 2.2;
    y4 = topOfBars + stepHeight * 3.4;

    drawSigBracket(b(1).XEndPoints(1), b(1).XEndPoints(2), y1, p_RSP_1st_vs_2nd, 'k');
    drawSigBracket(b(2).XEndPoints(1), b(2).XEndPoints(2), y2, p_VISp_1st_vs_2nd, [0.4 0.4 0.4]);
    drawSigBracket(b(1).XEndPoints(1), b(2).XEndPoints(1), y3, p_First_RSPvsVISp, [0.2 0.2 0.2]);
    drawSigBracket(b(1).XEndPoints(2), b(2).XEndPoints(2), y4, p_Second_RSPvsVISp, [0.2 0.2 0.2]);

    set(gca, 'Box', 'off', 'XTick', 1:2, ...
             'XTickLabel', {sprintf('First Half (%d-%d cm)', trackStart, midPoint), ...
                            sprintf('Second Half (%d-%d cm)', midPoint, trackEnd)}, ...
             'FontName', 'Arial', 'FontSize', 11);
    ylabel('% of ROIs (mean +/- SD, dots = sessions)', 'FontName', 'Arial', 'FontSize', 12);
    ylim([0, y4 + stepHeight*1.5]);

    title(sprintf('SMI RespP bin (%s)', thrLabel), 'FontWeight', 'normal', 'FontSize', 11, 'FontName', 'Arial');

    legend(b, {sprintf('RSP (n=%d sessions)', length(rspSession.pctFirst)), ...
               sprintf('VISp (n=%d sessions)', length(vispSession.pctFirst))}, ...
        'Location', 'southoutside', 'Box', 'off', 'FontName', 'Arial');

    axis square; box off;
    if exist('defaultAxesProperties', 'file') == 2, defaultAxesProperties(gca); end
    if exist('offsetAxes', 'file') == 2, offsetAxes(gca); end
    hold off;

    if isempty(smiThreshold)
        fileSuffix = 'allcells';
    else
        fileSuffix = sprintf('SMIgt%s', strrep(num2str(smiThreshold), '.', 'p'));
    end
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.0\peak_lap-halves_rsp_visp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, sprintf('peakbin_halves_within_lap_rsp_visp_%s', fileSuffix)));

    %% Per-session breakdown — check if bias is driven by one session/animal
    fprintf('\n=== Per-Session Peak Half Analysis (%s) ===\n', thrLabel);
    printSessionHalfStats('RSP', rspSession);
    printSessionHalfStats('VISp', vispSession);

end

%% Draw a significance bracket with stars (or "n.s.") between two x positions
function drawSigBracket(x1, x2, y, pVal, lineColor)
    tickHeight = 0.6;
    plot([x1, x1, x2, x2], [y, y+tickHeight, y+tickHeight, y], ...
        'Color', lineColor, 'LineWidth', 1, 'HandleVisibility', 'off');
    text(mean([x1, x2]), y + tickHeight*1.3, pToStars(pVal), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontName', 'Arial', 'Color', lineColor);
end

%% Convert a p-value to significance stars, Bonferroni-corrected for 4 comparisons
function starStr = pToStars(pVal)
    numComparisons = 4;
    alpha001 = 0.001 / numComparisons;
    alpha01  = 0.01  / numComparisons;
    alpha05  = 0.05  / numComparisons;

    if pVal <= alpha001
        starStr = '***';
    elseif pVal <= alpha01
        starStr = '**';
    elseif pVal <= alpha05
        starStr = '*';
    else
        starStr = 'n.s.';
    end
end

%% Use the RespP bin already computed for SMI analyses.
%% Added smiThreshold — a cell is only included if smiThreshold is empty
%% (no filtering) or its SMI exceeds smiThreshold.
function [validPeaks] = extractSavedStructPeaks(RegionData, trackStart, trackEnd, smiThreshold)
    if nargin < 4, smiThreshold = []; end
    validPeaks = [];
    for s = 1:length(RegionData)
        sess = RegionData(s);
        if ~isfield(sess, 'SMI') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end
        RespP = sess.SMI.RpBin;
        smiVals = sess.SMI.SMI;
        excludeFlags = sess.SMI.ExcludeEdgePeakCells;
        roisToAnalyze = sess.FilteredROIs;
        for i = 1:length(roisToAnalyze)
            roiIdx = roisToAnalyze(i);
            if excludeFlags(roiIdx)
                continue;
            end
            if ~isempty(smiThreshold) && (isnan(smiVals(roiIdx)) || smiVals(roiIdx) <= smiThreshold)
                continue;
            end
            cellPeakBin = RespP(roiIdx);
            if cellPeakBin >= trackStart && cellPeakBin <= trackEnd
                validPeaks = [validPeaks; cellPeakBin];
            end
        end
    end
end

%% Per-session first-half/second-half stats, shared by plotting and the printed table.
%% Added smiThreshold, same filtering logic as extractSavedStructPeaks.
function stats = computeSessionHalfStats(Data, trackStart, midPoint, trackEnd, smiThreshold)
    if nargin < 5, smiThreshold = []; end
    pctFirst  = [];
    pctSecond = [];
    nCellsArr = [];
    labels    = {};

    for s = 1:length(Data)
        sess = Data(s);
        if ~isfield(sess, 'SMI') || ~isfield(sess, 'FilteredROIs') || isempty(sess.FilteredROIs)
            continue;
        end

        RespP        = sess.SMI.RpBin;
        smiVals      = sess.SMI.SMI;
        excludeFlags = sess.SMI.ExcludeEdgePeakCells;
        roisToAnalyze = sess.FilteredROIs;

        sessionPeaks = [];
        for i = 1:length(roisToAnalyze)
            roiIdx = roisToAnalyze(i);
            if excludeFlags(roiIdx), continue; end
            if ~isempty(smiThreshold) && (isnan(smiVals(roiIdx)) || smiVals(roiIdx) <= smiThreshold)
                continue;
            end
            cellPeakBin = RespP(roiIdx);
            if cellPeakBin >= trackStart && cellPeakBin <= trackEnd
                sessionPeaks = [sessionPeaks; cellPeakBin];
            end
        end

        if isempty(sessionPeaks), continue; end

        nCells  = length(sessionPeaks);
        nFirst  = sum(sessionPeaks <= midPoint);
        nSecond = sum(sessionPeaks > midPoint);

        pctFirst(end+1)  = (nFirst / nCells) * 100;  %#ok<AGROW>
        pctSecond(end+1) = (nSecond / nCells) * 100; %#ok<AGROW>
        nCellsArr(end+1) = nCells;                   %#ok<AGROW>

        if isfield(sess, 'sessionName')
            labels{end+1} = sess.sessionName; %#ok<AGROW>
        elseif isfield(sess, 'animal')
            labels{end+1} = sprintf('%s_s%d', sess.animal, s); %#ok<AGROW>
        else
            labels{end+1} = sprintf('Session_%d', s); %#ok<AGROW>
        end
    end

    stats.pctFirst  = pctFirst;
    stats.pctSecond = pctSecond;
    stats.nCells    = nCellsArr;
    stats.labels    = labels;
end

%% Print the per-session breakdown table for one region
function printSessionHalfStats(regionName, stats)
    fprintf('\n--- %s ---\n', regionName);
    fprintf('%-20s %-12s %-12s %-12s %-12s\n', 'Session', 'N Cells', 'First Half', 'Second Half', '% First');
    for s = 1:length(stats.labels)
        nFirst  = round(stats.pctFirst(s)  / 100 * stats.nCells(s));
        nSecond = round(stats.pctSecond(s) / 100 * stats.nCells(s));
        fprintf('%-20s %-12d %-12d %-12d %-12.1f\n', ...
            stats.labels{s}, stats.nCells(s), nFirst, nSecond, stats.pctFirst(s));
    end
end