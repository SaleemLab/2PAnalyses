function [incIdx_o2, decIdx_o2, deltas_o2, incTable_o2, decTable_o2, ...
          incIdx_o3, decIdx_o3, deltas_o3, incTable_o3, decTable_o3, ...
          figHandle] = findLandmarkModulatedROIs(RSPData, baseName, omitName2, omitName3, varargin)

% Gemini: 
% For each bouton, compute the mean activity in a window around the landmark position (±15 cm by default) — 
% separately for Base and Omit conditions, averaged across all laps
% Delta = mean(Omit window) - mean(Base window) per bouton
% Threshold = 1 std of the delta distribution across all boutons
%Increase = delta > +threshold (more active at landmark position when it's omitted)
%Decrease = delta < -threshold (less active at landmark position when it's omitted)

    p = inputParser;
    addParameter(p, 'windowCm',    15,          @isnumeric);
    addParameter(p, 'threshold',   1,           @isnumeric);
    addParameter(p, 'trackLength', 200,         @isnumeric);
    addParameter(p, 'climAbs',     [0 1], @isnumeric);
    addParameter(p, 'SavePath',    '',          @(x) ischar(x) || isstring(x));
    addParameter(p, 'FontName',    'Arial',     @(x) ischar(x) || isstring(x));
    addParameter(p, 'smoothWin',   11,          @isnumeric);
    parse(p, varargin{:});

    windowCm    = p.Results.windowCm;
    threshold   = p.Results.threshold;
    trackLength = p.Results.trackLength;
    climAbs     = p.Results.climAbs;
    targetFont  = p.Results.FontName;
    smoothWin   = p.Results.smoothWin;

    landmarkPos = [80, 120];
    omitNames   = {omitName2, omitName3};
    omitColors  = {[0.541, 0.012, 0.012], [0.824, 0.016, 0.176]};
    omitLabels  = {strrep(omitName2,'_',' '), strrep(omitName3,'_',' ')};

    %% Collect data across sessions
    allBaseAct    = cell(1,2);
    allBaseOdd    = cell(1,2);
    allOmitAct    = cell(1,2);
    allDeltas     = cell(1,2);
    allSessionIDs = cell(1,2);
    allROIIDs     = cell(1,2);

    for oi = 1:2
        omitName = omitNames{oi};
        lPos     = landmarkPos(oi);

        baseAct_all    = [];
        baseOdd_all    = [];
        omitAct_all    = [];
        deltas_all     = [];
        sessionIDs_all = [];
        roiIDs_all     = [];

        for s = 1:length(RSPData)
            sess = RSPData(s);

            if ~isfield(sess.ConditionData, baseName) || ~isfield(sess.ConditionData, omitName)
                fprintf('  Session %d: missing %s — skipping\n', s, omitName);
                continue;
            end

            if isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
                idx = sess.FilteredROIs;
            else
                idx = 1:size(sess.ConditionData.(baseName).LapActivity, 1);
            end

            nBins   = size(sess.ConditionData.(baseName).LapActivity, 3);
            posBins = linspace(0, trackLength, nBins);
            winMask = posBins >= (lPos - windowCm) & posBins <= (lPos + windowCm);

            % Smooth
            baseLapSm = smoothLapActivity(sess.ConditionData.(baseName).LapActivity(idx, :, :), smoothWin);
            omitLapSm = smoothLapActivity(sess.ConditionData.(omitName).LapActivity(idx, :, :), smoothWin);

            nLaps = size(baseLapSm, 2);

            % Delta
            baseLaps = squeeze(mean(baseLapSm(:, :, winMask), 3, 'omitnan'));
            omitLaps = squeeze(mean(omitLapSm(:, :, winMask), 3, 'omitnan'));
            if size(baseLaps,1) ~= numel(idx), baseLaps = baseLaps'; omitLaps = omitLaps'; end
            sessDeltas = mean(omitLaps, 2, 'omitnan') - mean(baseLaps, 2, 'omitnan');

            % Mean across all laps
            baseAct = squeeze(mean(baseLapSm, 2, 'omitnan'));
            omitAct = squeeze(mean(omitLapSm, 2, 'omitnan'));
            if size(baseAct,1) ~= numel(idx), baseAct = baseAct'; omitAct = omitAct'; end

            % Odd laps only
            baseOdd = squeeze(mean(baseLapSm(:, 1:2:nLaps, :), 2, 'omitnan'));
            if size(baseOdd,1) ~= numel(idx), baseOdd = baseOdd'; end

            baseAct_all    = [baseAct_all;    baseAct];                   %#ok<AGROW>
            baseOdd_all    = [baseOdd_all;    baseOdd];                   %#ok<AGROW>
            omitAct_all    = [omitAct_all;    omitAct];                   %#ok<AGROW>
            deltas_all     = [deltas_all;     sessDeltas];                %#ok<AGROW>
            sessionIDs_all = [sessionIDs_all; repmat(s, numel(idx), 1)]; %#ok<AGROW>
            roiIDs_all     = [roiIDs_all;     idx(:)];                   %#ok<AGROW>
        end

        allBaseAct{oi}    = baseAct_all;
        allBaseOdd{oi}    = baseOdd_all;
        allOmitAct{oi}    = omitAct_all;
        allDeltas{oi}     = deltas_all;
        allSessionIDs{oi} = sessionIDs_all;
        allROIIDs{oi}     = roiIDs_all;
    end

    %% Threshold
    thresh_o2   = threshold * std(allDeltas{1}, 'omitnan');
    incIdx_o2   = find(allDeltas{1} >  thresh_o2);
    decIdx_o2   = find(allDeltas{1} < -thresh_o2);
    deltas_o2   = allDeltas{1};

    thresh_o3   = threshold * std(allDeltas{2}, 'omitnan');
    incIdx_o3   = find(allDeltas{2} >  thresh_o3);
    decIdx_o3   = find(allDeltas{2} < -thresh_o3);
    deltas_o3   = allDeltas{2};

    incTable_o2 = table(allSessionIDs{1}(incIdx_o2), allROIIDs{1}(incIdx_o2), allDeltas{1}(incIdx_o2), ...
        'VariableNames', {'SessionIdx','ROI_idx','delta'});
    decTable_o2 = table(allSessionIDs{1}(decIdx_o2), allROIIDs{1}(decIdx_o2), allDeltas{1}(decIdx_o2), ...
        'VariableNames', {'SessionIdx','ROI_idx','delta'});
    incTable_o3 = table(allSessionIDs{2}(incIdx_o3), allROIIDs{2}(incIdx_o3), allDeltas{2}(incIdx_o3), ...
        'VariableNames', {'SessionIdx','ROI_idx','delta'});
    decTable_o3 = table(allSessionIDs{2}(decIdx_o3), allROIIDs{2}(decIdx_o3), allDeltas{2}(decIdx_o3), ...
        'VariableNames', {'SessionIdx','ROI_idx','delta'});

    fprintf('\n--- Omit 2 (80cm) --- Inc: %d | Dec: %d\n', numel(incIdx_o2), numel(decIdx_o2));
    fprintf('\n--- Omit 3 (120cm) --- Inc: %d | Dec: %d\n', numel(incIdx_o3), numel(decIdx_o3));

    %% Normalise
    normMat = @(M, refM) (M - min(refM,[],2)) ./ max(max(refM,[],2) - min(refM,[],2), eps);
    % Z-score per bouton using baseline mean and std
    %normMat = @(M, refM) (M - mean(refM, 2, 'omitnan')) ./ std(refM, 0, 2, 'omitnan');



    baseInc_o2    = normMat(allBaseAct{1}(incIdx_o2,:), allBaseAct{1}(incIdx_o2,:));
    omitInc_o2    = normMat(allOmitAct{1}(incIdx_o2,:), allBaseAct{1}(incIdx_o2,:));
    diffInc_o2    = omitInc_o2 - baseInc_o2;
    baseIncOdd_o2 = normMat(allBaseOdd{1}(incIdx_o2,:), allBaseOdd{1}(incIdx_o2,:));

    baseDec_o2    = normMat(allBaseAct{1}(decIdx_o2,:), allBaseAct{1}(decIdx_o2,:));
    omitDec_o2    = normMat(allOmitAct{1}(decIdx_o2,:), allBaseAct{1}(decIdx_o2,:));
    diffDec_o2    = omitDec_o2 - baseDec_o2;
    baseDecOdd_o2 = normMat(allBaseOdd{1}(decIdx_o2,:), allBaseOdd{1}(decIdx_o2,:));

    baseInc_o3    = normMat(allBaseAct{2}(incIdx_o3,:), allBaseAct{2}(incIdx_o3,:));
    omitInc_o3    = normMat(allOmitAct{2}(incIdx_o3,:), allBaseAct{2}(incIdx_o3,:));
    diffInc_o3    = omitInc_o3 - baseInc_o3;
    baseIncOdd_o3 = normMat(allBaseOdd{2}(incIdx_o3,:), allBaseOdd{2}(incIdx_o3,:));

    baseDec_o3    = normMat(allBaseAct{2}(decIdx_o3,:), allBaseAct{2}(decIdx_o3,:));
    omitDec_o3    = normMat(allOmitAct{2}(decIdx_o3,:), allBaseAct{2}(decIdx_o3,:));
    diffDec_o3    = omitDec_o3 - baseDec_o3;
    baseDecOdd_o3 = normMat(allBaseOdd{2}(decIdx_o3,:), allBaseOdd{2}(decIdx_o3,:));

    posBins = linspace(0, trackLength, size(allBaseAct{1}, 2));

    %% Sort by peak in odd laps
    [~, pkI2] = max(baseIncOdd_o2,[],2); [~, sI2] = sort(pkI2);
    [~, pkD2] = max(baseDecOdd_o2,[],2); [~, sD2] = sort(pkD2);
    [~, pkI3] = max(baseIncOdd_o3,[],2); [~, sI3] = sort(pkI3);
    [~, pkD3] = max(baseDecOdd_o3,[],2); [~, sD3] = sort(pkD3);

    col1 = omitColors{1};
    col2 = omitColors{2};

    %% Figure layout — manual axes positions
    figHandle = figure('Color', 'w', 'Position', [50 50 1200 1000]);

    % Row heights and positions (normalized) — heatmap rows 2.5x taller than trace rows
    hmH  = 0.20;   % heatmap height
    trH  = 0.08;   % trace height
    gap  = 0.03;   % gap between rows
    bot  = 0.06;   % bottom margin
    left = 0.06;   % left margin
    colW = 0.19;   % column width
    colG = 0.015;  % column gap
    cbW  = 0.015;  % colorbar width

    % Row bottoms from bottom up: row4(trace dec), row3(hm dec), row2(trace inc), row1(hm inc)
    r4bot = bot;
    r3bot = r4bot + trH + gap;
    r2bot = r3bot + hmH + gap;
    r1bot = r2bot + trH + gap;

    colL = left + (0:3) * (colW + colG);

    %% Helper functions
    function ax = makeHM(row, col, data, pb, lmMain, lmOther, lmMainCol, ylabStr, titleStr, titleCol, showCB)
        rowBots = [r1bot r2bot r3bot r4bot];
        rowHs   = [hmH   trH   hmH   trH  ];
        ax = axes(figHandle, 'Position', [colL(col), rowBots(row), colW, rowHs(row)]); 
        imagesc(ax, pb, 1:size(data,1), data);
        colormap(ax, flipud(gray)); clim(climAbs);
        hold(ax, 'on');
        xline(ax, lmMain,  '--', 'Color', lmMainCol,     'LineWidth', 1.5);
        xline(ax, lmOther, ':',  'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
        set(ax, 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out', 'FontName', targetFont, 'FontSize', 8);
        if ~isempty(ylabStr)
            ylabel(ax, ylabStr, 'FontName', targetFont, 'FontSize', 9);
        else
            set(ax, 'YTickLabel', '');
        end
        if ~isempty(titleStr)
            title(ax, titleStr, 'Color', titleCol, 'FontName', targetFont, 'FontSize', 10, 'FontWeight', 'bold');
        end
        if showCB
            cb = colorbar(ax, 'Position', [colL(col)+colW+0.005, rowBots(row), cbW, rowHs(row)]);
            cb.Ticks = climAbs; cb.TickLabels = {'low','high'};
            cb.TickDirection = 'out'; cb.Box = 'off';
            cb.Label.String = 'Activity (norm.)';
            cb.Label.FontName = targetFont; cb.Label.FontSize = 9;
        end
        if row == 3 || row == 4
            xlabel(ax, 'Position (cm)', 'FontName', targetFont, 'FontSize', 9);
        else
            set(ax, 'XTickLabel', '');
        end
    end

    function [ax1, ax2] = makeTrace(row, col, pb, baseM, omitM, diffM, lmMain, lmOther, lmCol, ylabStr, omitLab)
        rowBots = [r1bot r2bot r3bot r4bot];
        rowHs   = [hmH   trH   hmH   trH  ];
        ax1 = axes(figHandle, 'Position', [colL(col), rowBots(row), colW, rowHs(row)]); %#ok<LAXES>
        hold(ax1, 'on');
        plot(ax1, pb, mean(baseM,1,'omitnan'), 'k-', 'LineWidth', 1.2, 'DisplayName', 'Base');
        plot(ax1, pb, mean(omitM,1,'omitnan'), '-',  'LineWidth', 1.2, 'Color', lmCol, 'DisplayName', omitLab);
        xline(ax1, lmMain,  '--', 'Color', lmCol,         'LineWidth', 1.0);
        xline(ax1, lmOther, ':',  'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
        ylabel(ax1, ylabStr, 'FontName', targetFont, 'FontSize', 8);
        legend(ax1, 'Location', 'best', 'Box', 'off', 'FontSize', 6, 'FontName', targetFont);
        set(ax1, 'Box', 'off', 'TickDir', 'out', 'FontName', targetFont, 'FontSize', 8);
        if row == 4, xlabel(ax1, 'Position (cm)', 'FontName', targetFont, 'FontSize', 9); end

        ax2 = axes(figHandle, 'Position', [colL(col+1), rowBots(row), colW, rowHs(row)]); %#ok<LAXES>
        hold(ax2, 'on');
        plot(ax2, pb, mean(diffM,1,'omitnan'), '-', 'LineWidth', 1.2, 'Color', lmCol);
        yline(ax2, 0, 'k:', 'LineWidth', 1);
        xline(ax2, lmMain,  '--', 'Color', lmCol,         'LineWidth', 1.0);
        xline(ax2, lmOther, ':',  'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
        ylabel(ax2, '\Delta Activity', 'FontName', targetFont, 'FontSize', 8);
        title(ax2, ['\Delta ' omitLab], 'Color', lmCol, 'FontName', targetFont, 'FontSize', 9);
        set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontName', targetFont, 'FontSize', 8);
        if row == 4, xlabel(ax2, 'Position (cm)', 'FontName', targetFont, 'FontSize', 9); end
    end

    %% ROW 1 — Increase heatmaps
    makeHM(1, 1, baseInc_o2(sI2,:), posBins, 80, 120, [0 0 0], ...
        sprintf('Increase Omit2\n(n=%d)', numel(incIdx_o2)), 'Base', [0 0 0], false);
    makeHM(1, 2, omitInc_o2(sI2,:), posBins, 80, 120, col1, ...
        '', omitLabels{1}, col1, false);
    makeHM(1, 3, baseInc_o3(sI3,:), posBins, 120, 80, [0 0 0], ...
        sprintf('Increase Omit3\n(n=%d)', numel(incIdx_o3)), 'Base', [0 0 0], false);
    makeHM(1, 4, omitInc_o3(sI3,:), posBins, 120, 80, col2, ...
        '', omitLabels{2}, col2, true);

    %% ROW 2 — Increase mean traces
    makeTrace(2, 1, posBins, baseInc_o2, omitInc_o2, diffInc_o2, 80,  120, col1, '\DeltaF/F', omitLabels{1});
    makeTrace(2, 3, posBins, baseInc_o3, omitInc_o3, diffInc_o3, 120, 80,  col2, '\DeltaF/F', omitLabels{2});

    %% ROW 3 — Decrease heatmaps
    makeHM(3, 1, baseDec_o2(sD2,:), posBins, 80, 120, [0 0 0], ...
        sprintf('Decrease Omit2\n(n=%d)', numel(decIdx_o2)), 'Base', [0 0 0], false);
    makeHM(3, 2, omitDec_o2(sD2,:), posBins, 80, 120, col1, ...
        '', omitLabels{1}, col1, false);
    makeHM(3, 3, baseDec_o3(sD3,:), posBins, 120, 80, [0 0 0], ...
        sprintf('Decrease Omit3\n(n=%d)', numel(decIdx_o3)), 'Base', [0 0 0], false);
    makeHM(3, 4, omitDec_o3(sD3,:), posBins, 120, 80, col2, ...
        '', omitLabels{2}, col2, true);

    %% ROW 4 — Decrease mean traces
    makeTrace(4, 1, posBins, baseDec_o2, omitDec_o2, diffDec_o2, 80,  120, col1, '\DeltaF/F', omitLabels{1});
    makeTrace(4, 3, posBins, baseDec_o3, omitDec_o3, diffDec_o3, 120, 80,  col2, '\DeltaF/F', omitLabels{2});

    if ~isempty(p.Results.SavePath)
        saveFigureFormats(figHandle, p.Results.SavePath);
    end
end