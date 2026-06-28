function [suppressedROIs, suppressionTable, figHandle] = findLandmarkSuppressedROIs(allData, targetArea, varargin)
% findLandmarkSuppressedROIs
%
% WHAT THIS FUNCTION DOES:
% Identifies boutons/ROIs that are spatially suppressed at landmark positions
% during baseline laps. The approach is:
%
%   1. Load baseline LapActivity [nROIs x nLaps x nPosBins] for each session
%   2. Average across all laps to get mean spatial activity [nROIs x nPosBins]
%   3. Z-score each ROI's activity across position bins (mean=0, std=1)
%      This normalises for differences in baseline firing rate across ROIs
%      NOTE: swap zAct = allBaseAct if data is already z-scored across position
%   4. For each landmark position, compute mean z-score in a window (±windowCm)
%   5. An ROI is "suppressed" at a landmark if its z-score there < zThresh
%   6. Categorise ROIs by stimulus type:
%      - Grating only  (positions 40 and/or 120, NOT 80 or 160)
%      - Plaid only    (positions 80 and/or 160, NOT 40 or 120)
%      - Both types    (at least one grating AND one plaid position)
%
% OUTPUTS:
%   suppressedROIs   — struct with fields: gratingOnly, plaidOnly, bothTypes, any
%   suppressionTable — table: SessionIdx, ROI_idx, nLandmarksSuppressed,
%                      meanZ, per-landmark flags, stimType
%   figHandle        — figure with heatmaps and mean traces per category
%
% USAGE:
%   [suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPData, 'RSP');
%   [suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPData, 'RSP', ...
%       'zThresh', -0.5, 'windowCm', 10);

    p = inputParser;
    addRequired(p, 'allData',    @isstruct);
    addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DaysToPlot',     [1,2,3,4,5],    @isnumeric);
    addParameter(p, 'TypeToPlot',     'Boutons',       @(x) ischar(x)||isstring(x));
    addParameter(p, 'windowCm',       10,              @isnumeric);
    addParameter(p, 'zThresh',        -1,              @isnumeric);
    addParameter(p, 'landmarks',      [40 80 120 160], @isnumeric);
    addParameter(p, 'gratingPos',     [40 120],        @isnumeric);
    addParameter(p, 'plaidPos',       [80 160],        @isnumeric);
    addParameter(p, 'trackLength',    200,             @isnumeric);
    addParameter(p, 'ApplySmoothing', true,            @islogical);
    addParameter(p, 'FontName',       'Arial',         @(x) ischar(x)||isstring(x));
    addParameter(p, 'SavePath',       '',              @(x) ischar(x)||isstring(x));
    parse(p, allData, targetArea, varargin{:});

    windowCm    = p.Results.windowCm;
    zThresh     = p.Results.zThresh;
    landmarks   = p.Results.landmarks;
    gratingPos  = p.Results.gratingPos;
    plaidPos    = p.Results.plaidPos;
    trackLength = p.Results.trackLength;
    targetFont  = p.Results.FontName;
    daysToPlot  = p.Results.DaysToPlot;
%     daysToPlot(daysToPlot == 200) = [];
    cellType    = p.Results.TypeToPlot;

    if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
    if isfield(allData,'Type') && ~isfield(allData,'TypeImaged')
        [allData.TypeImaged] = allData.Type;
    elseif all(cellfun(@isempty, {allData.TypeImaged}))
        [allData.TypeImaged] = deal(char(cellType));
    end

    %% baseline activity across all sessions
    allBaseAct    = [];
    allSessionIDs = [];
    allROIIDs     = [];

    for d = 1:length(daysToPlot)
        day = daysToPlot(d);
        if day == 5
            dayMask = ([allData.Day] == 5 | [allData.Day] == 200);
        else
            dayMask = ([allData.Day] == day);
        end

        daySessions = allData(dayMask & ...
            strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
            strcmpi(string({allData.TypeImaged}), string(cellType)));
        if isempty(daySessions), continue; end

        sessGlobalIdx = find(dayMask & ...
            strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
            strcmpi(string({allData.TypeImaged}), string(cellType)));

        for s = 1:length(daySessions)
            sess = daySessions(s);
            if ~isfield(sess, 'ConditionData') || isempty(sess.ConditionData), continue; end

            condNames = fieldnames(sess.ConditionData);
            baseIdx   = find(contains(lower(condNames), 'baseline') | contains(lower(condNames), 'default'), 1);
            if isempty(baseIdx), baseIdx = 1; end
            baseName  = condNames{baseIdx};

            if isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
                idx = sess.FilteredROIs;
            else
                idx = 1:size(sess.ConditionData.(baseName).LapActivity, 1);
            end

            lapAct = sess.ConditionData.(baseName).LapActivity(idx, :, :);
            if p.Results.ApplySmoothing, lapAct = smoothLapActivity(lapAct); end

            %  Mean across laps [nROI x nBins]
            baseAct = squeeze(mean(lapAct, 2, 'omitnan'));
            if size(baseAct,1) ~= numel(idx), baseAct = baseAct'; end

            allBaseAct    = [allBaseAct;    baseAct];                                  
            allSessionIDs = [allSessionIDs; repmat(sessGlobalIdx(s), numel(idx), 1)];
            allROIIDs     = [allROIIDs;     idx(:)];                                   
        end
    end

    if isempty(allBaseAct)
        warning('No baseline data found.');
        suppressedROIs = []; suppressionTable = table(); figHandle = [];
        return;
    end

    %%Z-score each ROI across position
    nBins   = size(allBaseAct, 2);
    posBins = linspace(0, trackLength, nBins);
    nROIs   = size(allBaseAct, 1);

    muAct  = mean(allBaseAct, 2, 'omitnan');
    stdAct = std(allBaseAct,  0, 2, 'omitnan');
    stdAct(stdAct == 0) = 1;
%     zAct   = (allBaseAct - muAct) ./ stdAct;  % [nROI x nBins]
    % NOTE: swap to line below if data already z-scored across position:
    zAct = allBaseAct;

    %% Check suppression at each landmark
    nLandmarks   = numel(landmarks);
    isSuppressed = false(nROIs, nLandmarks);

    for li = 1:nLandmarks
        lm      = landmarks(li);
        winMask = posBins >= (lm - windowCm) & posBins <= (lm + windowCm);
        zAtLM   = mean(zAct(:, winMask), 2, 'omitnan');
        isSuppressed(:, li) = zAtLM < zThresh;
    end

    % Which columns correspond to grating vs plaid landmarks
    gratingCols = ismember(landmarks, gratingPos);
    plaidCols   = ismember(landmarks, plaidPos);

    % suppressed at any grating / any plaid
    suppAnyGrating = any(isSuppressed(:, gratingCols), 2);
    suppAnyPlaid   = any(isSuppressed(:, plaidCols),   2);
    suppAnything   = any(isSuppressed, 2);

    %% Categorise by stimulus type
    isGratingOnly = suppAnyGrating & ~suppAnyPlaid;
    isPlaidOnly   = suppAnyPlaid   & ~suppAnyGrating;
    isBothTypes   = suppAnyGrating &  suppAnyPlaid;
    nSuppressed   = sum(isSuppressed, 2);

    %% 
    fprintf('\n--- Landmark-Suppressed ROIs (z < %.1f, window ±%d cm) ---\n', zThresh, windowCm);
    fprintf('  Total ROIs:                          %d\n', nROIs);
    fprintf('  Any suppression:                     %d (%.1f%%)\n', sum(suppAnything), 100*sum(suppAnything)/nROIs);
    fprintf('  ------------------------------------------------\n');
    fprintf('  Grating only  (40 and/or 120):       %d (%.1f%%)\n', sum(isGratingOnly), 100*sum(isGratingOnly)/nROIs);
    fprintf('    of which suppressed at pos 40:     %d\n', sum(isSuppressed(:, landmarks==40)));
    fprintf('    of which suppressed at pos 120:    %d\n', sum(isSuppressed(:, landmarks==120)));
    fprintf('    of which suppressed at both:       %d\n', sum(isSuppressed(:,landmarks==40) & isSuppressed(:,landmarks==120)));
    fprintf('  Plaid only    (80 and/or 160):       %d (%.1f%%)\n', sum(isPlaidOnly), 100*sum(isPlaidOnly)/nROIs);
    fprintf('    of which suppressed at pos 80:     %d\n', sum(isSuppressed(:, landmarks==80)));
    fprintf('    of which suppressed at pos 160:    %d\n', sum(isSuppressed(:, landmarks==160)));
    fprintf('    of which suppressed at both:       %d\n', sum(isSuppressed(:,landmarks==80) & isSuppressed(:,landmarks==160)));
    fprintf('  Both types    (grating + plaid):     %d (%.1f%%)\n', sum(isBothTypes), 100*sum(isBothTypes)/nROIs);
    fprintf('  ------------------------------------------------\n');
    fprintf('  Breakdown of Both by position count:\n');
    for n = 2:4
        fprintf('    Suppressed at %d positions:        %d\n', n, sum(isBothTypes & nSuppressed == n));
    end

    %% 
    suppressedROIs.gratingOnly = find(isGratingOnly);
    suppressedROIs.plaidOnly   = find(isPlaidOnly);
    suppressedROIs.bothTypes   = find(isBothTypes);
    suppressedROIs.any         = find(suppAnything);

    %% Mean z-score at all landmark positions per ROI
    lmBinMask = false(1, nBins);
    for li = 1:nLandmarks
        winMask   = posBins >= (landmarks(li) - windowCm) & posBins <= (landmarks(li) + windowCm);
        lmBinMask = lmBinMask | winMask;
    end
    meanZatLandmarks = mean(zAct(:, lmBinMask), 2, 'omitnan');

    % Stimulus type label per ROI
    stimType = repmat({'none'}, nROIs, 1);
    stimType(isGratingOnly) = {'gratingOnly'};
    stimType(isPlaidOnly)   = {'plaidOnly'};
    stimType(isBothTypes)   = {'both'};

    suppressionTable = table(...
        allSessionIDs, allROIIDs, nSuppressed, meanZatLandmarks, stimType, ...
        isSuppressed(:,1), isSuppressed(:,2), isSuppressed(:,3), isSuppressed(:,4), ...
        'VariableNames', {'SessionIdx','ROI_idx','nLandmarksSuppressed','meanZ_atLandmarks','stimType', ...
                          'suppLM_40','suppLM_80','suppLM_120','suppLM_160'});

    %%  3 columns x 2 rows (heatmap + trace)
    groups      = {find(isGratingOnly), find(isPlaidOnly), find(isBothTypes)};
    groupNames  = {'Grating only (40 and/or 120)', 'Plaid only (80 and/or 160)', 'Both types (grating + plaid)'};
    groupColors = {[0.2 0.5 0.8], [0.8 0.4 0.1], [0.4 0.2 0.6]};
    nGroups     = 3;

    % Sort each group by mean z-score at landmarks
    for g = 1:nGroups
        gIdx = groups{g};
        if ~isempty(gIdx)
            [~, sortOrd] = sort(meanZatLandmarks(gIdx), 'ascend');
            groups{g}    = gIdx(sortOrd);
        end
    end

    figHandle = figure('Color', 'w', 'Position', [50 50 1100 820]);

    hmH  = 0.36; trH = 0.17;
    bot1 = 0.50; bot2 = 0.10;
    colW = 0.26; colG = 0.04;
    colL = 0.07 + (0:nGroups-1) * (colW + colG);

    gratingColor = [0.2 0.5 0.8];
    plaidColor   = [0.8 0.4 0.1];

    for col = 1:nGroups
        gIdx = groups{col};
        zM   = zAct(gIdx, :);
        gCol = groupColors{col};

        %% Heatmap
        axHM = axes('Position', [colL(col), bot1, colW, hmH]);

        if isempty(zM)
            text(axHM, 0.5, 0.5, 'No boutons', 'HorizontalAlignment', 'center', ...
                'FontName', targetFont, 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
            title(axHM, groupNames{col}, 'Color', gCol, 'FontName', targetFont, 'FontSize', 10, 'FontWeight', 'bold');
            axis(axHM, 'off');
        else
            imagesc(axHM, posBins, 1:size(zM,1), zM);
            colormap(axHM, flipud(gray));
            clim([-2 2]);
            hold(axHM, 'on');
            for lm = landmarks
                if ismember(lm, gratingPos)
                    xline(axHM, lm, '--', 'Color', gratingColor, 'LineWidth', 1.5);
                else
                    xline(axHM, lm, '--', 'Color', plaidColor,   'LineWidth', 1.5);
                end
            end
            set(axHM, 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out', ...
                'XTick', landmarks, 'FontName', targetFont, 'FontSize', 9);
            title(axHM, sprintf('%s\n(n=%d)', groupNames{col}, numel(gIdx)), ...
                'Color', gCol, 'FontName', targetFont, 'FontSize', 10, 'FontWeight', 'bold');
            if col == 1
                ylabel(axHM, 'Boutons', 'FontName', targetFont, 'FontSize', 10);
            else
                set(axHM, 'YTickLabel', '');
            end
            if col == nGroups
                cb = colorbar(axHM, 'eastoutside');
                cb.Label.String  = 'Z-score';
                cb.TickDirection = 'out'; cb.Box = 'off';
                cb.FontName      = targetFont; cb.FontSize = 8;
            end
        end

        %% Mean trace
        axTR = axes('Position', [colL(col), bot2, colW, trH]);

        if ~isempty(zM)
            hold(axTR, 'on');
            muZ  = mean(zM, 1, 'omitnan');
            semZ = std(zM,  0, 1, 'omitnan') / sqrt(size(zM,1));
            fill(axTR, [posBins, fliplr(posBins)], [muZ+semZ, fliplr(muZ-semZ)], ...
                gCol, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            plot(axTR, posBins, muZ, '-', 'Color', gCol, 'LineWidth', 1.5);
            yline(axTR, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            for lm = landmarks
                if ismember(lm, gratingPos)
                    xline(axTR, lm, '--', 'Color', gratingColor, 'LineWidth', 1.5);
                else
                    xline(axTR, lm, '--', 'Color', plaidColor,   'LineWidth', 1.5);
                end
            end
            set(axTR, 'Box', 'off', 'TickDir', 'out', 'XLim', [0 trackLength], ...
                'XTick', landmarks, 'FontName', targetFont, 'FontSize', 9);
            xlabel(axTR, 'Position (cm)', 'FontName', targetFont, 'FontSize', 10);
            if col == 1
                ylabel(axTR, 'Mean Z-score ± SEM', 'FontName', targetFont, 'FontSize', 10);
            else
                set(axTR, 'YTickLabel', '');
            end
        end
    end

    %% axes in bottom right
    axLeg = axes('Position', [0.70, 0.01, 0.28, 0.06], 'Parent', figHandle); 
    hold(axLeg, 'on');
    h1 = plot(axLeg, nan, nan, '--', 'Color', gratingColor, 'LineWidth', 1.5);
    h2 = plot(axLeg, nan, nan, '--', 'Color', plaidColor,   'LineWidth', 1.5);
    legend(axLeg, [h1 h2], {'Grating landmarks (40, 120)', 'Plaid landmarks (80, 160)'}, ...
        'Orientation', 'horizontal', 'Box', 'off', ...
        'FontName', targetFont, 'FontSize', 8);
    axis(axLeg, 'off');

    sgtitle(sprintf('Landmark-Suppressed ROIs — %s  |  z < %.1f  |  window ±%d cm', ...
        targetArea, zThresh, windowCm), ...
        'FontName', targetFont, 'FontSize', 12, 'FontWeight', 'bold');

    if ~isempty(p.Results.SavePath)
        saveFigureFormats(figHandle, p.Results.SavePath);
    end
end



% function [suppressedROIs, suppressionTable, figHandle] = findLandmarkSuppressedROIs(allData, targetArea, varargin)
% % findLandmarkSuppressedROIs
% %
% % WHAT THIS FUNCTION DOES:
% % Identifies boutons/ROIs that are spatially suppressed at landmark positions
% % during baseline trials. The approach is:
% %
% % OUTPUTS:
% %   suppressedROIs   — struct with fields: exact1, exact2, exact3, exact4
% %   suppressionTable — table: SessionIdx, ROI_idx, nLandmarksSuppressed, meanZ, per-landmark flags
% %   figHandle        — figure with heatmaps and mean traces per category
% %
% % :
% %   [suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPData, 'RSP');
% %   [suppROIs, suppTable, fig] = findLandmarkSuppressedROIs(RSPData, 'RSP', ...
% %       'zThresh', -0.5, 'windowCm', 10);
% 
%     p = inputParser;
%     addRequired(p, 'allData',    @isstruct);
%     addRequired(p, 'targetArea', @(x) ischar(x) || isstring(x));
%     addParameter(p, 'DaysToPlot',     [1,2,3,4,5],    @isnumeric);
%     addParameter(p, 'TypeToPlot',     'Boutons',       @(x) ischar(x)||isstring(x));
%     addParameter(p, 'windowCm',       10,              @isnumeric);
%     addParameter(p, 'zThresh',        -1,              @isnumeric);
%     addParameter(p, 'landmarks',      [40 80 120 160], @isnumeric);
%     addParameter(p, 'trackLength',    200,             @isnumeric);
%     addParameter(p, 'ApplySmoothing', true,            @islogical);
%     addParameter(p, 'FontName',       'Arial',         @(x) ischar(x)||isstring(x));
%     addParameter(p, 'SavePath',       '',              @(x) ischar(x)||isstring(x));
%     parse(p, allData, targetArea, varargin{:});
% 
%     windowCm    = p.Results.windowCm;
%     zThresh     = p.Results.zThresh;
%     landmarks   = p.Results.landmarks;
%     trackLength = p.Results.trackLength;
%     targetFont  = p.Results.FontName;
%     daysToPlot  = p.Results.DaysToPlot;
% %     daysToPlot(daysToPlot == 200) = [];
%     cellType    = p.Results.TypeToPlot;
% 
%     if all(cellfun(@isempty, {allData.TargetArea})), [allData.TargetArea] = deal(char(targetArea)); end
%     if isfield(allData,'Type') && ~isfield(allData,'TypeImaged')
%         [allData.TypeImaged] = allData.Type;
%     elseif all(cellfun(@isempty, {allData.TypeImaged}))
%         [allData.TypeImaged] = deal(char(cellType));
%     end
% 
%     %% baseline activity across all sessions
%     allBaseAct    = [];
%     allSessionIDs = [];
%     allROIIDs     = [];
% 
%     for d = 1:length(daysToPlot)
%         day = daysToPlot(d);
%         if day == 5
%             dayMask = ([allData.Day] == 5 | [allData.Day] == 200);
%         else
%             dayMask = ([allData.Day] == day);
%         end
% 
%         daySessions = allData(dayMask & ...
%             strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
%             strcmpi(string({allData.TypeImaged}), string(cellType)));
%         if isempty(daySessions), continue; end
% 
%         sessGlobalIdx = find(dayMask & ...
%             strcmpi(string({allData.TargetArea}), string(targetArea)) & ...
%             strcmpi(string({allData.TypeImaged}), string(cellType)));
% 
%         for s = 1:length(daySessions)
%             sess = daySessions(s);
%             if ~isfield(sess, 'ConditionData') || isempty(sess.ConditionData), continue; end
% 
%             condNames = fieldnames(sess.ConditionData);
%             baseIdx   = find(contains(lower(condNames), 'baseline') | contains(lower(condNames), 'default'), 1);
%             if isempty(baseIdx), baseIdx = 1; end
%             baseName  = condNames{baseIdx};
% 
%             if isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
%                 idx = sess.FilteredROIs;
%             else
%                 idx = 1:size(sess.ConditionData.(baseName).LapActivity, 1);
%             end
% 
%             lapAct = sess.ConditionData.(baseName).LapActivity(idx, :, :);
%             if p.Results.ApplySmoothing, lapAct = smoothLapActivity(lapAct); end
% 
%             % Mean across laps [nROI x nBins]
%             baseAct = squeeze(mean(lapAct, 2, 'omitnan'));
%             if size(baseAct,1) ~= numel(idx), baseAct = baseAct'; end
% 
%             allBaseAct    = [allBaseAct;    baseAct];                                  
%             allSessionIDs = [allSessionIDs; repmat(sessGlobalIdx(s), numel(idx), 1)];
%             allROIIDs     = [allROIIDs;     idx(:)];                                   
%         end
%     end
% 
%     if isempty(allBaseAct)
%         warning('No baseline data found.');
%         suppressedROIs = []; suppressionTable = table(); figHandle = [];
%         return;
%     end
% 
%     %%  Z-score each ROI across position
%     nBins   = size(allBaseAct, 2);
%     posBins = linspace(0, trackLength, nBins);
%     nROIs   = size(allBaseAct, 1);
% 
%     muAct  = mean(allBaseAct, 2, 'omitnan');
%     stdAct = std(allBaseAct,  0, 2, 'omitnan');
%     stdAct(stdAct == 0) = 1;
%     zAct   = (allBaseAct - muAct) ./ stdAct;  % [nROI x nBins]
%     % Replace zAct with the raw data
%     %zAct = allBaseAct;  % already z-scored from preprocessing
% 
%     %%  Check suppression at each landmark
%     nLandmarks   = numel(landmarks);
%     isSuppressed = false(nROIs, nLandmarks);
% 
%     for li = 1:nLandmarks
%         lm      = landmarks(li);
%         winMask = posBins >= (lm - windowCm) & posBins <= (lm + windowCm);
%         zAtLM   = mean(zAct(:, winMask), 2, 'omitnan');
%         isSuppressed(:, li) = zAtLM < zThresh;
%     end
% 
%     %%  Categorise by number of suppressed landmarks
%     nSuppressed = sum(isSuppressed, 2);  % [nROIs x 1]
% 
%     isExact1 = nSuppressed == 1;
%     isExact2 = nSuppressed == 2;
%     isExact3 = nSuppressed == 3;
%     isExact4 = nSuppressed == 4;
%     isAny1p  = nSuppressed >= 1;
%     isAny2p  = nSuppressed >= 2;
%     isAny3p  = nSuppressed >= 3;
% 
%     %% 
%     fprintf('\n--- Landmark-Suppressed ROIs (z < %.1f, window ±%d cm) ---\n', zThresh, windowCm);
%     fprintf('  Total ROIs:                        %d\n',   nROIs);
%     fprintf('  ------------------------------------------------\n');
%     fprintf('  Suppressed at exactly 1 landmark:  %d (%.1f%%)\n', sum(isExact1), 100*sum(isExact1)/nROIs);
%     fprintf('  Suppressed at exactly 2 landmarks: %d (%.1f%%)\n', sum(isExact2), 100*sum(isExact2)/nROIs);
%     fprintf('  Suppressed at exactly 3 landmarks: %d (%.1f%%)\n', sum(isExact3), 100*sum(isExact3)/nROIs);
%     fprintf('  Suppressed at exactly 4 landmarks: %d (%.1f%%)\n', sum(isExact4), 100*sum(isExact4)/nROIs);
%     fprintf('  ------------------------------------------------\n');
%     fprintf('  Suppressed at >= 1 landmark:       %d (%.1f%%)\n', sum(isAny1p), 100*sum(isAny1p)/nROIs);
%     fprintf('  Suppressed at >= 2 landmarks:      %d (%.1f%%)\n', sum(isAny2p), 100*sum(isAny2p)/nROIs);
%     fprintf('  Suppressed at >= 3 landmarks:      %d (%.1f%%)\n', sum(isAny3p), 100*sum(isAny3p)/nROIs);
%     fprintf('  Suppressed at all 4 landmarks:     %d (%.1f%%)\n', sum(isExact4), 100*sum(isExact4)/nROIs);
% 
%     %% Output struct
%     suppressedROIs.exact1 = find(isExact1);
%     suppressedROIs.exact2 = find(isExact2);
%     suppressedROIs.exact3 = find(isExact3);
%     suppressedROIs.exact4 = find(isExact4);
%     suppressedROIs.any2p  = find(isAny2p);
%     suppressedROIs.any3p  = find(isAny3p);
% 
%     %% Mean z-score at landmark positions per ROI
%     lmBinMask        = false(1, nBins);
%     for li = 1:nLandmarks
%         lm          = landmarks(li);
%         winMask     = posBins >= (lm - windowCm) & posBins <= (lm + windowCm);
%         lmBinMask   = lmBinMask | winMask;
%     end
%     meanZatLandmarks = mean(zAct(:, lmBinMask), 2, 'omitnan');
% 
%     suppressionTable = table(...
%         allSessionIDs, allROIIDs, nSuppressed, meanZatLandmarks, ...
%         isSuppressed(:,1), isSuppressed(:,2), isSuppressed(:,3), isSuppressed(:,4), ...
%         'VariableNames', {'SessionIdx','ROI_idx','nLandmarksSuppressed','meanZ_atLandmarks', ...
%                           'suppLM1','suppLM2','suppLM3','suppLM4'});
% 
%     %% Figure — 5 columns (exact1, exact2, exact3, exact4) x 2 rows (heatmap + trace)
%     groups     = {find(isExact1), find(isExact2), find(isExact3), find(isExact4)};
%     groupNames = {'Exactly 1', 'Exactly 2', 'Exactly 3', 'All 4'};
%     nGroups    = 4;
% 
%     % Sort each group by mean z-score at landmarks (most suppressed first)
%     for g = 1:nGroups
%         gIdx = groups{g};
%         if ~isempty(gIdx)
%             [~, sortOrd] = sort(meanZatLandmarks(gIdx), 'ascend');
%             groups{g}    = gIdx(sortOrd);
%         end
%     end
% 
%     figHandle = figure('Color', 'w', 'Position', [50 50 1400 800]);
% 
%     hmH  = 0.38; trH = 0.18;
%     bot1 = 0.52; bot2 = 0.08;
%     colW = 0.20; colG = 0.025;
%     colL = 0.06 + (0:nGroups-1) * (colW + colG);
% 
%     for col = 1:nGroups
%         gIdx = groups{col};
%         zM   = zAct(gIdx, :);
% 
%         %% Heatmap
%         axHM = axes('Position', [colL(col), bot1, colW, hmH]);
% 
%         if isempty(zM)
%             text(axHM, 0.5, 0.5, 'No boutons', 'HorizontalAlignment', 'center', ...
%                 'FontName', targetFont, 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
%             title(axHM, groupNames{col}, 'FontName', targetFont, 'FontSize', 11, 'FontWeight', 'bold');
%             axis(axHM, 'off');
%         else
%             imagesc(axHM, posBins, 1:size(zM,1), zM);
%             colormap(axHM, flipud(gray));
%             clim([-2 2]);
%             hold(axHM, 'on');
%             for lm = landmarks
%                 xline(axHM, lm, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
%             end
%             set(axHM, 'YDir', 'normal', 'Box', 'off', 'TickDir', 'out', ...
%                 'XTick', landmarks, 'FontName', targetFont, 'FontSize', 9);
%             title(axHM, sprintf('%s\n(n=%d)', groupNames{col}, numel(gIdx)), ...
%                 'FontName', targetFont, 'FontSize', 10, 'FontWeight', 'bold');
%             if col == 1
%                 ylabel(axHM, 'Boutons', 'FontName', targetFont, 'FontSize', 10);
%             else
%                 set(axHM, 'YTickLabel', '');
%             end
%             if col == nGroups
%                 cb = colorbar(axHM, 'eastoutside');
%                 cb.Label.String  = 'Z-score';
%                 cb.TickDirection = 'out'; cb.Box = 'off';
%                 cb.FontName      = targetFont; cb.FontSize = 8;
%             end
%         end
% 
%         %% Mean trace
%         axTR = axes('Position', [colL(col), bot2, colW, trH]);
% 
%         if ~isempty(zM)
%             hold(axTR, 'on');
%             muZ  = mean(zM, 1, 'omitnan');
%             semZ = std(zM,  0, 1, 'omitnan') / sqrt(size(zM,1));
% 
%             fill(axTR, [posBins, fliplr(posBins)], [muZ+semZ, fliplr(muZ-semZ)], ...
%                 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
%             plot(axTR, posBins, muZ, 'k-', 'LineWidth', 1.5);
%             yline(axTR, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
% 
%             for lm = landmarks
%                 xline(axTR, lm, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
%             end
% 
%             set(axTR, 'Box', 'off', 'TickDir', 'out', 'XLim', [0 trackLength], ...
%                 'XTick', landmarks, 'FontName', targetFont, 'FontSize', 9);
%             xlabel(axTR, 'Position (cm)', 'FontName', targetFont, 'FontSize', 10);
%             if col == 1
%                 ylabel(axTR, 'Mean Z-score ± SEM', 'FontName', targetFont, 'FontSize', 10);
%             else
%                 set(axTR, 'YTickLabel', '');
%             end
%         end
%     end
% 
%     sgtitle(sprintf('Landmark-Suppressed ROIs — %s  |  z < %.1f  |  window ±%d cm', ...
%         targetArea, zThresh, windowCm), ...
%         'FontName', targetFont, 'FontSize', 12, 'FontWeight', 'bold');
% 
%     if ~isempty(p.Results.SavePath)
%         saveFigureFormats(figHandle, p.Results.SavePath);
%     end
% end