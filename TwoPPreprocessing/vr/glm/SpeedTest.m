%% Load selected ROIs from their respective session files, plot 4 panels each, and save individually
selectedROIs = {
     'M25132_20260228', 140,   'early_nospeedcorr';   % keep
    'M26003_20260325', 83,   'early_speedcorr';  % also kept 
 %   'M25132_20260228', 17,   'tbd';
%     'M25132_20260228', 6,   'early_nospeedcorr';
%    'M25132_20260228', 17,  'early_speedcorr';
  %  'M25132_20260228', 31,  'late_nospeedcorr';
 %    'M25132_20260228', 217, 'late_speedcorr';
%     'M25132_20260313', 197, 'early_nospeedcorr';
%    'M25132_20260313', 199, 'early_speedcorr';
%     'M25132_20260313', 51,  'late_nospeedcorr';
 %   'M25132_20260313', 108,  'late_speedcorr';
%     'M26003_20260325', 297, 'speedcorr';
   % 'M26003_20260325', 278, 'speedcorr';
    'M26003_20260324', 39, 'speedcorr_neg'; 
    %'M25132_20260228', 107, 'check' 
   % 'M26003_20260324', 58, 'speedcorr';
   'M26003_20260325', 4, 'speedcorr'; % KEEP 
    %'M26003_20260324', 58, 'speedcorr';
};

sessionFiles = containers.Map();
sessionFiles('M25132_20260228') = "Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260228\M25132_20260228_Response_M25132_LandManipCorridor_20260228_CombinedRuns.mat";
sessionFiles('M25132_20260313') = "Z:\ibn-vision\DATA\SUBJECTS\M25132\Analysis\20260313\M25132_20260313_Response_M25132_LandManipCorridor_20260313_CombinedRuns.mat";
sessionFiles('M26003_20260325') = "Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat";
sessionFiles('M26003_20260324') = "Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns.mat";

outDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.11\speedExamples\';
if ~exist(outDir, 'dir'), mkdir(outDir); end

landmarks = [40 80 120 160];
landmark_halfwidth = 4;  % spans 8cm total, centred on each landmark
window_cm = 10; % for mean speed
local_window_cm = 15; % for peak search 
zoom_halfwidth_cm = 40; % for plot 

uniqueSessions = unique(selectedROIs(:,1), 'stable');

for s = 1:numel(uniqueSessions)
    sessName = uniqueSessions{s};
    fprintf('Loading %s...\n', sessName);
    response = load(sessionFiles(sessName));

    data = response.lapPositionActivity.spks;
    speedData = response.lapPositionRunningSpeed;
    baseline_idx = response.trialIndicesByCondition.Baseline;

    theseRows = find(strcmp(selectedROIs(:,1), sessName));
    for ri = 1:numel(theseRows)
        cellnum = selectedROIs{theseRows(ri), 2};
        category = selectedROIs{theseRows(ri), 3};

        cell_spks_all = squeeze(data(cellnum, :, :));
        cell_spks = cell_spks_all(baseline_idx, :);
        cell_speed = speedData(baseline_idx, :);

        xpos = linspace(0, 200, size(cell_spks, 2));
        bin_width = xpos(2) - xpos(1);
        n_laps = numel(baseline_idx);

        mean_trace = mean(cell_spks, 1, 'omitnan');
        [peakVal, peak_bin] = max(mean_trace);
        peak_pos = xpos(peak_bin);
        [~, lmIdx] = min(abs(peak_pos - landmarks));
        lm = landmarks(lmIdx);

        window_bins = max(1, round(window_cm / bin_width));
        pre_start = max(1, peak_bin - window_bins);
        pre_end = max(pre_start, peak_bin - 1);
        speed_per_lap = mean(cell_speed(:, pre_start:pre_end), 2, 'omitnan');

        local_window_bins = round(local_window_cm / bin_width);
        lo = max(1, peak_bin - local_window_bins);
        hi = min(size(cell_spks, 2), peak_bin + local_window_bins);
        [~, local_peak_bin_rel] = max(cell_spks(:, lo:hi), [], 2, 'omitnan');
        local_peak_bin = local_peak_bin_rel + lo - 1;
        local_peak_pos = xpos(local_peak_bin)';

        [r, p] = corr(speed_per_lap, local_peak_pos, 'Type', 'Spearman', 'rows', 'complete');

        [sorted_speed, sortIdx] = sort(speed_per_lap);
        sorted_spks = cell_spks(sortIdx, :);

        %% 
        figA = figure('Color','w','Position',[100 100 1500 350]);
        tl = tiledlayout(figA, 1, 12, 'TileSpacing','compact','Padding','compact');
        sgtitle(sprintf('%s cell#%d — dist=%+.1fcm, r=%.2f, p=%.3f', ...
            sessName, cellnum, peak_pos-lm, r, p), 'Interpreter','none');

        ax1 = nexttile(tl, [1 4]);
        imagesc(ax1, xpos, 1:n_laps, cell_spks);
        colormap(ax1, flipud(gray));
        set(ax1, 'YDir', 'normal'); hold(ax1, 'on');
        xlim(ax1, [xpos(1) xpos(end)]); ylim(ax1, [1 n_laps]);
        addLandmarkBands(ax1, landmarks, landmark_halfwidth, 'x');
        xline(ax1, peak_pos, 'r-', 'LineWidth', 1);
        ylabel(ax1, 'Laps (natural order)'); xlabel(ax1, 'Position (cm)');
        title(ax1, 'Full raster (Baseline)');

        ax2 = nexttile(tl, [1 4]);
        imagesc(ax2, xpos, 1:n_laps, sorted_spks);
        colormap(ax2, flipud(gray));
        set(ax2, 'YDir', 'normal'); hold(ax2, 'on');
        xlim(ax2, [peak_pos-zoom_halfwidth_cm, peak_pos+zoom_halfwidth_cm]);
        ylim(ax2, [1 n_laps]);
        addLandmarkBands(ax2, landmarks, landmark_halfwidth, 'x');
        xline(ax2, peak_pos, 'r-', 'LineWidth', 1.5);
        ylabel(ax2, 'Laps (sorted by speed)'); xlabel(ax2, 'Position (cm)');
        title(ax2, 'Raster (Baseline, sorted by speed)');

        ax3 = nexttile(tl, [1 2]);
        plot(ax3, sorted_speed, 1:n_laps, 'b-', 'LineWidth', 1); hold(ax3,'on');
        scatter(ax3, sorted_speed, 1:n_laps, 30, 'b', 'filled');
        ylim(ax3, [1 n_laps]);
        xlabel(ax3, sprintf('Speed (cm/s)\n(%dcm preceding peak)', window_cm));
        title(ax3, 'Speed per lap');

        ax4 = nexttile(tl, [1 2]); hold(ax4,'on');
        scatter(ax4, speed_per_lap, local_peak_pos, 30, 'filled', 'MarkerFaceAlpha', 0.6);
        addLandmarkBands(ax4, landmarks, landmark_halfwidth, 'y');
        xlabel(ax4, 'Speed (cm/s)'); ylabel(ax4, 'Local peak position (cm)');
        title(ax4, sprintf('r=%.2f, p=%.3f', r, p));
        

        % Instead of relying on imagesc's default color scaling:
        clim_max = prctile(cell_spks(:), 97);   % adjust percentile as needed
        set(ax1, 'CLim', [0, clim_max]);
        set(ax2, 'CLim', [0, clim_max]);
        %% Apply axes properties to EVERY axes in the figure
        allAxes = findall(figA, 'Type', 'Axes');
        for a = 1:numel(allAxes)
            defaultAxesProperties(allAxes(a), 1);
        end

        %% Save individually, matching your naming convention
        saveName = sprintf('%s_#%d_%s', sessName, cellnum, category);
        saveFigureFormats(figA, fullfile(outDir, saveName));
        close(figA);
        fprintf('  Saved %s\n', saveName);
    end
end

fprintf('\nAll %d selected ROIs saved to %s\n', size(selectedROIs,1), outDir);

function addLandmarkBands(ax, landmarks, halfwidth, orientation)
% Draws a translucent band, centred on each landmark, spanning
% landmark +/- halfwidth. orientation='x': vertical bands (position on
% x-axis, e.g. rasters). orientation='y': horizontal bands (position on
% y-axis, e.g. the scatter panel).
    xl = xlim(ax); yl = ylim(ax);
    for L = landmarks
        if strcmp(orientation, 'x')
            if L+halfwidth < xl(1) || L-halfwidth > xl(2), continue; end
            patch(ax, [L-halfwidth, L+halfwidth, L+halfwidth, L-halfwidth], ...
                      [yl(1), yl(1), yl(2), yl(2)], [0.2 0.4 0.9], ...
                      'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        else
            if L+halfwidth < yl(1) || L-halfwidth > yl(2), continue; end
            patch(ax, [xl(1), xl(2), xl(2), xl(1)], ...
                      [L-halfwidth, L-halfwidth, L+halfwidth, L+halfwidth], [0.2 0.4 0.9], ...
                      'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
end



% %% Fast population-level speed-confound check across RF=-5 candidates
% % Reuses the exact per-cell logic already validated on #108/#367/#209,
% % just looped, so it stays fast and simple rather than implementing the
% % full grid-search delay estimator.
% RFpos = EXP_all.GLMs{1}.Perf.bestRFpos(:);
% rfm5_idx = find(T.goodcells & RFpos == -5);
% candidateCellNames = EXP_all.Spk.CellListString(rfm5_idx);
% fprintf('%d RF=-5 candidates\n', numel(candidateCellNames));
% results = table();
% 
% for i = 1:numel(candidateCellNames)   % your ranked RF=-5 list from earlier
%     thisName = candidateCellNames{i};
%     parts = strsplit(thisName, '_');
%     mousenumber = parts{1};
%     sessionName = parts{2};
%     neuronIdx = str2double(regexprep(parts{3}, 'cell#', ''));
% 
%     try
%         infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
%         loadedInfo = load(infoPath, 'sessionFileInfo');
%         sessionFileInfo = loadedInfo.sessionFileInfo;
%         stimNames = {sessionFileInfo.stimFiles.name};
%         targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);
%         if isempty(targetIdx), continue; end
% 
%         response = load(sessionFileInfo.stimFiles(targetIdx).Response);
%         data = response.lapPositionActivity.spks;
%         speedData = response.lapPositionRunningSpeed;
%         baseline_idx = response.trialIndicesByCondition.Baseline;
% 
%         cell_spks_all = squeeze(data(neuronIdx, :, :));
%         cell_spks = cell_spks_all(baseline_idx, :);
%         cell_speed = speedData(baseline_idx, :);
% 
%         xpos = linspace(0, 200, size(cell_spks, 2));
%         bin_width = xpos(2) - xpos(1);
% 
%         mean_trace = mean(cell_spks, 1, 'omitnan');
%         [~, peak_bin] = max(mean_trace);
%         peak_pos = xpos(peak_bin);
% 
%         window_bins = max(1, round(10 / bin_width));
%         pre_start = max(1, peak_bin - window_bins);
%         speed_per_lap = mean(cell_speed(:, pre_start:max(pre_start,peak_bin-1)), 2, 'omitnan');
% 
%         local_bins = round(30 / bin_width);
%         lo = max(1, peak_bin - local_bins);
%         hi = min(size(cell_spks, 2), peak_bin + local_bins);
%         [~, local_bin_rel] = max(cell_spks(:, lo:hi), [], 2, 'omitnan');
%         local_peak_pos = xpos(local_bin_rel + lo - 1)';
% 
%         [r, p] = corr(speed_per_lap, local_peak_pos, 'Type', 'Spearman', 'rows', 'complete');
% 
%         results = [results; table({thisName}, peak_pos, r, p, numel(baseline_idx), ...
%             'VariableNames', {'Cell','PeakPos','R','P','nLaps'})];
%     catch
%         continue
%     end
% end
% 
% nSig = sum(results.P < 0.05);
% fprintf('%d/%d (%.0f%%) RF=-5 candidates show a significant speed-position correlation (p<0.05)\n', ...
%     nSig, height(results), 100*nSig/height(results));
% disp(results);
% 
% 
%% Speed-confound check for EVERY cell in this session, exported as one PDF
% Adds RF position pulled from EXP_all — ASSUMPTION FLAGGED BELOW, since
% I can't verify this lookup matches without you confirming EXP_all is
% loaded and uses the same 'MOUSE_SESSION_cell#NUM' naming convention
% seen elsewhere tonight (e.g. EXP_all.Spk.CellListString).

% response = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat");
% data = response.lapPositionActivity.spks;           % [nCells x nLaps x nBins]
% speedData = response.lapPositionRunningSpeed;       % [nLaps x nBins]
% baseline_idx = response.trialIndicesByCondition.Baseline;
% 
% mousenumber = 'M26003';
% sessionName = '20260325_2';
% nCells = size(data, 1);
% 
% xpos = linspace(0, 200, size(data, 3));
% bin_width = xpos(2) - xpos(1);
% n_laps = numel(baseline_idx);
% 
% window_cm = 10;
% window_bins = max(1, round(window_cm / bin_width));
% local_window_cm = 30;
% local_window_bins = round(local_window_cm / bin_width);
% 
% pdfName = sprintf('%s_%s_speedConfoundCheck_allCells.pdf', mousenumber, sessionName);
% if isfile(pdfName), delete(pdfName); end   % start fresh, since exportgraphics APPENDS
% 
% resultsTable = table();
% 
% for cellnum = 1:nCells
%     try
%         cell_spks_all = squeeze(data(cellnum, :, :));
%         cell_spks = cell_spks_all(baseline_idx, :);
%         cell_speed = speedData(baseline_idx, :);
% 
%         mean_trace = mean(cell_spks, 1, 'omitnan');
%         [peakVal, peak_bin] = max(mean_trace);
%         peak_pos = xpos(peak_bin);
% 
%         pre_start = max(1, peak_bin - window_bins);
%         pre_end = max(pre_start, peak_bin - 1);
%         speed_per_lap = mean(cell_speed(:, pre_start:pre_end), 2, 'omitnan');
% 
%         lo = max(1, peak_bin - local_window_bins);
%         hi = min(size(cell_spks, 2), peak_bin + local_window_bins);
%         [~, local_peak_bin_rel] = max(cell_spks(:, lo:hi), [], 2, 'omitnan');
%         local_peak_bin = local_peak_bin_rel + lo - 1;
%         local_peak_pos = xpos(local_peak_bin)';
% 
%         [r, p] = corr(speed_per_lap, local_peak_pos, 'Type', 'Spearman', 'rows', 'complete');
% 
%         % ---- RF position lookup — ASSUMPTION: EXP_all is loaded in this
%         % workspace and uses this exact naming convention. If not found,
%         % RFpos_this is left as NaN and the page still generates. ----
%         RFpos_this = NaN;
%         if exist('EXP_all', 'var')
%             thisName = sprintf('%s_%s_cell#%d', mousenumber, sessionName, cellnum);
%             matchIdx = find(strcmp(EXP_all.Spk.CellListString, thisName), 1);
%             if ~isempty(matchIdx)
%                 RFpos_this = EXP_all.GLMs{1}.Perf.bestRFpos(matchIdx);
%             end
%         end
% 
%         resultsTable = [resultsTable; table(cellnum, peak_pos, r, p, RFpos_this, ...
%             'VariableNames', {'Cell','PeakPos','R','P','RFpos'})];
% 
%         %% ---- Plotting ----
%         [sorted_speed, sortIdx] = sort(speed_per_lap);
%         sorted_spks = cell_spks(sortIdx, :);
% 
%         fig = figure('Color','w','Position',[100 100 1200 350], 'Visible','off');
%         tl = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
%         sgtitle(sprintf('Cell #%d — RF = %.1f deg — speed-sorted raster & retro/prospective check', ...
%             cellnum, RFpos_this));
% 
%         ax1 = nexttile(tl);
%         imagesc(ax1, xpos, 1:n_laps, sorted_spks);
%         colormap(ax1, flipud(gray));
%         set(ax1, 'YDir', 'normal'); hold(ax1, 'on');
%         xline(ax1, peak_pos, 'r-', 'LineWidth', 1);
%         xlim(ax1, [max(xpos(1),peak_pos-50), min(xpos(end),peak_pos+40)]);
%         ylabel(ax1, 'Laps (sorted by speed)');
%         xlabel(ax1, 'Position (cm)');
%         title(ax1, 'Raster (Baseline, sorted by speed)');
%         box(ax1, 'off'); set(ax1, 'TickDir', 'out');
% 
%         ax2 = nexttile(tl); hold(ax2, 'on');
%         plot(ax2, sorted_speed, 1:n_laps, 'k.-', 'MarkerSize', 10);
%         ylim(ax2, [1 n_laps]);
%         xlabel(ax2, sprintf('Speed (cm/s)\n(%dcm preceding peak)', window_cm));
%         title(ax2, 'Speed per lap');
%         box(ax2, 'off'); set(ax2, 'TickDir', 'out');
% 
%         ax3 = nexttile(tl); hold(ax3, 'on');
%         scatter(ax3, speed_per_lap, local_peak_pos, 30, 'filled', 'MarkerFaceAlpha', 0.6);
%         xlabel(ax3, 'Speed (cm/s)');
%         ylabel(ax3, 'Local peak position (cm)');
%         title(ax3, sprintf('r=%.2f, p=%.3f', r, p));
%         box(ax3, 'off'); set(ax3, 'TickDir', 'out');
% 
%         exportgraphics(fig, pdfName, 'Append', true);
%         close(fig);
% 
%         fprintf('Cell #%d done (%d/%d)\n', cellnum, cellnum, nCells);
%     catch ME
%         fprintf('Cell #%d skipped — %s\n', cellnum, ME.message);
%         continue
%     end
% end
% 
% fprintf('\nSaved all pages to %s\n', pdfName);
% disp(resultsTable);