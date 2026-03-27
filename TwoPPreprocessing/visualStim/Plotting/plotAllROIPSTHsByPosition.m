function plotAllROIPSTHsByPosition(sessionFileInfo, response, timeWindow, baselineWindow, peakTimeWindow)
% plotAllROIPSTHsByPosition  Loop over every ROI and save one page per ROI
%   into a single multi-page PDF. Each page tiles that ROI's normalized,
%   per-trial baseline-subtracted PSTHs across spatial positions.

%% Defaults
if nargin<3 || isempty(timeWindow)
    timeWindow = [0 2];
end
if nargin<4 || isempty(baselineWindow)
    baselineWindow = [-1.25 0];
end
if nargin<5 || isempty(peakTimeWindow)
    peakTimeWindow = [0.25 2]; 
end

psth   = response.psthData;
nPos   = numel(psth);
stimVs = vertcat(psth.stimValue);

% Extract unique coordinates
azimuth   = sort(unique(stimVs(:,1)), 'descend'); % Match visual field convention
elevation = sort(unique(stimVs(:,2)), 'ascend'); 
nAz = numel(azimuth);
nEl = numel(elevation);

% Prepare output folder & PDF
outDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
pdfPath = fullfile(outDir,'All_ROI_PSTH_byPosition.pdf');
if exist(pdfPath,'file'), delete(pdfPath); end

% Number of ROIs
nROI = size(psth(1).alignedResponses, 1);

%% Loop over each ROI
for roiIdx = 1:nROI
    fprintf('Processing ROI %d/%d...\n', roiIdx, nROI);
    
    % --- Step 1: Find best position for this ROI using per-trial subtraction ---
    peaks = nan(nPos,1);
    for k = 1:nPos
        tVec = psth(k).timeVector(:);
        
        % Robust extraction: Force [Time x Trials] even if nTrials = 1
        rawResp = psth(k).alignedResponses(roiIdx,:,:);
        data2d = reshape(rawResp, size(rawResp,2), size(rawResp,3)); 
        
        % Ensure masks are within actual data bounds to avoid indexing errors
        actualMin = min(tVec);
        actualMax = max(tVec);
        
        preMask = tVec >= max(baselineWindow(1), actualMin) & ...
                  tVec <= min(baselineWindow(2), actualMax);
              
        if ~any(preMask)
            % Fallback: if baseline is completely missing, use the first sample
            trialF0 = data2d(1, :);
        else
            trialF0 = mean(data2d(preMask, :), 1, 'omitnan');
        end
        
        zeroedTrials = data2d - trialF0; 
        meanTr = mean(zeroedTrials, 2, 'omitnan');
        
        postMask = tVec >= max(peakTimeWindow(1), actualMin) & ...
                   tVec <= min(peakTimeWindow(2), actualMax);
        
        if any(postMask)
            peaks(k) = max(meanTr(postMask), [], 'omitnan');
        else
            peaks(k) = -Inf;
        end
    end
    
    [absMax, bestIdx] = max(peaks);
    if isempty(absMax) || isnan(absMax) || absMax <= 0, absMax = 1; end
    
    % --- Step 2: Create figure for this ROI ---
    fig = figure('Visible','off', 'Color','w', 'Position',[100 100 1400 900]);
    til = tiledlayout(nEl, nAz, 'TileSpacing','compact', 'Padding','compact');
    sgtitle(sprintf('ROI %d PSTHs (Per-Trial Baseline Subtracted)', roiIdx), 'FontWeight','bold');
    
    % --- Step 3: Plot each spatial position ---
    for k = 1:nPos
        pos = psth(k).stimValue;
        col = find(azimuth == pos(1), 1);
        row = find(elevation == pos(2), 1);
        
        idx = (row-1)*nAz + col;
        ax = nexttile(til, idx); hold(ax,'on');
        
        tVec = psth(k).timeVector(:);
        rawResp = psth(k).alignedResponses(roiIdx,:,:);
        data2d = reshape(rawResp, size(rawResp,2), size(rawResp,3));
        
        % Per-trial baseline subtraction (Redoing for plotting)
        actualMin = min(tVec);
        actualMax = max(tVec);
        preMask = tVec >= max(baselineWindow(1), actualMin) & ...
                  tVec <= min(baselineWindow(2), actualMax);
              
        if ~any(preMask)
            trialF0 = data2d(1, :);
        else
            trialF0 = median(data2d(preMask, :), 1, 'omitnan');
        end
        zeroedTrials = data2d - trialF0; 
        
        meanTr = mean(zeroedTrials, 2, 'omitnan');
        semTr  = std(zeroedTrials, 0, 2, 'omitnan') ./ sqrt(size(zeroedTrials, 2));
        
        normMean = meanTr / absMax;
        normSEM  = semTr  / absMax;
        
        maskPlot = tVec >= max(baselineWindow(1), actualMin) & ...
                   tVec <= min(timeWindow(2), actualMax);
        
        t_plot = tVec(maskPlot);
        m_plot = normMean(maskPlot);
        s_plot = normSEM(maskPlot);
        
        if ~isempty(t_plot)
            % Shaded SEM
            fill(ax, [t_plot; flipud(t_plot)], [m_plot + s_plot; flipud(m_plot - s_plot)], ...
                [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.5);
            % Mean line
            plot(ax, t_plot, m_plot, 'k', 'LineWidth', 1.2);
        end
        
        % Reference lines
        xline(ax, 0, 'r--', 'LineWidth', 1); 
        yline(ax, 0, 'k:', 'Alpha', 0.3);    
        
        if k == bestIdx
            set(ax, 'LineWidth', 2, 'Box', 'on', 'XColor', 'r', 'YColor', 'r');
        end
        
        title(ax, sprintf('Az %d, El %.1f', pos(1), pos(2)), 'FontSize', 8);
        xlim(ax, [baselineWindow(1), timeWindow(2)]);
        ylim(ax, [-0.6, 1.2]); 
        set(ax, 'FontSize', 7);
        
        if col ~= 1, ylabel(ax, ''); else ylabel(ax, 'norm \DeltaF/F'); end
        if row ~= nEl, xlabel(ax, ''); else xlabel(ax, 'Time (s)'); end
    end
    
    % --- Step 4: Save page ---
    exportgraphics(fig, pdfPath, 'Append', true, 'ContentType','vector');
    close(fig);
end
fprintf('Processing complete. PDF saved: %s\n', pdfPath);
end