
%% model figs 
close all; clear axImg axStrip axAz axHold figHold patches_data

outputDirFig1 = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_7_Section4';
FOV_deg = 160;

% generate the landmark layout figure and extract patches
plotFinalFigures_Sonali(EXP_all, 'Texture-layout-38');
figSrc = gcf;

allAxes = findobj(figSrc, 'Type', 'axes');
azIdx = [];
for i = 1:numel(allAxes)
    t = get(get(allAxes(i),'Title'),'String');
    fprintf('%d: %s\n', i, t);
    if contains(lower(t), 'azimuth')
        azIdx = i;
    end
end
axAz = allAxes(azIdx);

figHold = figure('Visible','off');
axHold = axes('Parent', figHold);

kidsAll = get(axAz,'Children');
isPatch = strcmp(get(kidsAll,'Type'),'patch');
srcPatches = kidsAll(isPatch);
srcPatches = flipud(srcPatches(:));   % preserve original front-to-back order
patches_data = gobjects(numel(srcPatches),1);
for k = 1:numel(srcPatches)
    patches_data(k) = copyobj(srcPatches(k), axHold);
end
close(figSrc);

% build the real figure: axImg (right) with your styling applied
figBModel=figure('Color','w');
axImg = subplot(1,1,1); hold on;

img_right = imread(fullfile(outputDirFig1,'pos76.png'));
image(axImg, [-FOV_deg/2 FOV_deg/2], [1 0], img_right);
set(axImg,'YDir','normal');
xlim(axImg, [-FOV_deg/2 FOV_deg/2]); ylim(axImg, [0 1]);
xline(axImg, 0, 'k--', 'LineWidth', 1);              % straight-ahead reference line
xticks(axImg, -round(FOV_deg/2, -1) : 5 : round(FOV_deg/2, -1));
yticks(axImg, []); yticklabels(axImg, []);
xlabel(axImg, 'Visual azimuth (deg)');
title(axImg, 'View at corridor position = 76');
box(axImg,'off'); set(axImg,'TickDir','out');
grid(axImg,'on');
set(axImg, 'GridColor','w', 'GridLineStyle','--', 'GridAlpha',0.8, 'Layer','top');

% create the thin strip axes directly above axImg
drawnow;
imgPos = get(axImg, 'Position');
axStrip = axes('Position', [imgPos(1), imgPos(2)+imgPos(4)+0.01, imgPos(3), 0.08]);
hold(axStrip, 'on');

%
%  EDIT these two to control how much of each side you want shown 
includeRightSide = true;
includeLeftSide  = true;
rightMaxDeg = FOV_deg/2;   % e.g. 80 -- how far right to show
leftMaxDeg  = FOV_deg/2;   % e.g. 80 -- how far left to show (mirrored)

if includeRightSide
    for k = 1:numel(patches_data)
        copyobj(patches_data(k), axStrip);
    end
end

if includeLeftSide
    % Mirror each patch to the left side: negate XData so it lands on
    % the opposite side of 0 deg, since both corridor walls are symmetric.
    for k = 1:numel(patches_data)
        pMirror = copyobj(patches_data(k), axStrip);
        pMirror.XData = -pMirror.XData;
    end
end
close(figHold);


if includeRightSide && rightMaxDeg < FOV_deg/2
    fprintf('Note: rightMaxDeg < FOV_deg/2 -- adjust xlim below if you want asymmetric limits.\n');
end

xlim(axStrip, [-leftMaxDeg rightMaxDeg]); ylim(axStrip, [0 1]);
set(axStrip,'XTick',[]); set(axStrip,'YTick',[]);
set(axStrip,'Box','off'); set(axStrip,'XColor','none'); set(axStrip,'YColor','none');


fprintf('\n--- Patches in axStrip ---\n');
kids = get(axStrip, 'Children');
for k = 1:numel(kids)
    if strcmp(get(kids(k),'Type'), 'patch')
        fc = get(kids(k),'FaceColor');
        xd = get(kids(k),'XData');
        fprintf('color=[%.2f %.2f %.2f]  xrange=[%.2f %.2f]\n', fc, min(xd), max(xd));
    end
end

saveFigureFormats(figBModel, fullfile(outputDirFig1, 'corridor_inVisualAngle'));


%% load niave glm kerels 
close all; 
% load example session 
%EXP = load('SpkLin_EXP_M26004_20260318.mat');
plotFinalFigures_Sonali(EXPToTest, 'GLM-kernels', 260)
plotFinalFigures_Sonali(EXPToTest, 'GLM-kernels', [20 39])
figCDEF = gcf;
saveFigureFormats(figCDEF, fullfile(outputDirFig1, 'eg_nativeGLMKernel'));


EXP = load("C:\Users\sonali.sriranga\Desktop\V1\analyzed\M25131\20260318\SpkLin_EXP_M25131_20260318.mat");
%%
outputDirFig2 = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_8_Section4';
%
plotFinalFigures_Sonali(EXP_all, 'Resp-completeSnake-base-swap23-omit2-omit3') % any subset of conditions works
% plotFinalFigures_Sonali(EXP_all, 'Resp-completeSnake') % this will only include rois where all conditons were present
figA = gcf;

axAll = findall(figA, 'Type', 'Axes');
pos = cell2mat(get(axAll, 'Position'));
[~, sortIdx] = sortrows(pos, [-2 1]);
axAll = axAll(sortIdx);
pos = pos(sortIdx,:);
yvals = round(pos(:,2), 3);
uniqueY = unique(yvals, 'stable');

topY = uniqueY(1);
row1_axes = axAll(yvals == topY);

bottomY = uniqueY(end);
row3_axes = axAll(yvals == bottomY);
pos3 = cell2mat(get(row3_axes,'Position'));
[~, sIdx3] = sort(pos3(:,1));
row3_axes = row3_axes(sIdx3);

% 
ncond = numel(row1_axes) / 3;   % row1_axes = [DATA x ncond, VS x ncond, VSP x ncond]
assert(numel(row3_axes) == 1 + 2*ncond, ...
    'Unexpected number of bottom-row axes - check the figure layout assumptions.');

% FIX: strip any "(n=...)" suffix from the DATA-row titles before reusing
% them in Figure 2, since that n reflects the ALL-CELLS population (row 1,
% igroup==1), not the spatially-selective population used in Figure 2.
condTitles = cell(1,ncond);
for i = 1:ncond
    rawTitle = row1_axes(i).Title.String;
    condTitles{i} = regexprep(rawTitle, '\s*\(n=\d+\)', '');
end

% FIGURE 1
figAEdited = figure('Position', [100 100 1400 700]);
tlNew = tiledlayout(figAEdited, 3, ncond, 'TileSpacing','compact','Padding','compact');  
colormap(figAEdited, flipud(gray(256)));

for i = 1:numel(row1_axes)
    axNew = nexttile(tlNew);
    copyobj(allchild(row1_axes(i)), axNew);
    axNew.XLim = row1_axes(i).XLim;
    axNew.YLim = row1_axes(i).YLim;
    box(axNew,'off'); set(axNew,'TickDir','out');
    set(axNew, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
    set(axNew, 'YTick', [], 'YTickLabel',[]);
    axNew.Title.String = row1_axes(i).Title.String;

    % DATA row = first ncond tiles, regardless of ncond 
    if i <= ncond
        hold(axNew, 'on');
        for posBin = [40 80 120 160]
            xline(axNew, posBin, '--', 'Color', [1 1 1], 'LineWidth', 1);
        end
    end

    if i == numel(row1_axes)
        cb = colorbar(axNew);
        cb.TickDirection = 'out';
        cb.Box = 'off';
    end
end

allNewAxes1 = findall(figAEdited, 'Type', 'Axes');
set(allNewAxes1, 'CLim', [0 1]);

% FIGURE 2 
figKernelResid = figure('Position', [100 100 1600 500]);
tlResid = tiledlayout(figKernelResid, 2, ncond+1, 'TileSpacing','compact','Padding','compact');  % <-- ncond+1, not 7

axKernel = nexttile(tlResid, [2 1]);
kernelSrc = row3_axes(1);
copyobj(allchild(kernelSrc), axKernel);
axKernel.XLim = kernelSrc.XLim;
axKernel.YLim = kernelSrc.YLim;
axKernel.CLim = kernelSrc.CLim;
colormap(axKernel, RedWhiteBlue);
% FIX: append the ACTUAL spatially-selective n (from this panel's own row
% count / YLim), not anything copied from Figure 1's all-cells titles
nSpatialTotal = round(kernelSrc.YLim(2));
axKernel.Title.String = sprintf('%s (n=%d)', kernelSrc.Title.String, nSpatialTotal);
set(axKernel, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
set(axKernel, 'YTick', [], 'YTickLabel', []);
box(axKernel,'off'); set(axKernel,'TickDir','out');
hold(axKernel, 'on');
for posBin = [40 80 120 160]
    xline(axKernel, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end

%  sized based on ncond, not hardcoded 2:7/8:13
dataVSaxes  = row3_axes(2 : 1+ncond);
dataVSPaxes = row3_axes(2+ncond : 1+2*ncond);

for i = 1:ncond
    axNew = nexttile(tlResid);
    copyobj(allchild(dataVSaxes(i)), axNew);
    axNew.XLim = dataVSaxes(i).XLim;
    axNew.YLim = dataVSaxes(i).YLim;
    axNew.CLim = dataVSaxes(i).CLim;
    colormap(axNew, RedWhiteBlue);
    % FIX: compute the actual spatially-selective n for THIS panel,
    % rather than reusing the all-cells count from condTitles
    nSpatial_i = round(dataVSaxes(i).YLim(2));
    axNew.Title.String = sprintf('%s (n=%d)', condTitles{i}, nSpatial_i);
    set(axNew, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
    set(axNew, 'YTick', [], 'YTickLabel', []);
    box(axNew,'off'); set(axNew,'TickDir','out');
    if i == 1, ylabel(axNew, 'DATA - VS'); end
    if i == ncond
        cb = colorbar(axNew);
        cb.TickDirection = 'out';
        cb.Box = 'off';
    end
end

for i = 1:ncond
    axNew = nexttile(tlResid);
    copyobj(allchild(dataVSPaxes(i)), axNew);
    axNew.XLim = dataVSPaxes(i).XLim;
    axNew.YLim = dataVSPaxes(i).YLim;
    axNew.CLim = dataVSPaxes(i).CLim;
    colormap(axNew, RedWhiteBlue);
    % FIX: same, computed from THIS panel's own row count
    nSpatial_i2 = round(dataVSPaxes(i).YLim(2));
    axNew.Title.String = sprintf('%s (n=%d)', condTitles{i}, nSpatial_i2);
    set(axNew, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
    set(axNew, 'YTick', [], 'YTickLabel', []);
    box(axNew,'off'); set(axNew,'TickDir','out');
    if i == 1, ylabel(axNew, 'DATA - VSP'); end
    if i == ncond
        cb = colorbar(axNew);
        cb.TickDirection = 'out';
        cb.Box = 'off';
    end
end

sgtitle(tlResid, 'Spatially selective cells: spatial kernel & model residuals');

%
saveFigureFormats(figAEdited, fullfile(outputDirFig2, 'VS_VSP_resp_CompleteShake\VS_VSP_resp_CompleteShake_omit23_swap23'));
saveFigureFormats(figKernelResid, fullfile(outputDirFig2, 'Kernel_Residuals_CompleteSnake\Kernel_Residuals_CompleteSnake_omit23_swap23'));


%% make select examples to plot with the completesnake fig 

T = plotFinalFigures_Sonali(EXP_all, 'LLHi-w/oSpace');

% Pick the "spatial" example cell: spatially selective, strongest LLHrel
spatial_idx_pool = find(T.spatialcells);
[~, best_i] = max(T.LLHrel(spatial_idx_pool));
spatialCellIdx = spatial_idx_pool(best_i);

% Pick the "visual" example cell: good QC, vision-driven, not spatially selective
visual_idx_pool = find(T.goodcells & ~T.signicells);
[~, best_j] = max(T.LLHi_vis(visual_idx_pool));
visualCellIdx = visual_idx_pool(best_j);

fprintf('Spatial example cell: %s (idx %d, LLHrel=%.3f)\n', ...
    EXP_all.Spk.CellListString{spatialCellIdx}, spatialCellIdx, T.LLHrel(spatialCellIdx));
fprintf('Visual example cell:  %s (idx %d, LLHi_vis=%.3f)\n', ...
    EXP_all.Spk.CellListString{visualCellIdx}, visualCellIdx, T.LLHi_vis(visualCellIdx));

% Now plot both cells together
% plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', [spatialCellIdx, visualCellIdx]);
plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', [spatialCellIdx, visualCellIdx]);


% plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', {'M26004_20260318_cell#367'})
%% LLH VS vs VSP (scatter) 
plotFinalFigures_Sonali(EXP_all, 'LLHi-VSvsVSP-scatter')
figB = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figB, fullfile(outputDirFig2, 'LLH_improvement\LLH_improvement'));


%%
plotFinalFigures_Sonali(EXP_all, 'RF-LLHrel-summary')
figSummary = gcf;

axAll = findall(figSummary, 'Type', 'Axes');

for k = 1:numel(axAll)
    axSrc = axAll(k);

    % Use the title as both subfolder name and filename base
    titleStr = axSrc.Title.String;
    if isempty(titleStr)
        nameBase = sprintf('panel_%02d', k);
    else
        if iscell(titleStr), titleStr = strjoin(titleStr, '_'); end
        nameBase = matlab.lang.makeValidName(titleStr);  % sanitize for folder/file use
    end

    % Subfolder named after the panel, file inside named the same
    saveSubfolder = fullfile(outputDirFig2, nameBase);
    if ~exist(saveSubfolder, 'dir')
        mkdir(saveSubfolder);
    end

    % Copy this axes into its own standalone figure
    figB = figure;
    axNew = copyobj(axSrc, figB);
    axNew.Position = [0.15 0.15 0.75 0.75];

    % Re-attach legend if this panel had one
    oldLeg = get(axSrc, 'Legend');
    if ~isempty(oldLeg)
        legend(axNew, 'Location', oldLeg.Location);
        legend(axNew, 'boxoff');
    end

    defaultAxesProperties(axNew, 1);

    saveFigureFormats(figB, fullfile(saveSubfolder, nameBase));

    close(figB);
end


%% plot_position76_azimuth.m
% Two-panel figure:
close all 

plotFinalFigures_Sonali(EXP_all, 'Texture-layout-38')

leftImgPath  = 'Screenshot_2026-07-31_142304.png';
rightImgPath = "Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_7_Section4\pos76.png";

corridor_cm = 200;  
FOV_deg     = 160;   % <-- EDIT: total horizontal field of view spanned by rightImgPath

% position strip
figure('Color','w');
subplot(1,2,1); hold on;
img_left = imread(leftImgPath);
% Draws the strip with 0 cm at bottom, 200 cm at top.
% Flip the image argument order below if your PNG is stored upside down.
image([0 1], [corridor_cm 0], img_left);
set(gca,'YDir','normal');
xlim([0 1]); ylim([0 corridor_cm]);
set(gca,'XTick',[]);
ylabel('Position along corridor (cm)');
title('Corridor texture strip');
box off; set(gca,'TickDir','out');

 
subplot(1,2,2); hold on;
xticks()


img_right = imread(rightImgPath);
% Centre the image horizontally on 0 deg: spans -FOV/2 .. +FOV/2
image([-FOV_deg/2 FOV_deg/2], [1 0], img_right);
set(gca,'YDir','normal');
xlim([-FOV_deg/2 FOV_deg/2]); ylim([0 1]);
xline(0, 'k--', 'LineWidth', 1);   % straight-ahead reference line

xticks(-round(FOV_deg/2, -1) : 5 : round(FOV_deg/2, -1));
yticks([]); yticklabels([]);
xlabel('Visual azimuth (deg)');
title('View at corridor position = 76 %');
box off; set(gca,'TickDir','out');
grid on;
set(gca, 'GridColor', 'w', 'GridLineStyle', '--', 'GridAlpha', 0.8, 'Layer', 'top');



%% spatial componenets
close all 
% plotFinalFigures_Sonali(EXP_all, 'Ordered-Kernels-Space')

plotFinalFigures_Sonali(EXP_all, 'ExampleTraces-Kernels-Space')
% -> shows Space kernel ordered by Space
FigSpaceKer = gca;
FigSpaceKer = ancestor(FigSpaceKer, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigSpaceKer, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_10_Section4\spaceKernel_withEgs\', 'space_with_measured_kernal_egs'));


%% landmark and bg components % change to 10 exampels 
close all 
plotFinalFigures_Sonali(EXP_all, 'ScaledExamples-Kernels-Landmarks-BG')
FigLandBGKer = gca;
FigLandBGKer = ancestor(FigLandBGKer, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigLandBGKer, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_11_Section4\LandBGKernel_withEgs\', 'landmarks_BG_measured_kernal_egs'));

%% plot the 6 background segments 
% this plots the 6 bg segments starting from the very beginning - therefore
% remember to reorder them due to the offset. 
plotVisualFeatures


%% plot the revealed omitted background texture: l2 l3 l4 

plotVisualFeatures_omittedbackground

%% check for omit cells
plotFinalFigures_Sonali(EXP_rsp, 'OmitCells-summary')

%% look up cells
plotFinalFigures_Sonali(EXP_all, 'CellCategory-table');   % builds + saves the table
lookupCell([], 'M25131_20260318_cell#5');             % query one cell
T = CellCategoryTable;
T(T.IsSpatial & T.Animal==4, :)    
%% pca and cluserting test 
% shows Landmarks + BG + Data (measured response), all THREE panels,
%    ordered by Landmarks (the first tag) - NOT forced to Space
% %
% % testing cluserting : 
% 
% EXP = EXP_all;  % <-- your loaded batch variable
% 
% p_th = 0.05; LLHrel_th = 0.015;
% iPos = numel(EXP.GLMs{1}.Tuning);
% itex = 3;
% 
% goodcells  = EXP.Maps{1}.Tuning.nlogL_pval <= p_th & EXP.GLMs{1}.Tuning(itex).pval <= p_th;
% signicells = EXP.GLMs{1}.Tuning(iPos).pval <= p_th;
% 
% nlogL   = EXP.GLMs{1}.Perf.nlogL;
% nSpikes = EXP.GLMs{1}.Perf.nSpikes(:,1);
% LLHi_vis = -(nlogL(:,121) - nlogL(:,1)) ./ nSpikes / log(2);
% LLHi_pos = -(nlogL(:,end) - nlogL(:,1)) ./ nSpikes / log(2);
% LLHrel = LLHi_pos ./ LLHi_vis;
% goodLLHcells = LLHrel >= 1 + LLHrel_th;
% 
% spatial_mask = goodcells & signicells & goodLLHcells;
% 
% Space = squeeze(EXP.GLMs{1}.Tuning(iPos).meanrespModel(:,1,:));
% Xref  = Space(spatial_mask,:);
% 
% mu = mean(Xref, 2, 'omitnan');
% sd = std(Xref, 0, 2, 'omitnan');
% Xref_z = (Xref - mu) ./ sd;
% valid_idx = find(~isnan(mean(Xref_z,2,'omitnan')));
% Xvalid = Xref_z(valid_idx,:);
% 
% fprintf('N cells going into clustering: %d\n', numel(valid_idx));
% 
% %% =======================================================================
% %  (1) JUSTIFY NUMBER OF PCA COMPONENTS - variance explained
% %      (separate question from cluster count; do this FIRST)
% % =========================================================================
% [~, score, ~, ~, explained] = pca(Xvalid);
% cumulative_explained = cumsum(explained);
% 
% nShow = min(20, numel(explained));
% 
% figure('Name','PCA variance justification');
% subplot(1,2,1);
% plot(1:nShow, explained(1:nShow), 'o-', 'LineWidth', 1.2);
% xlabel('Principal component'); ylabel('% variance explained');
% title('Scree plot');
% box off; set(gca,'TickDir','out');
% 
% subplot(1,2,2);
% plot(1:nShow, cumulative_explained(1:nShow), 'o-', 'LineWidth', 1.2); hold on;
% yline(90, 'r--', '90%'); yline(95, 'k--', '95%');
% xlabel('Number of components'); ylabel('Cumulative % variance explained');
% title('Cumulative variance');
% box off; set(gca,'TickDir','out');
% sgtitle('PCA component justification');
% 
% % Printed summary: components needed to reach common thresholds
% n_for_80 = find(cumulative_explained >= 80, 1, 'first');
% n_for_90 = find(cumulative_explained >= 90, 1, 'first');
% n_for_95 = find(cumulative_explained >= 95, 1, 'first');
% 
% fprintf('\n--- PCA variance-explained summary ---\n');
% fprintf('Components needed for 80%% variance: %d\n', n_for_80);
% fprintf('Components needed for 90%% variance: %d\n', n_for_90);
% fprintf('Components needed for 95%% variance: %d\n', n_for_95);
% for i = 1:min(10, numel(explained))
%     fprintf('  PC%2d: %5.2f%% (cumulative %5.2f%%)\n', i, explained(i), cumulative_explained(i));
% end
% 
% % ---- SET n_PCs_justified here ----
% % Standard convention: retain components explaining >=90% cumulative
% % variance. 25 components meet this threshold.
% n_PCs_justified = n_for_90;
% fprintf('\nUsing n_PCs_justified = %d (%.1f%% variance explained, >=90%% threshold)\n', ...
%     n_PCs_justified, cumulative_explained(n_PCs_justified));
% 
% %% =========================================================================
% %  (2) JUSTIFY NUMBER OF CLUSTERS - elbow (CH) + silhouette, evaluated on
% %      the PCA-reduced space (PCA(4)), using OUR chosen distance/linkage.
% %
% %      NOTE: evalclusters(...,'linkage',...) does NOT respect custom
% %      distance/linkage settings - it silently uses its own internal
% %      defaults regardless of dist_metric/dist_method set above. This
% %      was caught because results were identical before/after switching
% %      to Ward/Euclidean. Fixed by building the tree ourselves with
% %      linkage(), then computing silhouette/CH manually at each K.
% % =========================================================================
% dist_metric = 'correlation';
% dist_method = 'average';
% 
% score_reduced = score(:, 1:n_PCs_justified);
% maxK = 20;
% 
% Z = linkage(score_reduced, dist_method, dist_metric);
% 
% Klist = 2:maxK;
% sil_scores = nan(size(Klist));
% ch_scores  = nan(size(Klist));
% 
% for i = 1:numel(Klist)
%     Ktest = Klist(i);
%     Ttest = cluster(Z, 'maxclust', Ktest);
%     if numel(unique(Ttest)) < 2
%         continue;  % skip degenerate cases
%     end
%     s = silhouette(score_reduced, Ttest, dist_metric);
%     sil_scores(i) = mean(s, 'omitnan');
%     try
%         ch_eva = evalclusters(score_reduced, Ttest, 'CalinskiHarabasz');
%         ch_scores(i) = ch_eva.CriterionValues;
%     catch
%         ch_scores(i) = NaN;
%     end
% end
% 
% figure('Name','Cluster count: Calinski-Harabasz');
% plot(Klist, ch_scores, 'o-', 'LineWidth', 1.2);
% xlabel('Number of clusters'); ylabel('Calinski-Harabasz score');
% title(sprintf('CH: %s dist, %s linkage, PCA(%d)', dist_metric, dist_method, n_PCs_justified));
% box off; set(gca,'TickDir','out');
% 
% figure('Name','Cluster count: Silhouette');
% plot(Klist, sil_scores, 'o-', 'LineWidth', 1.2);
% xlabel('Number of clusters'); ylabel('Silhouette score');
% title(sprintf('Silhouette: %s dist, %s linkage, PCA(%d)', dist_metric, dist_method, n_PCs_justified));
% box off; set(gca,'TickDir','out');
% 
% fprintf('\n--- Cluster count summary (%s / %s linkage) ---\n', dist_metric, dist_method);
% fprintf('Silhouette scores by K:\n');
% for i = 1:numel(Klist)
%     fprintf('  K=%2d: silhouette=%.4f, CH=%.2f\n', Klist(i), sil_scores(i), ch_scores(i));
% end
% 
% % ---- SET n_HClusters_justified here, after inspecting the CH/silhouette
% % plots above for this specific PCA(25)/correlation/average run ----
% n_HClusters_justified = 6;   % <-- UPDATE based on where the elbow/peak falls in THIS run's output
% fprintf('\nUsing n_HClusters_justified = %d\n', n_HClusters_justified);
% 
% %% =========================================================================
% %  (3) Quick visual sanity check at the chosen K, using PCA(4) space
% % =========================================================================
% Z = linkage(score_reduced, dist_method, dist_metric);
% T = cluster(Z, 'maxclust', n_HClusters_justified);
% 
% figure('Name','Example cluster mean shapes');
% nplot = n_HClusters_justified;
% for g = 1:nplot
%     subplot(ceil(nplot/2), 2, g);
%     idx = valid_idx(T == g);
%     M = Xref(idx,:);
%     mn = min(M,[],2,'omitnan'); mx = max(M,[],2,'omitnan');
%     Mn = (M - mn) ./ max(mx-mn, eps);   % normalize each cell before averaging
%     plot(mean(Mn,1,'omitnan'), 'LineWidth', 1.5);
%     title(sprintf('Group %d (n=%d)', g, numel(idx)));
%     box off; set(gca,'TickDir','out');
% end
% sgtitle(sprintf('K=%d, PCA(%d), %s dist, %s linkage', n_HClusters_justified, n_PCs_justified, dist_metric, dist_method));