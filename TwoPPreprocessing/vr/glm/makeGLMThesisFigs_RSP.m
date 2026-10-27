
%% model figs 
close all; clear axImg axStrip axAz axHold figHold patches_data

outputDirFig1 = 'ZZ:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.9\';
EXP_all = load("C:\Users\sonali.sriranga\Desktop\RSP\batch\SpkLinRSP_batch_SpkLinRSP_VisSpkLinParams.mat"); 
EXPToTest = load("C:\Users\sonali.sriranga\Desktop\RSP\analyzed\M25132\20260228\SpkLin_EXP_M25132_20260228.mat");
outputDirFig2 = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.9\';
%% load niave glm kerels 
close all; 
% load example session 
% %EXP = load('SpkLin_EXP_M26004_20260318.mat');
% plotFinalFigures_Sonali(EXPToTest, 'GLM-kernels', 260)
% plotFinalFigures_Sonali(EXPToTest, 'GLM-kernels', [20 39])
% figCDEF = gcf;
%saveFigureFormats(figCDEF, fullfile(outputDirFig1, 'eg_nativeGLMKernel'));
%%
% this is the version of the figure in the the thesis:
plotFinalFigures_Sonali_JF(EXP_all, 'Resp-snake-base-omit2')

% this is the older version of the figure
% plotFinalFigures_Sonali(EXP_all, 'Resp-completeSnake-base-swap23-omit2-omit3') % any subset of conditions works
% plotFinalFigures_Sonali(EXP_all, 'Resp-completeSnake') % this will only include rois where all conditons were present
figA = gcf;

axAll = findall(figA, 'Type', 'Axes');
pos = cell2mat(get(axAll, 'Position'));
[~, sortIdx] = sortrows(pos, [-2 1]);
axAll = axAll(sortIdx);
pos = pos(sortIdx,:);
yvals = round(pos(:,2), 3);
uniqueY = unique(yvals, 'stable');

% FIX: explicitly identify all three rows rather than only top/bottom.
% Row 1 (top)    = "All cells" raw responses -- titles ("base"/"omit2") live here
% Row 2 (middle) = "Spatially selective cells" raw responses -- the data we
%                  actually want for Figure 1, but with no titles of its own
% Row 3 (bottom) = spatial kernels + DATA-VS / DATA-VSP residuals
assert(numel(uniqueY) == 3, ...
    'Expected exactly 3 rows of axes in the source figure - check the figure layout assumptions.');

topY = uniqueY(1);
row1_axes = axAll(yvals == topY);      % All cells (titles live here)

midY = uniqueY(2);
row2_axes = axAll(yvals == midY);      % Spatially selective cells, raw responses
pos2 = cell2mat(get(row2_axes, 'Position'));
[~, sIdx2] = sort(pos2(:,1));
row2_axes = row2_axes(sIdx2);          % sort left-to-right, same as row3 below

bottomY = uniqueY(end);
row3_axes = axAll(yvals == bottomY);   % Kernels + residuals
pos3 = cell2mat(get(row3_axes,'Position'));
[~, sIdx3] = sort(pos3(:,1));
row3_axes = row3_axes(sIdx3);

% 
ncond = numel(row1_axes) / 3;   % row1_axes = [DATA x ncond, VS x ncond, VSP x ncond]
assert(numel(row2_axes) == 3*ncond, ...
    'Unexpected number of middle-row axes - check the figure layout assumptions.');
assert(numel(row3_axes) == 1 + 2*ncond, ...
    'Unexpected number of bottom-row axes - check the figure layout assumptions.');

% FIX: strip any "(n=...)" suffix from the DATA-row titles before reusing
% them in Figure 1/2, since that n reflects the ALL-CELLS population (row 1),
% not the spatially-selective population these figures actually plot.
condTitles = cell(1,ncond);
for i = 1:ncond
    rawTitle = row1_axes(i).Title.String;
    condTitles{i} = regexprep(rawTitle, '\s*\(n=\d+\)', '');
end

% FIGURE 1 -- now built from ROW 2 (spatially selective raw responses),
% using row1's titles for condition labels since row2 has none of its own.
figAEdited = figure('Position', [100 100 1400 700]);
tlNew = tiledlayout(figAEdited, 3, ncond, 'TileSpacing','compact','Padding','compact');  
colormap(figAEdited, flipud(gray(256)));

for i = 1:numel(row2_axes)
    axNew = nexttile(tlNew);
    copyobj(allchild(row2_axes(i)), axNew);
    axNew.XLim = row2_axes(i).XLim;
    axNew.YLim = row2_axes(i).YLim;
    box(axNew,'off'); set(axNew,'TickDir','out');
    set(axNew, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
    set(axNew, 'YTick', [], 'YTickLabel',[]);
    
    % Title text doesn't exist on row2 -- reuse the (n-stripped) condition
    % label from row1, and append THIS row's own spatially-selective n.
    condIdx = mod(i-1, ncond) + 1;
    nSpatial_i = round(row2_axes(i).YLim(2));
    axNew.Title.String = sprintf('%s (n=%d)', condTitles{condIdx}, nSpatial_i);

    % DATA row = first ncond tiles, regardless of ncond 
    if i <= ncond
        hold(axNew, 'on');
        for posBin = [40 80 120 160]
            xline(axNew, posBin, '--', 'Color', [1 1 1], 'LineWidth', 1);
        end
    end

    if i == numel(row2_axes)
        cb = colorbar(axNew);
        cb.TickDirection = 'out';
        cb.Box = 'off';
    end
end

allNewAxes1 = findall(figAEdited, 'Type', 'Axes');
set(allNewAxes1, 'CLim', [0 1]);

% FIGURE 2  (unchanged -- already correctly used row3)
figKernelResid = figure('Position', [100 100 1200 700]);
tlResid = tiledlayout(figKernelResid, 2, ncond+1, 'TileSpacing','compact','Padding','compact');  

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

% saveFigureFormats(figAEdited, fullfile(outputDirFig2, 'VS_VSP_resp_Shake\VS_VSP_resp_snake_base_omit2'));
% saveFigureFormats(figKernelResid, fullfile(outputDirFig2, 'Kernel_Residuals_Snake\Kernel_Residuals_Snake_base-omit2'));

saveFigureFormats(figAEdited, fullfile(outputDirFig2, 'SuppFig_VS_VSP_resp_Shake\VS_VSP_resp_snake_allcondition'));
saveFigureFormats(figKernelResid, fullfile(outputDirFig2, 'SuppFig_Kernel_Residuals_Snake\Kernel_Residuals_Snake_allcondition'));

%
%
%% V2 of the above figure; this version is used in the thesis 



% this is the version of the figure in the the thesis:
plotFinalFigures_Sonali_JF(EXP_all, 'Resp-snake-base-omit2')

figA = gcf;

axAll = findall(figA, 'Type', 'Axes');
pos = cell2mat(get(axAll, 'Position'));
[~, sortIdx] = sortrows(pos, [-2 1]);
axAll = axAll(sortIdx);
pos = pos(sortIdx,:);
yvals = round(pos(:,2), 3);
uniqueY = unique(yvals, 'stable');

topY = uniqueY(1);
row1_axes = axAll(yvals == topY);      % All cells (igroup==1)

midY = uniqueY(2);
row2_axes = axAll(yvals == midY);      % Spatially selective (igroup==2)

bottomY = uniqueY(end);
row3_axes = axAll(yvals == bottomY);
pos3 = cell2mat(get(row3_axes,'Position'));
[~, sIdx3] = sort(pos3(:,1));
row3_axes = row3_axes(sIdx3);

ncond = numel(row1_axes) / 3;   % row1_axes = [DATA x ncond, VS x ncond, VSP x ncond]
assert(numel(row2_axes) == numel(row1_axes), ...
    'Row 1 and row 2 axis counts do not match - check figure layout assumptions.');
assert(numel(row3_axes) == 1 + 2*ncond, ...
    'Unexpected number of bottom-row axes - check the figure layout assumptions.');

% FIX: strip any "(n=...)" suffix from the DATA-row titles before reusing
% them below, since that n reflects the ALL-CELLS population (row 1,
% igroup==1), not the spatially-selective population used elsewhere.
condTitles = cell(1,ncond);
for i = 1:ncond
    rawTitle = row1_axes(i).Title.String;
    condTitles{i} = regexprep(rawTitle, '\s*\(n=\d+\)', '');
end

% ===================== FIGURE 1: ALL CELLS =====================
figAllCells = figure('Position', [100 100 1400 700]);
tlAll = tiledlayout(figAllCells, 3, ncond, 'TileSpacing','compact','Padding','compact');
colormap(figAllCells, flipud(gray(256)));

for i = 1:numel(row1_axes)
    axNew = nexttile(tlAll);
    copyobj(allchild(row1_axes(i)), axNew);
    axNew.XLim = row1_axes(i).XLim;
    axNew.YLim = row1_axes(i).YLim;
    box(axNew,'off'); set(axNew,'TickDir','out');
    set(axNew, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
    set(axNew, 'YTick', [], 'YTickLabel',[]);
    axNew.Title.String = row1_axes(i).Title.String;

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

allAxesAllCells = findall(figAllCells, 'Type', 'Axes');
set(allAxesAllCells, 'CLim', [0 1]);
sgtitle(tlAll, 'All cells');

% FIGURE 1b: SPATIALLY SELECTIVE CELLS 
% Matched to the same population used for the kernel/residual figure
% below, so Figure 1b and Figure 2 describe the same neurons.
figSpatialCells = figure('Position', [100 100 1400 700]);
tlSpatial = tiledlayout(figSpatialCells, 3, ncond, 'TileSpacing','compact','Padding','compact');
colormap(figSpatialCells, flipud(gray(256)));

for i = 1:numel(row2_axes)
    axNew = nexttile(tlSpatial);
    copyobj(allchild(row2_axes(i)), axNew);
    axNew.XLim = row2_axes(i).XLim;
    axNew.YLim = row2_axes(i).YLim;
    box(axNew,'off'); set(axNew,'TickDir','out');
    set(axNew, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
    set(axNew, 'YTick', [], 'YTickLabel',[]);
    axNew.Title.String = row2_axes(i).Title.String;

    if i <= ncond
        hold(axNew, 'on');
        for posBin = [40 80 120 160]
            xline(axNew, posBin, '--', 'Color', [1 1 1], 'LineWidth', 1);
        end
    end

    if i == numel(row2_axes)
        cb = colorbar(axNew);
        cb.TickDirection = 'out';
        cb.Box = 'off';
    end
end

allAxesSpatial = findall(figSpatialCells, 'Type', 'Axes');
set(allAxesSpatial, 'CLim', [0 1]);
sgtitle(tlSpatial, 'Spatially selective cells');

% row3_axes was already the spatially-selective population (built once,
% after the igroup loop, from ss_idx_sorted), so this is unchanged.
figKernelResid = figure('Position', [100 100 1200 700]);
tlResid = tiledlayout(figKernelResid, 2, ncond+1, 'TileSpacing','compact','Padding','compact');

axKernel = nexttile(tlResid, [2 1]);
kernelSrc = row3_axes(1);
copyobj(allchild(kernelSrc), axKernel);
axKernel.XLim = kernelSrc.XLim;
axKernel.YLim = kernelSrc.YLim;
axKernel.CLim = kernelSrc.CLim;
colormap(axKernel, RedWhiteBlue);
nSpatialTotal = round(kernelSrc.YLim(2));
axKernel.Title.String = sprintf('%s (n=%d)', kernelSrc.Title.String, nSpatialTotal);
set(axKernel, 'XTick', [40 80 120 160], 'XTickLabel', {'40','80','120','160'});
set(axKernel, 'YTick', [], 'YTickLabel', []);
box(axKernel,'off'); set(axKernel,'TickDir','out');
hold(axKernel, 'on');
for posBin = [40 80 120 160]
    xline(axKernel, posBin, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
end

dataVSaxes  = row3_axes(2 : 1+ncond);
dataVSPaxes = row3_axes(2+ncond : 1+2*ncond);

for i = 1:ncond
    axNew = nexttile(tlResid);
    copyobj(allchild(dataVSaxes(i)), axNew);
    axNew.XLim = dataVSaxes(i).XLim;
    axNew.YLim = dataVSaxes(i).YLim;
    axNew.CLim = dataVSaxes(i).CLim;
    colormap(axNew, RedWhiteBlue);
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


saveFigureFormats(figAllCells,    fullfile(outputDirFig2, 'VS_VSP_resp_Shake-allcells\VS_VSP_resp_snake_base_omit2'));
saveFigureFormats(figSpatialCells, fullfile(outputDirFig2, 'VS_VSP_resp_Shake\VS_VSP_resp_snake_base_omit2'));
saveFigureFormats(figKernelResid, fullfile(outputDirFig2, 'Kernel_Residuals_Snake\Kernel_Residuals_Snake_base_omit2'));


%% make select examples to plot with the completesnake fig - use this version for thesis 
T = plotFinalFigures_Sonali(EXP_all, 'LLHi-w/oSpace');


% Pick the top 10 "spatial" example cells: spatially selective, strongest LLHrel
spatial_idx_pool = find(T.spatialcells);
[~, sortIdx] = sort(T.LLHrel(spatial_idx_pool), 'descend');
spatialCellIdx = spatial_idx_pool(sortIdx(1:10));

% Pick the top 10 "visual" example cells: good QC, vision-driven, not spatially selective
visual_idx_pool = find(T.goodcells & ~T.signicells);
[~, sortIdx2] = sort(T.LLHi_vis(visual_idx_pool), 'descend');
visualCellIdx = visual_idx_pool(sortIdx2(1:10));

for k = 1:10
    fprintf('Spatial example cell %d: %s (idx %d, LLHrel=%.3f)\n', k, ...
        EXP_all.Spk.CellListString{spatialCellIdx(k)}, spatialCellIdx(k), T.LLHrel(spatialCellIdx(k)));
end
for k = 1:10
    fprintf('Visual example cell %d: %s (idx %d, LLHi_vis=%.3f)\n', k, ...
        EXP_all.Spk.CellListString{visualCellIdx(k)}, visualCellIdx(k), T.LLHi_vis(visualCellIdx(k)));
end

% Now plot each group of 10 cells (two separate calls to the same existing function)
%plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', spatialCellIdx);
plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', {'M26003_20260325_cell#8', 'M26003_20260324_cell#29', 'M25132_20260228_cell#104'});
FigResoSingleCells = gca;
FigResoSingleCells = ancestor(FigResoSingleCells, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigResoSingleCells, fullfile(outputDirFig2, 'VS_VSP_Data_examplefits\VS_VSP_Data_examplefits_allconditions-SPATIALCELLS'));

close all; 

%plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', visualCellIdx); % run
%this to select examples 
plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', {'M25132_20260228_cell#50', 'M26003_20260325_cell#278', 'M26003_20260325_cell#44'});
FigResoSingleCells_vs = gca;
FigResoSingleCells_vs = ancestor(FigResoSingleCells_vs, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigResoSingleCells_vs, fullfile(outputDirFig2, 'VS_VSP_Data_examplefits\VS_VSP_Data_examplefits_allconditions-VISUALCELLS'));

% plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', {'M26004_20260318_cell#367'})
%% LLH VS vs VSP (scatter) 
plotFinalFigures_Sonali(EXP_all, 'LLHi-VSvsVSP-scatter')
figB = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figB, fullfile(outputDirFig2, 'LLH_improvement\LLH_improvement'));

%% RF centres


% plotFinalFigures_Sonali(EXP_all, 'LLHrel-RFposition-spatial-discrete')
plotFinalFigures_Sonali(EXP_all, 'LLHrel-RFposition-spatialOnly')
figC = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figC, fullfile(outputDirFig2, 'RFcentres_spatialcellsonly\RFcentres_spatialcellsonly'));


%% incldued in thesis 
close all; 
plotFinalFigures_Sonali(EXP_all, 'RFposition-Scatter-SpatialVsNonspatial')
figd = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figd, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.11\rfs_scatter_space_vis\', 'rfs_scatter_space_vis'));


%% included in thesis 
close all; 
plotFinalFigures_Sonali(EXP_all, 'RFposition-Latency-AllGoodcells')
fige = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(fige, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.11\rfs_split_visual_latency\', 'rfs_split_visual_latency'));

%% rf centre vs LLHrel 

plotFinalFigures_Sonali(EXP_all, 'RFposition-SplitViolin-SpatialVsNonspatial')
figC = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figC, fullfile(outputDirFig2, 'RFcentres_spatialAndNonSpatialcells\RFcentres_bothcelltypes'));
%% latency - compare spatial and non-spatial 
plotFinalFigures_Sonali(EXP_all, 'Latency-NonspatialVsSpatial-Corrected')
figD = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figD, fullfile(outputDirFig2, 'latency_spatial_non-spatial\latency_spatial_non-spatial'));
% Corrected Wilcoxon rank-sum test (non-spatial vs spatial latency): p = 2.32e-99


%%
plotFinalFigures_Sonali(EXP_all, 'Latency-RFm5VsRest-Corrected')
figD = gcf;
defaultAxesProperties(gca, 1);
saveFigureFormats(figD, fullfile(outputDirFig2, 'latency_RFm5_vs_rest\latency_RFm5_vs_rest'));
%%
plotFinalFigures_Sonali(EXP_all, 'RF-LLHrel-summary')
figSummary = gcf;

savepathtemp= 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.9\LLHSummary\'
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
    saveSubfolder = fullfile(savepathtemp, nameBase);
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

%%
RFpos = EXP_all.GLMs{1}.Perf.bestRFpos(:);
uniqueRF = unique(RFpos);
fprintf('Number of distinct RF centres: %d\n', numel(uniqueRF));
disp(uniqueRF);


%% spatial componenets
close all 
% plotFinalFigures_Sonali(EXP_all, 'Ordered-Kernels-Space')

plotFinalFigures_Sonali(EXP_all, 'ExampleTraces-Kernels-Space')
% -> shows Space kernel ordered by Space
FigSpaceKer = gca;
FigSpaceKer = ancestor(FigSpaceKer, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigSpaceKer, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.10\spaceKernel_withEgs\', 'space_with_measured_kernal_egs_mean3'));


%% Identify and pull up cells with negative RF positions
RFpos = EXP_all.GLMs{1}.Perf.bestRFpos(:);
negRF_cells = find(T.goodcells & RFpos < 0);

fprintf('%d cells with RFpos < 0\n', numel(negRF_cells));
negRF_cellNames = EXP_all.Spk.CellListString(negRF_cells);
for i = 1:numel(negRF_cellNames)
    fprintf('  %s\n', negRF_cellNames{i});
end
%% correaltion between fitted kernel and measured responses 
close all 
plotFinalFigures_Sonali(EXP_all, 'SpaceKernel-MeasuredVsEstimated-Corr')
FigcorrSpaceKernel_response = gca;
FigcorrSpaceKernel_response = ancestor(FigcorrSpaceKernel_response, 'figure');

defaultAxesProperties(gca, 1);
saveFigureFormats(FigcorrSpaceKernel_response, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_10_Section4\pearsonR_shapecorrelation\shapecorrelation_spacekernel_response'));


%% landmark and bg components % change to 10 exampels 
close all 
plotFinalFigures_Sonali(EXP_all, 'ScaledExamples-Kernels-Landmarks-BG')
FigLandBGKer = gca;
FigLandBGKer = ancestor(FigLandBGKer, 'figure');
% defaultAxesProperties(gca, 1);
%saveFigureFormats(FigLandBGKer, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_11_Section4\LandBGKernel_withEgs\', 'landmarks_BG_measured_kernal_egs'));


%%  in out rfs 
plotFinalFigures_Sonali(EXP_all, 'ScaledExamples-Kernels-Landmarks-BG-RFout')
T=plotFinalFigures_Sonali(EXP_all, 'ScaledExamples-Kernels-Landmarks-BG-RFin')


%% Identify and pull up cells with negative RF positions
RFpos = EXP_all.GLMs{1}.Perf.bestRFpos(:);
negRF_cells = find(T.goodcells & RFpos < 0);
fprintf('%d cells with RFpos < 0 (of %d goodcells)\n', numel(negRF_cells), sum(T.goodcells));
fprintf('RFpos range for these cells: %.1f to %.1f deg\n', ...
    min(RFpos(negRF_cells)), max(RFpos(negRF_cells)));

% Actual cell name IDs, not raw indices
negRF_cellNames = EXP_all.Spk.CellListString(negRF_cells);
fprintf('Cell names:\n');
for i = 1:numel(negRF_cellNames)
    fprintf('  %s\n', negRF_cellNames{i});
end

% Pull up full per-predictor kernel detail for the first several,
% using your existing single-cell plotting infrastructure
nExamples = min(5, numel(negRF_cellNames));
plotFinalFigures_Sonali(EXPToTest, 'GLM-kernels', 'M26003_20260325_cell#293')
%% 
close all 
plotFinalFigures_Sonali(EXP_all, 'ExampleTraces-KernelsOnly-Landmarks-BG')
FigLandBGKer_withS = gca;
FigLandBGKer_withS = ancestor(FigLandBGKer_withS, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigLandBGKer_withS, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.11\visualKernel_withEgs\', 'visual_kernal_egs_axesscaled_mean4'));



%%

close all 
plotFinalFigures_Sonali(EXP_all, 'BottomOfKernelHeatmap')
%%
close all 
plotFinalFigures_Sonali(EXP_all, 'ScaledExamples-Kernels-LLHr-Landmarks-BG')
FigLandBGKer_withS = gca;
FigLandBGKer_withS = ancestor(FigLandBGKer_withS, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigLandBGKer_withS, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.11\visualKernel_withEgs\', 'visual_kernal_egs_scaled_mean4'));




%%

% used the version used to scale seperately 
close all 
plotFinalFigures_Sonali(EXP_all, 'NotScaled-Kernels-LLHr-Landmarks-BG')
FigLandBGKer_withoutS = gca;
FigLandBGKer_withoutS = ancestor(FigLandBGKer_withoutS, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigLandBGKer_withoutS, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.11\visualKernel_withEgs\', 'visual_kernal_egs_mean4'));


close all; 
% this
plotFinalFigures_Sonali(EXP_all, 'ExampleTraces-Kernels-Landmarks-BG')
% plot examples to comapre estimated vs measured 
%plotFinalFigures_Sonali(EXP_all, 'LandmarkBGBalance-SpannedExamples')
FigLandBGKer_withmeasured = gca;
FigLandBGKer_withmeasured = ancestor(FigLandBGKer_withmeasured, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigLandBGKer_withmeasured, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_11_Section4\LandBGKernel_withEgs\', 'landmarks_BG_measured_kernal_egs_Scaled'));


% correlation between the landmark and bg kernels and the measured response
plotFinalFigures_Sonali(EXP_all, 'VisualKernel-MeasuredVsEstimated-Corr')
FigLandBGKer_corrmeanpresonse= gca;
FigLandBGKer_corrmeanpresonse = ancestor(FigLandBGKer_corrmeanpresonse, 'figure');
% defaultAxesProperties(gca, 1);
saveFigureFormats(FigLandBGKer_corrmeanpresonse, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_11_Section4\LandBGKernel_corr_with_measuredresponse\', 'LandBGKernel_corr_with_measuredresponse'));



% correlation with full model VSP 
plotFinalFigures_Sonali(EXP_all, 'FullModel-MeasuredVsEstimated-Corr')


%% are spatial components correalted with background and landmark components? 
% Spatial vs landmark kernel shape: median r = 0.043, MAD = 0.172, n = 3184
% Spatial vs background kernel shape: median r = -0.006, MAD = 0.171, n = 3184


plotFinalFigures_Sonali(EXP_all, 'SpatialVsVisual-ShapeCorr_withEgs')

%%
goodOmit_idx = find(T.goodOmitCells);
fprintf('Number of goodOmitCells: %d\n', numel(goodOmit_idx));
cellIDs_toShow = EXP_all.Spk.CellListString(goodOmit_idx(1:5));
plotFinalFigures_Sonali(EXP_all, 'Resp-singleCell', cellIDs_toShow);
%% plot the 6 background segments 
% this plots the 6 bg segments starting from the very beginning - therefore
% remember to reorder them due to the offset. 
plotVisualFeatures


%% plot the revealed omitted background texture: l2 l3 l4 

plotVisualFeatures_omittedbackground

%% check for omit cells
plotFinalFigures_Sonali(EXP_all, 'OmitCells-summary')

%% create csv with all cells 

%% find count per conditons across all sessions incldued in model fitting 



expt = [];
n = 0;
n = n + 1;
expt(n).animal = 'M25126';
expt(n).group = 'V1';
expt(n).BGtex = 'BG12';
expt(n).series = {20260311, 20260312, 20260313};
n = n + 1;
expt(n).animal = 'M25131';
expt(n).group = 'V1';
expt(n).BGtex = 'BG12';
expt(n).series = {20260318, 20260321, 20260322};
n = n + 1;
expt(n).animal = 'M26004';
expt(n).group = 'V1';
expt(n).BGtex = 'BG12';
expt(n).series = {20260318, 20260321, 20260322};
n = n + 1;
expt(n).animal = 'M26005';
expt(n).group = 'V1';
expt(n).BGtex = 'BG12';
expt(n).series = {20260318, 20260321, 20260322};


baseDir = 'C:\Users\sonali.sriranga\Desktop\V1\analyzed';

nSessions = sum(cellfun(@numel, {expt.series}));
condTrialCounts_perSession = nan(nSessions, 6);  % rows = sessions, cols = conditions
sessionsFound = 0;
sessionsMissing = {};

for a = 1:numel(expt)
    animal = expt(a).animal;
    for s = 1:numel(expt(a).series)
        seriesID = expt(a).series{s};
        sessionDir = fullfile(baseDir, animal, num2str(seriesID));

        d = dir(fullfile(sessionDir, 'SpkLin_EXP_*.mat'));

        if isempty(d)
            sessionsMissing{end+1} = sessionDir; %#ok<SAGROW>
            fprintf('MISSING (no matching file): %s\n', sessionDir);
            continue;
        end

        fpath = fullfile(d(1).folder, d(1).name);
        S = load(fpath, 'Nav');

        sessionsFound = sessionsFound + 1;
        for c = 1:6
            mask = S.Nav.condition == c;
            condTrialCounts_perSession(sessionsFound, c) = numel(unique(S.Nav.trialID(mask)));
        end
    end
end

condNames_report = {'base','swap23','swap34','omit2','omit3','omit4'};
fprintf('\n--- Per-condition trial counts (mean +/- SD across %d sessions) ---\n', sessionsFound);
for c = 1:6
    m = mean(condTrialCounts_perSession(1:sessionsFound, c), 'omitnan');
    sdv = std(condTrialCounts_perSession(1:sessionsFound, c), 'omitnan');
    tot = sum(condTrialCounts_perSession(1:sessionsFound, c), 'omitnan');
    fprintf('%s: %.1f +/- %.1f trials/session (pooled total: %d)\n', condNames_report{c}, m, sdv, tot);
end

% as these two were only presentred to three of the 4 mice
for c = [3 6]  % swap34, omit4
    vals = condTrialCounts_perSession(1:sessionsFound, c);
    vals_present = vals(vals > 0);   % exclude sessions/mice that didn't get this condition
    fprintf('Condition %d: %.1f +/- %.1f trials/session (n sessions = %d, pooled total = %d)\n', ...
        c, mean(vals_present), std(vals_present), numel(vals_present), sum(vals_present));
end


% sum manipulated-condition trial counts (columns 2-6) per session, then
% take mean/SD across sessions - this matches base's ± convention (SD
% across sessions), not SD across condition-means
manipTotal_perSession = sum(condTrialCounts_perSession(1:sessionsFound, 2:6), 2, 'omitnan');
fprintf('Manipulated conditions combined: %.1f +/- %.1f trials/session (N=%d sessions)\n', ...
    mean(manipTotal_perSession), std(manipTotal_perSession), sessionsFound);


%%

% count number of cells imaged in each session 
baseDir = 'C:\Users\sonali.sriranga\Desktop\V1\analyzed';

fprintf('%-10s %-12s %-10s\n', 'Animal', 'Series', 'nCells');
for a = 1:numel(expt)
    animal = expt(a).animal;
    for s = 1:numel(expt(a).series)
        seriesID = expt(a).series{s};
        sessionDir = fullfile(baseDir, animal, num2str(seriesID));
        d = dir(fullfile(sessionDir, 'SpkLin_EXP_*.mat'));
        if isempty(d)
            fprintf('%-10s %-12d %-10s\n', animal, seriesID, 'MISSING');
            continue;
        end
        fpath = fullfile(d(1).folder, d(1).name);
        S = load(fpath, 'Spk');
        nCells = numel(S.Spk.CellListString);
        fprintf('%-10s %-12d %-10d\n', animal, seriesID, nCells);
    end
end


%% landmark and bg kernel metrics

plotFinalFigures_Sonali(EXP_all, 'VisualKernel-MeasuredVsEstimated-Corr');
plotFinalFigures_Sonali(EXP_all,'VisualKernel-Amplitude-Compare');
plotFinalFigures_Sonali(EXP_all,'VisualKernel-Amplitude-Compare-Violin');


%%

plotFinalFigures_Sonali(EXP_all, 'RFcompare-LandmarkDist')
plotFinalFigures_Sonali(EXP_all, 'NegRF-LandmarkDist-Visualise')

plotFinalFigures_Sonali(EXP_all, 'RFm5-InOutSplit')
plotFinalFigures_Sonali(EXP_all, 'RFm5-InOutSplit-KernelAndData') % looks ok ish
%% load 
resp025 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_Response_M26003_LandManipCorridor_20260325_CombinedRuns.mat");
sfi025 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260325\M26003_20260325_sessionFileInfo.mat");


%
resp024 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_Response_M26003_LandManipCorridor_20260324_CombinedRuns.mat");
sfi024 = load("Z:\ibn-vision\DATA\SUBJECTS\M26003\Analysis\20260324\M26003_20260324_sessionFileInfo.mat");

% pushed name into the worspace
plotFinalFigures_Sonali(EXP_all, 'NegRF-GetCandidateNames')

% Loop through candidate cells, loading each session's response only ONCE
signalToUse = 'spks';
pauseSeconds = 2;

% Parse all names first
nCells = numel(cellNames);
mouseList = cell(nCells,1);
sessionList = cell(nCells,1);
neuronIdxList = nan(nCells,1);
for i = 1:nCells
    parts = strsplit(cellNames{i}, '_');
    mouseList{i} = parts{1};
    sessionList{i} = parts{2};
    neuronIdxList(i) = str2double(regexprep(parts{3}, 'cell#', ''));
end

% Group by unique mouse+session combo
sessionKey = strcat(mouseList, '_', sessionList);
[uniqueKeys, ~, keyIdx] = unique(sessionKey, 'stable');

for s = 1:numel(uniqueKeys)
    cellsInThisSession = find(keyIdx == s);
    mousenumber = mouseList{cellsInThisSession(1)};
    sessionName = sessionList{cellsInThisSession(1)};

    fprintf('\n--- Loading session %s_%s (%d cells) ---\n', ...
        mousenumber, sessionName, numel(cellsInThisSession));

    try
        infoPath = findSessionFileInfoFilePath(mousenumber, sessionName);
        if ~isfile(infoPath)
            fprintf('  Info missing — skipping session\n');
            continue
        end
        loadedInfo = load(infoPath, 'sessionFileInfo');
        sessionFileInfo = loadedInfo.sessionFileInfo;
        stimNames = {sessionFileInfo.stimFiles.name};

        targetIdx = find(contains(stimNames, "Corridor") & contains(stimNames, "CombinedRuns"), 1);
        if isempty(targetIdx)
            allCorridorIdx = find(contains(stimNames, "Corridor"));
            if isscalar(allCorridorIdx)
                targetIdx = allCorridorIdx;
            end
        end
        if isempty(targetIdx)
            fprintf('  No corridor stimulus found — skipping session\n');
            continue
        end

        response = load(sessionFileInfo.stimFiles(targetIdx).Response); % loaded ONCE per session
    catch ME
        fprintf('  Skipped session %s_%s — %s\n', mousenumber, sessionName, ME.message);
        continue
    end

    % Now loop through just this session's candidate cells, no reloading
    for c = 1:numel(cellsInThisSession)
        i = cellsInThisSession(c);
        neuronIdx = neuronIdxList(i);
        fprintf('  Cell %d (meanDist=%.1f cm)\n', neuronIdx, cellDistances(i));

        try
            figHandle = plotSingleROI_ConditionStack(sessionFileInfo, response, neuronIdx, signalToUse, false);
            pause(pauseSeconds);
            close(figHandle);
        catch ME
            fprintf('    Skipped cell %d — %s\n', neuronIdx, ME.message);
        end
    end
end


%% ran this to get an example of cells that have early or late peaks
plotFinalFigures_Sonali(EXP_all, 'ScaledExamples-Kernels-Landmarks-BG')


%% test if a cell exists 




%% Check whether M26003_20260325 cell#83 exists in the GLM/kernel population, and cross-check its properties
targetName = 'M26003_20260325_cell#83';

matchIdx = find(strcmp(EXP_all.Spk.CellListString, targetName));

if isempty(matchIdx)
    fprintf('%s was NOT found in EXP_all.Spk.CellListString.\n', targetName);
    fprintf('This could mean it failed goodcells/spatial-fit criteria, or was excluded upstream.\n');
else
    fprintf('%s found at index %d\n', targetName, matchIdx);

    isGood = T.goodcells(matchIdx);
    isSpatial = T.signicells(matchIdx) & goodLLHcells(matchIdx);
    RFpos = EXP_all.GLMs{1}.Perf.bestRFpos(matchIdx);

    fprintf('  goodcells: %d\n', isGood);
    fprintf('  spatially modulated: %d\n', isSpatial);
    fprintf('  RF position: %.1f deg\n', RFpos);

    % Cross-check landmark kernel peak position against the raw-data peak (112.6cm)
    mapLandmark_idx = (1:ncond) + 3*ncond;
    landmarkTrace = squeeze(EXP_all.Maps{mapLandmark_idx(1)}.Tuning.meanrespModel);
    xpos = linspace(0, 200, size(landmarkTrace, 2));

    kernelTrace = landmarkTrace(matchIdx, :);
    [~, peakBin] = max(kernelTrace);
    kernelPeakPos = xpos(peakBin);

    fprintf('  Kernel-derived peak position: %.1f cm (raw-data peak was 112.6 cm)\n', kernelPeakPos);
    fprintf('  Match: %s\n', string(abs(kernelPeakPos - 112.6) < 5));
end


%% 


%% Settings you may want to change
breakpt        = 80;    % x value where compression starts
xmax           = 125;   % max x value to display
compressFactor = 0.25;  % smaller = more squeeze above breakpt
 
x_data = linspace(0, xmax, 500);
y_data = sin(x_data/10) + 0.05*randn(size(x_data));
 
%% Transform the data (inline, no separate function)
xt_data = x_data;
above = x_data > breakpt;
xt_data(above) = breakpt + (x_data(above) - breakpt) * compressFactor;
 

figweirdticks = figure; 
hold on;
plot(xt_data, y_data, 'k-', 'LineWidth', 1.2);
ax = gca;
 

% Normal spacing below the breakpoint. Denser real values above it,
% so they end up visually packed together after transforming.
ticks_below = 0:20:breakpt;
ticks_above = (breakpt + 5):5:xmax;
allticks    = [ticks_below, ticks_above];
 
% Transform the tick positions the same way (inline again). Keep the
% labels showing the true, untransformed values.
allticks_t = allticks;
above_t = allticks > breakpt;
allticks_t(above_t) = breakpt + (allticks(above_t) - breakpt) * compressFactor;
 
xmax_t = xmax;
if xmax > breakpt
    xmax_t = breakpt + (xmax - breakpt) * compressFactor;
end
 
set(ax, 'XTick', allticks_t, 'XTickLabel', allticks);
xlim([0, xmax_t]);
 

minorticks = setdiff(0:5:breakpt, allticks);   % skip values already major
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = minorticks;         % positions below breakpt need no transform
 
xlabel('x');
ylabel('y');
title('Compressed axis: 0-80 normal, 80-125 compressed');
box off; set(ax, 'TickDir', 'out');


figweirdticks= gca;
figweirdticks = ancestor(figweirdticks, 'figure');
defaultAxesProperties(gca, 1);
saveFigureFormats(figweirdticks, fullfile('Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section4_Fig3.10\axesforschematic\', 'axes_for_fig'));

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