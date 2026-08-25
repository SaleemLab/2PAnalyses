%% fov_retinotopy_overlay.m
% Overlay per-cell preferred RF azimuth onto the actual imaging FOV
% (mean image) for one example session, ONE PLANE.
% Peak location chosen by RF amplitude; significance = p<0.05 (uncorrected)
% at that single peak location.

clear; close all;

matFile = "Z:\ibn-vision\DATA\SUBJECTS\M25131\Analysis\20260318\M25131_20260318_sessionROIData.mat";
opsFile = "Z:\ibn-vision\DATA\SUBJECTS\M25131\Analysis\20260318\M25131_20260318_2pData_M25131_SparseNoiseTexture_20260318_00001.mat";

selectedPlane = 3;

bgField    = 'max_proj';
mapTypeIdx = 4;

azRange = [-80, 20];
elRange = [-40, 60];

opsData = load(opsFile);
opsData = opsData.twoPData(selectedPlane);
ops = opsData.ops;
bgImg = double(ops.(bgField));
[Ly, Lx] = size(bgImg);

S = load(matFile);
roiInfo = S.roiInfo;
sparseNoiseRF = S.sparseNoiseRF;

iscell        = logical(roiInfo.iscell);
planeIdentity = double(roiInfo.roiPlaneIdentity);
xpixAll       = roiInfo.xpix;
ypixAll       = roiInfo.ypix;
initMap       = sparseNoiseRF.initMap;
initPMap      = sparseNoiseRF.initPMap;
gridSize      = sparseNoiseRF.gridSize;

nCells = numel(xpixAll);
planeMask = (planeIdentity == selectedPlane);

az = linspace(azRange(1), azRange(2), gridSize(2));
el = linspace(elRange(1), elRange(2), gridSize(1));

prefAzDeg = nan(nCells,1);
pValsAtPeak = nan(nCells,1);

for i = 1:nCells
    if ~planeMask(i) || isempty(initMap{i})
        continue
    end
    rfMap = initMap{i}(:, :, end, mapTypeIdx);
    pMap  = initPMap{i}(:, :, end, mapTypeIdx);

    [~, linIdx] = max(rfMap(:));
    [~, peakAzIdx] = ind2sub(size(rfMap), linIdx);

    prefAzDeg(i) = az(peakAzIdx);
    pValsAtPeak(i) = pMap(linIdx);
end

validCells = planeMask & ~isnan(prefAzDeg);

isSignificant = false(nCells,1);
isSignificant(validCells) = pValsAtPeak(validCells) < 0.05;
cellsToPlot = validCells & isSignificant;

fprintf('%d / %d ROIs (plane %d) significant at p<0.05 (uncorrected)\n', ...
    sum(cellsToPlot), sum(validCells), selectedPlane);

cmap = jet(256);
overlayRGB  = zeros(Ly, Lx, 3);
overlayMask = false(Ly, Lx);

for i = 1:nCells
    if ~cellsToPlot(i)
        continue
    end
    xp = double(xpixAll{i}(:)) + 1;
    yp = double(ypixAll{i}(:)) + 1;

    validPix = (xp >= 1 & xp <= Lx) & (yp >= 1 & yp <= Ly);
    xp = xp(validPix);
    yp = yp(validPix);
    if isempty(xp)
        continue
    end

    colorIdx = round((prefAzDeg(i) - azRange(1)) / diff(azRange) * 255) + 1;
    colorIdx = min(max(colorIdx, 1), 256);
    thisColor = cmap(colorIdx, :);

    linInd = sub2ind([Ly, Lx], yp, xp);
    for c = 1:3
        chan = overlayRGB(:,:,c);
        chan(linInd) = thisColor(c);
        overlayRGB(:,:,c) = chan;
    end
    overlayMask(linInd) = true;
end

figure('Units','centimeters','Position',[2 2 9 8]);

bgClipped = mat2gray(bgImg, [prctile(bgImg(:),1), prctile(bgImg(:),99)]);
bgRGB = repmat(bgClipped, [1 1 3]);

imshow(bgRGB); hold on;

overlayAlpha = 0.85 * overlayMask;
hOverlay = imshow(overlayRGB);
set(hOverlay, 'AlphaData', overlayAlpha);

axis image off;
title(sprintf('Preferred RF azimuth (plane %d)', selectedPlane), 'FontSize', 9);

cbAxes = axes('Position', [0.83 0.15 0.04 0.7]);
nSteps = 256;
gradientVals = linspace(azRange(1), azRange(2), nSteps)';
gradientImg = reshape(cmap, [nSteps, 1, 3]);

imagesc(cbAxes, [0 1], gradientVals, gradientImg);
set(cbAxes, 'YDir', 'normal');
set(cbAxes, 'XTick', []);
cbAxes.YAxisLocation = 'right';
ylabel(cbAxes, 'Azimuth (\circ)', 'FontSize', 8);
box(cbAxes, 'on');

set(gcf, 'Color', 'w');

%%
% %% population_rf_visual_field_map_multi_session.m
% % Pool significant cells' RF peak locations across multiple sessions/FOVs
% % and plot in visual-field space (azimuth x elevation).
% 
% clear; close all;
% 
% sessionFiles = {
%     "Z:\ibn-vision\DATA\SUBJECTS\M25131\Analysis\20260318\M25131_20260318_sessionROIData.mat"
%     "Z:\ibn-vision\DATA\SUBJECTS\M25131\Analysis\20260321\M25131_20260321_sessionROIData.mat"
%     "Z:\ibn-vision\DATA\SUBJECTS\M25131\Analysis\20260322\M25131_20260322_sessionROIData.mat"
%     "Z:\ibn-vision\DATA\SUBJECTS\M26005\Analysis\20260318\M26005_20260318_sessionROIData.mat"
%     "Z:\ibn-vision\DATA\SUBJECTS\M26005\Analysis\20260321\M26005_20260321_sessionROIData.mat"
%     "Z:\ibn-vision\DATA\SUBJECTS\M26005\Analysis\20260322\M26005_20260322_sessionROIData.mat"
%     % add more session paths here...
% };
% 
% mapTypeIdx = 4;
% azRange = [-90, 90];
% elRange = [-40, 60];
% pThresh = 0.05;
% 
% allAz = [];
% allEl = [];
% allMaps = [];
% 
% for s = 1:numel(sessionFiles)
%     S = load(sessionFiles{s});
%     roiInfo = S.roiInfo;
%     sparseNoiseRF = S.sparseNoiseRF;
% 
%     initMap  = sparseNoiseRF.initMap;
%     initPMap = sparseNoiseRF.initPMap;
%     gridSize = sparseNoiseRF.gridSize;
% 
%     az = linspace(azRange(1), azRange(2), gridSize(2));
%     el = linspace(elRange(1), elRange(2), gridSize(1));
% 
%     nCells = numel(initMap);
% 
%     for i = 1:nCells
%         if isempty(initMap{i})
%             continue
%         end
%         rfMap = initMap{i}(:, :, end, mapTypeIdx);
%         pMap  = initPMap{i}(:, :, end, mapTypeIdx);
% 
%         [~, linIdx] = max(rfMap(:));
%         [peakElIdx, peakAzIdx] = ind2sub(size(rfMap), linIdx);
% 
%         if pMap(linIdx) < pThresh
%             allAz(end+1,1) = az(peakAzIdx);
%             allEl(end+1,1) = el(peakElIdx);
%             allMaps = cat(3, allMaps, rfMap);
%         end
%     end
%     fprintf('Session %d (%s): pooled cells so far = %d\n', s, sessionFiles{s}, numel(allAz));
% end
% 
% fprintf('Total significant cells pooled across %d sessions: %d\n', numel(sessionFiles), numel(allAz));
% 
% %% --- Aggregate heatmap across ALL pooled cells ---
% meanRFMap = mean(allMaps, 3, 'omitnan');
% meanRFMapNorm = (meanRFMap - min(meanRFMap(:))) / (max(meanRFMap(:)) - min(meanRFMap(:)));
% 
% %% --- Plot ---
% figure('Units','centimeters','Position',[2 2 9 7]);
% 
% imagesc(az, el, meanRFMapNorm);
% set(gca, 'YDir', 'normal');
% colormap(parula);
% cb = colorbar;
% cb.Label.String = 'response';
% caxis([0 1]);
% hold on;
% 
% xlabel('azimuth (deg)');
% ylabel('elevation (deg)');
% xline(0, 'k:', 'Alpha', 0.5);
% yline(0, 'k:', 'Alpha', 0.5);
% axis equal tight;
% 
% plot(allAz, allEl, 'ro', 'MarkerSize', 4, 'LineWidth', 0.75);
% 
% title(sprintf('Population RF map (n = %d cells, %d sessions)', numel(allAz), numel(sessionFiles)));