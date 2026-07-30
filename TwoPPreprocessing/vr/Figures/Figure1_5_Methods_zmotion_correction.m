%% make_zregistration_plots.m
% MATLAB translation of the Python make_plots() function - generates the
% three z-registration diagnostic figures:
%   1. Reference image per plane, best reference highlighted
%   2. Final combined reference vs mean image after z-registration
%   3. Z-positions of reference images over time + correlation trace
%      with bad frames shaded
% From sylvia's pipeline converted to matlab

%% Run this file to align and plot the average piezo trace as a function of time within a frame 

 align_piezo_to_frames
%%
fallMatPath = "Z:\ibn-vision\DATA\SUBJECTS\M25132\Processed\20260326\suite2p\plane_z\Fall.mat";
corrThresh = 3;

%%
S = load(fallMatPath, 'ops');
ops = S.ops;

%%  reference images per plane
allPlanes = 0:(double(ops.nplanes) - 1);
ignorePlanes = double(ops.ignore_flyback_singleplanes);
planeIds = setdiff(allPlanes, ignorePlanes);

refImgs = ops.refImg_singleplanes;

if iscell(refImgs)
    nRefs = numel(refImgs);
    getRefImg = @(i) refImgs{i};
else
    % Stacked numeric array - assume first dimension indexes plane
    nRefs = size(refImgs, 1);
    getRefImg = @(i) squeeze(refImgs(i, :, :));
end

nRows = floor(sqrt(nRefs));
nCols = ceil(nRefs / nRows);

figb = figure('Name', 'Reference Images per Plane', 'Color', 'w', 'Position', [50 50 1200 800]);
for i = 1:nRefs
    subplot(nRows, nCols, i);
    imagesc(getRefImg(i));
    colormap(gray);
    axis image off;
    if planeIds(i) == double(ops.reference_plane)
        title(sprintf('Plane %d (best reference)', planeIds(i)), 'Color', 'r');
    else
        title(sprintf('Plane %d', planeIds(i)));
    end
end
sgtitle('Reference Images per Plane');

defaultAxesProperties(gca, 1);

baseFileName = "reference_image_m25132_20260326";
outputDir = "Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\Methods\Fig2.5_zmotioncCorrection\ref_image_per_plane";

if ~exist(outputDir, 'dir'), mkdir(outputDir); end


fullSavePath = char(fullfile(outputDir, baseFileName));
saveFigureFormats(figb, fullSavePath);


%%  Plot 2: final reference vs mean image after z-registration ----


figc = figure('Name', 'Mean Image after Z-Registration', 'Color', 'w', 'Position', [50 50 1400 400]);


if double(ops.align_by_chan) == 1
    meanImg = ops.meanImg;
    nColsFig2 = 3;
else
    meanImg = ops.meanImg_chan2;
    nColsFig2 = 4;
end

subplot(1, nColsFig2, 1);
imagesc(ops.refImg); colormap(gray); axis image off;
title('Final Reference Image (weighted neighbors)');

subplot(1, nColsFig2, 2);
lo = prctile(meanImg(:), 1);
hi = prctile(meanImg(:), 99);
imagesc(meanImg, [lo, hi]); colormap(gray); axis image off;
title('Mean Image after Z-Registration');

subplot(1, nColsFig2, 3);
refLo = prctile(ops.refImg(:), 1); refHi = prctile(ops.refImg(:), 99);
meanLo = prctile(meanImg(:), 1); meanHi = prctile(meanImg(:), 99);

refNorm = (ops.refImg - refLo) / (refHi - refLo);
meanNorm = (meanImg - meanLo) / (meanHi - meanLo);
refNorm = min(max(refNorm, 0), 1);     % clip to [0,1] after percentile normalization
meanNorm = min(max(meanNorm, 0), 1);

rgbImg = zeros([size(meanImg), 3]);
rgbImg(:,:,1) = refNorm;   % Red
rgbImg(:,:,2) = meanNorm;  % Green
imagesc(rgbImg); axis image off;
title('Mean Image (green) vs Reference (red)');

if double(ops.align_by_chan) == 2
    subplot(1, nColsFig2, 4);
    meanImgC1 = ops.meanImg;
    lo1 = prctile(meanImgC1(:), 1);
    hi1 = prctile(meanImgC1(:), 99);
    imagesc(meanImgC1, [lo1, hi1]); colormap(gray); axis image off;
    title('Mean Image Channel 1');
end

defaultAxesProperties(gca, 1);

baseFileName = "mean_image_zregistration_m25132_20260326";
outputDir = "Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\Methods\Fig2.5_zmotioncCorrection\mean_image_zregistration";

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = char(fullfile(outputDir, baseFileName));
saveFigureFormats(figc, fullSavePath);

%%  z-positions over time + correlation with best reference ----
bestReference = find(planeIds == double(ops.reference_plane), 1) - 1;   % 0-indexed to match Python logic

% corrs_time_refs_planes: [n_frames x n_refs x n_planes]
corrsTimeRefsPlanes = ops.corrs_time_refs_planes;
[nFrames, ~, nPlanesUsed] = size(corrsTimeRefsPlanes);

[~, planesTimeRefs] = max(corrsTimeRefsPlanes, [], 3, 'omitnan');
planesTimeRefs = planesTimeRefs - 1;   % convert to 0-indexed to match planeIds

% frames_per_folder may be a scalar (single experiment) or a vector
% (multiple experiments concatenated) - handle both.
framesPerFolder = double(ops.frames_per_folder);
if isscalar(framesPerFolder)
    expStarts = [];   % only one experiment - no internal boundaries to mark
else
    expStarts = cumsum(framesPerFolder);
    expStarts = expStarts(1:end-1);
end

figd = figure('Name', 'Z-Positions', 'Color', 'w', 'Position', [50 50 1200 700]);


planeColors = [
    0.20 0.20 0.60;   
    0.20 0.40 0.80;   
    0.20 0.70 0.90;   
    0.30 0.80 0.30;   
    0.80 0.80 0.20;  
    0.95 0.60 0.10;   
    0.90 0.20 0.20;  
    0.60 0.10 0.10;   
];

subplot(2,1,1);
hold on;
for p = 1:nPlanesUsed
    thisPlaneNum = planeIds(p);
    if (p - 1) == bestReference
        offset = 0;
        color = 'r';
    else
        offset = ((p - 1) - bestReference) * numel(planeIds);
        color = planeColors(thisPlaneNum + 1, :);
    end
    yVals = planeIds(planesTimeRefs(:, p) + 1) + offset;
    scatter(1:nFrames, yVals, 4, color, 'filled', 'DisplayName', sprintf('Plane %d', thisPlaneNum));
end
for es = expStarts'
    xline(es, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
title('Z-Positions for all reference images');
ylabel('Best imaging plane (+ offset for non-ref planes)');
legend('show', 'Location', 'eastoutside');
xlim([1, nFrames]);
box off;
defaultAxesProperties(gca, 1);


subplot(2,1,2);
hold on;

corrsTimePlanes = squeeze(corrsTimeRefsPlanes(:, bestReference + 1, :));
corrsTime = max(corrsTimePlanes, [], 2, 'omitnan');
threshold = mean(corrsTime, 'omitnan') - corrThresh * std(corrsTime, 'omitnan');

yPad = 0.1 * range(corrsTime);
yMin = min(corrsTime) - yPad;
yMax = max(corrsTime) + yPad;

badframes = corrsTime <= threshold;

plot(corrsTime, 'k-', 'LineWidth', 0.5, 'DisplayName', 'Correlation');
yline(threshold, 'k-', 'LineWidth', 1, 'DisplayName', 'Corr. thresh.');
for es = expStarts'
    xline(es, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');   % FIX: 'gray' is not a valid xline color spec
end

inBad = false;
startIdx = NaN;
for i = 1:numel(badframes)
    if badframes(i) && ~inBad
        startIdx = i;
        inBad = true;
    elseif ~badframes(i) && inBad
        patch([startIdx-0.5, i-0.5, i-0.5, startIdx-0.5], ...
            [yMin, yMin, yMax, yMax], ...
            [0.15 0.15 0.15], 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        inBad = false;
    end
end
if inBad
    patch([startIdx-0.5, numel(badframes)-0.5, numel(badframes)-0.5, startIdx-0.5], ...
          [yMin, yMin, yMax, yMax], ...
          [0.15 0.15 0.15], 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

title('Correlation: best reference with frames');
xlabel('Frame');
ylabel('Correlation coefficient');
legend('show', 'Location', 'eastoutside');
xlim([1, nFrames]);
ylim([yMin, yMax]);
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.3f');
box off;

defaultAxesProperties(gca, 1);

% Align both subplots' x-axes despite differing legend widths ----
% Must run AFTER both subplots and their legends are fully drawn, since
% creating a legend is what resizes each axes independently in the first
% place - doing this earlier would just get undone by a later legend call.
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');

newLeft  = max(pos1(1), pos2(1));
newWidth = min(pos1(3), pos2(3));

pos1(1) = newLeft; pos1(3) = newWidth;
pos2(1) = newLeft; pos2(3) = newWidth;

set(ax1, 'Position', pos1);
set(ax2, 'Position', pos2);


baseFileName = "z_positions_m25132_20260326";
outputDir = "Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\Methods\Fig2.5_zmotioncCorrection\z_positions";

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fullSavePath = char(fullfile(outputDir, baseFileName));
saveFigureFormats(figd, fullSavePath);