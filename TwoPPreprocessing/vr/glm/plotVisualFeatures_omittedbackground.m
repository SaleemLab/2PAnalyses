%% Background texture revealed by landmark omission
% Standalone script - loads BG2.jpg directly, shows the full corridor
% strip with omitted-landmark positions (L2, L3, L4 by sequence) boxed,
% then a wider context window around each with the exposed segment
% highlighted.

%% Parameters
basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';
saveSubfolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_7_Section3';
nameBase = 'Omission_BGreveal';

corridorL = 200;              % corridor length (cm)
landmark_width_pct = 4;       % landmark width, % of corridor
window_pct = 20;              % wider context window shown per omission, % of corridor

% Landmark sequence along corridor: L1@20%, L2@40%, L1@60%, L2@80%
% (sequence positions 1,2,3,4). Omitted landmarks are sequence positions 2,3,4.
seqPositions_pct = [40, 60, 80];
seqLabels = {'L2', 'L3', 'L4'};

col_omit = [12 95 196; 124 11 161; 11 124 161] / 255;   % Omit L2 / L3 / L4 colours

%% Load texture
img = imread(fullfile(basePath, 'BG2.jpg'));
if isa(img, 'uint8')
    img = double(img) / 255;
elseif isa(img, 'logical')
    img = double(img);
end
if size(img, 3) > 1
    img = rgb2gray(img);
end
[h, w] = size(img);
px_per_cm = w / corridorL;

hw = landmark_width_pct / 2;
winHalf = window_pct / 2;

pct2px = @(pct) round(pct/100 * corridorL * px_per_cm);

%% Plot
figOm = figure('Name', 'Omission - BG reveal', 'Position', [100 100 1100 750]);

% --- top: full corridor strip with omission positions boxed ---
axFull = subplot(2, 3, [1 2 3]);
imagesc([0 corridorL], [0 1], img, [0 1]);
colormap(axFull, gray); axis tight;
set(axFull, 'YTick', []);
xlabel('Corridor position (cm)');
title('Landmark positions omitted on omission trials');
hold on;
for k = 1:numel(seqPositions_pct)
    c = seqPositions_pct(k);
    s_pct = c - hw; e_pct = c + hw;
    rectangle('Position', [s_pct/100*corridorL, 0, (e_pct-s_pct)/100*corridorL, 1], ...
              'EdgeColor', col_omit(k,:), 'LineWidth', 2.5);
    text(c, 1.12, sprintf('Omit %s', seqLabels{k}), 'Color', col_omit(k,:), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end
hold off;

% --- bottom: wider context window per omission, exposed segment highlighted ---
for k = 1:numel(seqPositions_pct)
    c = seqPositions_pct(k);
    win_s = c - winHalf; win_e = c + winHalf;
    s_px = pct2px(win_s) + 1;
    e_px = pct2px(win_e);
    cropWin = img(:, s_px:e_px);

    ax = subplot(2, 3, 3 + k);
    imagesc([win_s win_e], [0 1], cropWin, [0 1]);
    colormap(ax, gray); axis tight;
    set(ax, 'YTick', []);
    xlabel('Corridor position (%)');
    title(sprintf('Omit %s\nexposed: %.0f-%.0f%%', seqLabels{k}, c-hw, c+hw), ...
          'Color', col_omit(k,:), 'FontWeight', 'bold');
    hold on;
    rectangle('Position', [c-hw, 0, landmark_width_pct, 1], ...
              'EdgeColor', col_omit(k,:), 'LineWidth', 3);
    hold off;
end

sgtitle('Background texture revealed by landmark omission', 'FontSize', 13);

%% Save
saveFigureFormats(figOm, fullfile(saveSubfolder, nameBase));