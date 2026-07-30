function plot_texture_layout(pos0_cm, basePath)
% PLOT_TEXTURE_LAYOUT
% Standalone function producing Figure 1D with Panel B in Retinotopic Coordinates.

if nargin < 1 || isempty(pos0_cm)
    pos0_cm = 60; % Default viewer position in cm
end
if nargin < 2
    basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';
end

%% --------------------------- PARAMETERS -------------------------------
corridorL = 200;       % Corridor length (cm)
corridorH = 12;        % Corridor height (cm)
texwidth  = 0.04;      % Landmark width fraction
halfW_cm  = corridorH / 2; % Half corridor width = 6 cm

az_min = 0.5;
az_max = 85;
screen_limit = 80;
az_bins = 500;
az_vec  = linspace(az_min, az_max, az_bins);

%% --------------------------- LOAD TEXTURES ----------------------------
textures = struct('matrix', []);

% Load Background 2
textures(3).matrix = imread_and_prep(fullfile(basePath, 'BG2.jpg'));

% Load Specific Landmarks: grey (1), grating_vertical (6), plaid (8)
landmarkFiles = {'grey.jpg', 'grating_vertical.jpg', 'plaid.jpg'};
landmarkIdx   = [1, 6, 8]; 
for i = 1:length(landmarkFiles)
    img = imread_and_prep(fullfile(basePath, landmarkFiles{i}));
    textures(landmarkIdx(i)).matrix = img;
end

%% ---------------- BUILD COMPOSITE CORRIDOR (CATALOGUE LOGIC) ----------------
targetBG      = textures(3).matrix;
finalBGlength = size(targetBG, 2);
finalBGheight = size(targetBG, 1);
final_width   = round(finalBGlength * texwidth);

grating_src = textures(6).matrix;
plaid_src   = textures(8).matrix;
end_src     = textures(1).matrix;

% Interpolate landmark sizes to match background height and width
[u_p, v_p]   = meshgrid(1:size(plaid_src,2), 1:size(plaid_src,1));
[uq_p, vq_p] = meshgrid(linspace(1, size(plaid_src,2), final_width), linspace(1, size(plaid_src,1), finalBGheight));
grating_res = interp2(u_p, v_p, grating_src, uq_p, vq_p);
plaid_res   = interp2(u_p, v_p, plaid_src, uq_p, vq_p);

% Interpolate End Wall
[u_e, v_e]   = meshgrid(1:size(end_src,2), 1:size(end_src,1));
[uq_e, vq_e] = meshgrid(linspace(1, size(end_src,2), final_width), linspace(1, size(end_src,1), finalBGheight));
end_res     = interp2(u_e, v_e, end_src, uq_e, vq_e);

% Assemble Composite Texture (Horizontal: 0 -> 200 cm)
tex_composite = targetBG; 
centers_px = [0.20, 0.40, 0.60, 0.80] * finalBGlength;

% Overlay Landmarks (G - P - G - P)
for idx = 1:4
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    if mod(idx,2) == 1
        tex_composite(:, start_col:start_col+final_width-1) = grating_res;
    else
        tex_composite(:, start_col:start_col+final_width-1) = plaid_res;
    end
end

% Overlay End Wall at the end of corridor (200 cm)
tex_composite(:, end-final_width+1:end) = end_res;

%% ---------------- RETINOTOPIC MAPPING MATRIX --------------------------
% Generate Meshgrid: Y = Mouse position along corridor (cm), X = Visual Azimuth (deg)
pos_vec = linspace(0, corridorL, finalBGlength);
[AZ_grid, POS_grid] = meshgrid(az_vec, pos_vec);

% Calculate wall hit column indices (X_wall)
X_wall_hit = POS_grid + halfW_cm * (1 ./ tand(AZ_grid));
col_coords = (X_wall_hit / corridorL) * (finalBGlength - 1) + 1;

% For full 2D retinotopic projection, sample across the vertical wall rows
row_indices = round(linspace(1, finalBGheight, finalBGlength))';
row_coords  = repmat(row_indices, 1, az_bins);

% Sample composite texture in 2D space using interp2
retinotopic_mat = interp2(1:finalBGlength, 1:finalBGheight, tex_composite, col_coords, row_coords, 'linear', 0);

%% ------------------------------ PLOTTING ------------------------------
figure('Color', 'w', 'Position', [100 100 850 650]);

% --- Left Subplot (Panel A): Corridor Profile (0-200 cm Vertical) ---
ax_left = axes('Position', [0.08 0.10 0.18 0.82]);

% Transpose composite image so position (0-200 cm) runs vertically along Y-axis
imagesc(ax_left, [0 1], [0 corridorL], tex_composite'); 
colormap(ax_left, 'gray');
set(ax_left, 'YDir', 'normal');
ylabel(ax_left, 'Position along corridor (cm)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel(ax_left, ''); set(ax_left, 'XTick', []);
ylim(ax_left, [0 corridorL]); 
box(ax_left, 'off'); set(ax_left, 'TickDir', 'out');
title(ax_left, 'Corridor Profile', 'FontSize', 11, 'FontWeight', 'bold');
hold(ax_left, 'on');

% Viewer position marker
yline(ax_left, pos0_cm, '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
plot(ax_left, 0.5, pos0_cm, 'v', 'MarkerFaceColor', [0 0.45 0.74], ...
    'MarkerEdgeColor', 'none', 'MarkerSize', 8);

% --- Right Subplot (Panel B): Visual Azimuth Retinotopic Projection ---
ax_main = axes('Position', [0.38 0.10 0.55 0.82]);

imagesc(ax_main, az_vec, pos_vec, retinotopic_mat, [0 1]);
colormap(ax_main, 'gray');
set(ax_main, 'YDir', 'normal');
xlabel(ax_main, 'Visual azimuth (°)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ax_main, ''); set(ax_main, 'YTick', []);
xlim(ax_main, [0 az_max]); ylim(ax_main, [0 corridorL]);
box(ax_main, 'off'); set(ax_main, 'TickDir', 'out');
hold(ax_main, 'on');

% Off-screen section shading (80° - 85°)
patch(ax_main, [screen_limit az_max az_max screen_limit], [0 0 corridorL corridorL], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xline(ax_main, screen_limit, ':k', 'LineWidth', 1.2);

% RF Window Box Overlay
rf_color = [0.85 0.10 0.55];
rectangle(ax_main, 'Position', [0.5, 2, screen_limit - 0.5, corridorL - 4], ...
    'EdgeColor', rf_color, 'LineWidth', 1.8, 'LineStyle', '--');
text(ax_main, 40, 12, 'RF_{win}', 'Color', rf_color, ...
    'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% --- Top Annotation Header Bar ---
axes('Position', [0.38 0.93 0.55 0.03]);
hold on; axis off; xlim([0 az_max]); ylim([0 1]);
patch([0 15 15 0], [0 0 1 1], [0.7 0.7 0.7], 'EdgeColor', 'none');
text(7.5, 0.5, 'End', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
patch([15 45 45 15], [0 0 1 1], [0.90 0.40 0.35], 'EdgeColor', 'none');
text(30, 0.5, 'L_2', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
patch([45 65 65 45], [0 0 1 1], [0.95 0.80 0.50], 'EdgeColor', 'none');
text(55, 0.5, 'BG', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
patch([65 85 85 65], [0 0 1 1], [0.45 0.70 0.90], 'EdgeColor', 'none');
text(75, 0.5, 'L_1', 'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

end

%% =============================================================================
% HELPER FUNCTIONS
% =============================================================================

function img = imread_and_prep(filepath)
if exist(filepath, 'file') == 2
    img = imread(filepath);
    if isa(img, 'uint8')
        img = double(img)/255; 
    elseif isa(img, 'logical')
        img = double(img);
    end
    if size(img, 3) > 1, img = rgb2gray(img); end
else
    % Fallback synthetic pattern if file is missing
    [X, Y] = meshgrid(1:200, 1:128);
    img = 0.5 + 0.4 * sin(2 * pi * X / 12);
end
end