%% plot_corridor_figure.m
% Fully self-contained script to generate the corridor texture figure
% Panels:
%   A - Base corridor (BG2, G-P-G-P) with annotations
%   B - Background periodicity (BG2, no landmarks, shaded chunks)
%   C - All 4 backgrounds in standard G-P-G-P configuration
%   D - Landmark manipulation conditions (BG2)

%%---------------------------PARAMETERS---------------------------------
corridorL    = 200;
corridorH    = 12;
texwidth     = 0.04;
fig_texwidth = 0.02;

basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';

%%---------------------------LOAD TEXTURES------------------------------
textures = struct('matrix', []);

bgNames = {'BG1.jpg', 'BG2.jpg', 'BG3.jpg', 'BG4.jpg'};
for k = 1:4
    img = imread(fullfile(basePath, bgNames{k}));
    if isa(img, 'uint8'),      img = double(img)/255;
    elseif isa(img,'logical'), img = double(img); end
    if size(img,3) > 1, img = rgb2gray(img); end
    textures(k+1).matrix = img;
end

landmarkFiles = {'grey.jpg', 'grating_vertical.jpg', 'plaid.jpg'};
landmarkIdx   = [1, 6, 8];
for i = 1:length(landmarkFiles)
    img = imread(fullfile(basePath, landmarkFiles{i}));
    if isa(img, 'uint8'),      img = double(img)/255;
    elseif isa(img,'logical'), img = double(img); end
    if size(img,3) > 1, img = rgb2gray(img); end
    textures(landmarkIdx(i)).matrix = img;
end

%%---------------------------DERIVED QUANTITIES-------------------------
finalBGlength   = size(textures(2).matrix, 2);
finalBGheight   = size(textures(2).matrix, 1);
final_width     = round(finalBGlength * texwidth);
fig_final_width = round(finalBGlength * fig_texwidth);
centers_px      = [0.20, 0.40, 0.60, 0.80] * finalBGlength;

% Resize landmarks to experiment width
grating_src = textures(6).matrix;
plaid_src   = textures(8).matrix;
[u,  v]  = meshgrid(1:size(grating_src,2), 1:size(grating_src,1));
[uq, vq] = meshgrid(linspace(1,size(grating_src,2),final_width), ...
                    linspace(1,size(grating_src,1),finalBGheight));
grating_res = interp2(u, v, grating_src, uq, vq);
plaid_res   = interp2(u, v, plaid_src,   uq, vq);

% Resize to figure display width
[u2,  v2]  = meshgrid(1:size(grating_res,2), 1:size(grating_res,1));
[uq2, vq2] = meshgrid(linspace(1,size(grating_res,2),fig_final_width), ...
                      linspace(1,size(grating_res,1),finalBGheight));
fig_grating = interp2(u2, v2, grating_res, uq2, vq2);
fig_plaid   = interp2(u2, v2, plaid_res,   uq2, vq2);

%%---------------------------FIGURE LAYOUT------------------------------
fig    = figure('Name','Corridor Figure','Position',[50 50 1500 950]);
left_m = 0.14;  right_m = 0.02;
top_m  = 0.05;  bot_m   = 0.06;
gap    = 0.025;
n_rows = 1 + 1 + 4 + 6;
strip_h = (1 - top_m - bot_m - gap*3) / n_rows;
W       = 1 - left_m - right_m;

standard = [1 2 1 2];
bg1      = textures(2).matrix;
bg2      = textures(3).matrix;

%%---------------------------PANEL A-----------------------------------
texA = place_landmarks_fn(bg2, standard, centers_px, fig_final_width, fig_grating, fig_plaid);
yA   = 1 - top_m - strip_h;
axA  = draw_strip(fig, [left_m, yA, W, strip_h*0.85], texA, [0 200]);
set(axA, 'XTick',[40 80 120 160],'XTickLabel',{'40','80','120','160'},...
    'XColor','k','FontSize',9);
xlabel(axA, 'Position (cm)', 'FontSize', 10);

lm_labels = {'G','P','G','P'};
pos_cm    = [40 80 120 160];
for ii = 1:4
    text(axA, pos_cm(ii), 14, lm_labels{ii}, ...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color','k');
end

x1n = left_m + (40/200)*W;
x2n = left_m + (80/200)*W;
yn  = yA - 0.015;
annotation(fig,'doublearrow',[x1n x2n],[yn yn],'Head1Length',6,'Head2Length',6,'Head1Width',6,'Head2Width',6,'Color','k');
annotation(fig,'textbox',[(x1n+x2n)/2-0.02, yn-0.025, 0.04, 0.02],...
    'String','40 cm','EdgeColor','none','HorizontalAlignment','center','FontSize',8);

text(axA, -25, 6, 'A', 'FontSize', 13, 'FontWeight', 'bold');

%%---------------------------PANEL B-----------------------------------
yB  = yA - strip_h - gap;
axB = draw_strip(fig, [left_m, yB, W, strip_h*0.85], bg2, [0 200]);
set(axB, 'XTick',[],'XColor','none');

% Auto-detect identical 24cm repeats in the background texture
chunk_cm = 24;
chunk_px = round(chunk_cm / corridorL * finalBGlength);
n_windows = floor(finalBGlength / chunk_px);

% Extract each full 24cm chunk
chunks = zeros(finalBGheight, chunk_px, n_windows);
for ww = 1:n_windows
    c1 = (ww-1)*chunk_px + 1;
    c2 = ww*chunk_px;
    chunks(:,:,ww) = bg2(:, c1:c2);
end

% Find first chunk that repeats (use as reference)
tol = 0.01;
ref = chunks(:,:,1);
identical = false(1, n_windows);
for ww = 1:n_windows
    diff_rms = sqrt(mean(mean((chunks(:,:,ww) - ref).^2)));
    identical(ww) = diff_rms < tol;
end
% If first chunk is partial, re-reference to first identical match
first_full = find(identical, 1, 'first');
if first_full > 1
    ref = chunks(:,:,first_full);
    for ww = 1:n_windows
        diff_rms = sqrt(mean(mean((chunks(:,:,ww) - ref).^2)));
        identical(ww) = diff_rms < tol;
    end
end

% Shade identical chunks
for ww = 1:n_windows
    if identical(ww)
        x1 = (ww-1)*chunk_cm;
        x2 = ww*chunk_cm;
        patch(axB,[x1 x2 x2 x1],[0 0 12 12],[0.2 0.5 0.9],...
            'FaceAlpha',0.25,'EdgeColor','none');
    end
end
n_identical = sum(identical);
text(axB, 100, 15, sprintf('%d identical 24 cm background repeats', n_identical), ...
    'HorizontalAlignment','center','FontSize',8,'Color',[0.1 0.3 0.7]);
text(axB, -25, 6, 'B', 'FontSize', 13, 'FontWeight', 'bold');

%%---------------------------PANEL C-----------------------------------
bg_indices = [2 3 4 5];
bg_labels  = {'BGLeft','BGRight','BGCeil','BGFloor'};
yC_top     = yB - strip_h - gap;

for kk = 1:4
    bg  = textures(bg_indices(kk)).matrix;
    tex = place_landmarks_fn(bg, standard, centers_px, fig_final_width, fig_grating, fig_plaid);
    yC  = yC_top - (kk-1)*strip_h;
    axC = draw_strip(fig, [left_m, yC, W, strip_h*0.85], tex, [0 200]);
    set(axC,'XTick',[],'XColor','none');
    text(axC, -25, 6, bg_labels{kk}, 'FontSize',8,'HorizontalAlignment','right');
    if kk == 1
        text(axC, -25, 15, 'C', 'FontSize',13,'FontWeight','bold');
    end
end

%%---------------------------PANEL D-----------------------------------
conditions = {
    [1 2 1 2], 'G - P - G - P';
    [1 1 2 2], 'G - G - P - P';
    [1 2 2 1], 'G - P - P - G';
    [1 0 1 2], 'G - o - G - P';
    [1 2 0 2], 'G - P - o - P';
    [1 2 1 0], 'G - P - G - o';
};

yD_top = yC_top - 4*strip_h - gap;

for kk = 1:6
    lm_types = conditions{kk,1};
    lbl      = conditions{kk,2};
    tex = place_landmarks_fn(bg1, lm_types, centers_px, fig_final_width, fig_grating, fig_plaid);
    yD  = yD_top - (kk-1)*strip_h;
    axD = draw_strip(fig, [left_m, yD, W, strip_h*0.85], tex, [0 200]);

    if kk == 6
        set(axD,'XTick',[40 80 120 160],'XTickLabel',{'40','80','120','160'},...
            'XColor','k','FontSize',9);
        xlabel(axD,'Position (cm)','FontSize',10);
    else
        set(axD,'XTick',[],'XColor','none');
    end

    text(axD, -25, 6, lbl, 'FontSize',8,'HorizontalAlignment','right');
    if kk == 1
        text(axD, -25, 15, 'D', 'FontSize',13,'FontWeight','bold');
    end
end


%
outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\Methods\Fig2.3_CorridorTextures\';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    saveFigureFormats(fig, fullfile(outputDir, 'Textures_base_manip_allwalls'));

    
    
    
 

%%---------------------------LOCAL FUNCTIONS----------------------------
function tex = place_landmarks_fn(bg, landmark_types, centers, lw, G, P)
    tex = bg;
    for ii = 1:length(landmark_types)
        if landmark_types(ii) == 0, continue; end
        c  = round(centers(ii));
        sc = max(c - round(lw/2) + 1, 1);
        ec = min(sc + lw - 1, size(tex,2));
        if landmark_types(ii) == 1
            tex(:, sc:ec) = G(:, 1:(ec-sc+1));
        elseif landmark_types(ii) == 2
            tex(:, sc:ec) = P(:, 1:(ec-sc+1));
        end
    end
end

function ax = draw_strip(fig, pos, tex, xlims)
    ax = axes(fig, 'Position', pos);
    imagesc(xlims, [0 12], tex, [0 1]);
    colormap(ax, gray);
    axis tight; box off;
    set(ax, 'YColor','none','YTick',[]);
end