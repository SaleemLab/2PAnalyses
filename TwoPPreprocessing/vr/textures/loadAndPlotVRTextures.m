%%---------------------------PARAMETERS---------------------------------
corridorL = 200; 
corridorH = 12; 
texwidth = 0.04; 
BG_contrast = 0.5; 
% Define your folder path here
basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';

%%---------------------------LOAD TEXTURES-------------------------------
textures = struct('matrix', []);
% Load Backgrounds
bgNames = {'BG1.jpg', 'BG2.jpg', 'BG3.jpg', 'BG4.jpg'};
for k = 1:4
    img = imread(fullfile(basePath, bgNames{k}));
    
    % Corrected Type Checking
    if isa(img, 'uint8')
        img = double(img)/255; 
    elseif isa(img, 'logical')
        img = double(img);
    end
    
    if size(img, 3) > 1, img = rgb2gray(img); end 
    textures(k+1).matrix = img;
end

% Load Specific Landmarks
landmarkFiles = {'grey.jpg', 'grating_vertical.jpg', 'plaid.jpg'};
landmarkIdx   = [1, 6, 8]; 
for i = 1:length(landmarkFiles)
    img = imread(fullfile(basePath, landmarkFiles{i}));
    
    if isa(img, 'uint8')
        img = double(img)/255; 
    elseif isa(img, 'logical')
        img = double(img);
    end
    
    if size(img, 3) > 1, img = rgb2gray(img); end
    textures(landmarkIdx(i)).matrix = img;
end
% Create horizontal grating from vertical for the subplot
textures(7).matrix = textures(6).matrix';

%%---------------------------PLOTTING SECTION----------------------------
finalBGlength = size(textures(2).matrix, 2);
finalBGheight = size(textures(2).matrix, 1);

% --- FIGURE 1: Catalog ---
fig1 = figure('Name', 'Loaded Textures Catalog');
titles = {'BG1', 'BG2', 'BG3', 'BG4', 'Vertical Grating', 'Horizontal Grating', 'Plaid'};
plotIdx = [2, 3, 4, 5, 6, 7, 8];
for i = 1:7
    subplot(7,1,i)
    tex = textures(plotIdx(i)).matrix;
    imagesc(tex, [0 1]);
    title([titles{i}, ' (', num2str(size(tex,2)), 'x', num2str(size(tex,1)), ')'])
    colormap gray; axis equal; axis off;
end

% --- FIGURE 2: Virtual Corridor ---
fig2 = figure('Name', 'Final Corridor Layout');
grating_src = textures(6).matrix;
plaid_src = textures(8).matrix;
final_width = round(finalBGlength * texwidth); 
[u, v] = meshgrid(1:size(plaid_src,2), 1:size(plaid_src,1));
[uq, vq] = meshgrid(linspace(1, size(plaid_src,2), final_width), linspace(1, size(plaid_src,1), finalBGheight));
grating_res = interp2(u, v, grating_src, uq, vq);
plaid_res = interp2(u, v, plaid_src, uq, vq);

for k = 2:5
    subplot(4,1,k-1)
    tex = textures(k).matrix;
    centers_px = [0.20, 0.40, 0.60, 0.80] * finalBGlength;
    
    for idx = 1:4
        c = round(centers_px(idx));
        start_col = c - round(final_width/2) + 1;
        
        if mod(idx,2) == 1
            tex(:, start_col:start_col+final_width-1) = grating_res;
        else
            tex(:, start_col:start_col+final_width-1) = plaid_res;
        end
    end
    
    imagesc([0, 200], [0, 12], tex, [0 1]);
    ylabel('');
    set(gca, 'YColor', 'none');

    if (k-1) == 4
        xlabel('Position (cm)');
        set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'FontSize', 11); 
        set(gca, 'XTick', [40, 80, 120, 160, 200], ...
                 'XTickLabel', {'40', '80', '120', '160', '200'}, ...
                 'XColor', 'k');
    else
        xlabel('');
        set(gca, 'XTick', [], 'XTickLabel', []); 
        set(gca, 'XColor', 'none');             
    end
    
    colormap gray; axis tight; box off;
end

% --- EXPORT FIGURE 1 & 2 (JPG & SVG) ---
exportgraphics(fig1, fullfile(basePath, 'Figure1_Catalog.jpg'), 'Resolution', 300);
exportgraphics(fig1, fullfile(basePath, 'Figure1_Catalog.svg'), 'ContentType', 'vector');

exportgraphics(fig2, fullfile(basePath, 'Figure2_Corridors.jpg'), 'Resolution', 300);
exportgraphics(fig2, fullfile(basePath, 'Figure2_Corridors.svg'), 'ContentType', 'vector');

%%---------------------------SETUP BACKGROUND 2--------------------------
targetBG = textures(3).matrix;
centers_px = [0.20, 0.40, 0.60, 0.80] * finalBGlength;

% Create figure
fig3 = figure('Name', 'Experimental Conditions - BG2', 'Position', [100, 50, 1200, 1100]);

%% 1. Original Layout (G - P - G - P)
subplot(6,1,1)
tex = targetBG; 
for idx = 1:4
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    if mod(idx,2) == 1, tex(:, start_col:start_col+final_width-1) = grating_res;
    else, tex(:, start_col:start_col+final_width-1) = plaid_res; end
end
imagesc([0, 200], [0, 12], tex, [0 1]);
colormap gray; axis tight; box off; 
set(gca, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []); 

%% 2. Swap 2 and 3 (G - G - P - P)
subplot(6,1,2)
tex = targetBG; 
types = {grating_res, grating_res, plaid_res, plaid_res}; 
for idx = 1:4
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    tex(:, start_col:start_col+final_width-1) = types{idx};
end
imagesc([0, 200], [0, 12], tex, [0 1]);
colormap gray; axis tight; box off; 
set(gca, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);

%% 3. Swap 3 and 4 (G - P - P - G) - THE BOTTOM PLOT
subplot(6,1,3)
tex = targetBG; 
types = {grating_res, plaid_res, plaid_res, grating_res}; 
for idx = 1:4
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    tex(:, start_col:start_col+final_width-1) = types{idx};
end
imagesc([0, 200], [0, 12], tex, [0 1]);

set(gca, 'YColor', 'none', 'YTick', [], 'YTickLabel', []);
ylabel('');
set(gca, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);
colormap gray; axis tight; box off;

%% 2. Omit Landmark 2 (G - X - G - P)
subplot(6,1,4)
tex = targetBG; 
for idx = [1, 3, 4] 
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    if mod(idx,2) == 1, tex(:, start_col:start_col+final_width-1) = grating_res;
    else, tex(:, start_col:start_col+final_width-1) = plaid_res; end
end
imagesc([0, 200], [0, 12], tex, [0 1]);
colormap gray; axis tight; box off; 
set(gca, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);

%% 3. Omit Landmark 3 (G - P - X - P)
subplot(6,1,5)
tex = targetBG; 
for idx = [1, 2, 4] 
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    if mod(idx,2) == 1, tex(:, start_col:start_col+final_width-1) = grating_res;
    else, tex(:, start_col:start_col+final_width-1) = plaid_res; end
end
imagesc([0, 200], [0, 12], tex, [0 1]);
colormap gray; axis tight; box off; 
set(gca, 'XColor', 'none', 'YColor', 'none', 'XTick', [], 'YTick', []);

%% 4. Omit Landmark 4 (G - P - G - X)
subplot(6,1,6)
tex = targetBG; 
for idx = [1, 2, 3] 
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    if mod(idx,2) == 1, tex(:, start_col:start_col+final_width-1) = grating_res;
    else, tex(:, start_col:start_col+final_width-1) = plaid_res; end
end
imagesc([0, 200], [0, 12], tex, [0 1]);
colormap gray; axis tight; box off; 
xlabel('Position (cm)', 'FontSize', 12);
set(gca, 'XColor', 'k', 'FontSize', 11);
set(gca, 'XTick', [40, 80, 120, 160], ...
         'XTickLabel', {'40', '80', '120', '160'}, 'YTick', [], 'YColor', 'none');
ylabel('');


% --- EXPORT BG2 CONDITIONS (JPG & SVG) ---
exportgraphics(fig3, fullfile(basePath, 'BG2_Conditions.jpg'), 'Resolution', 300);
exportgraphics(fig3, fullfile(basePath, 'BG2_Conditions.svg'), 'ContentType', 'vector');