%%---------------------------PARAMETERS---------------------------------
corridorL = 200; 
corridorH = 12; 
texwidth = 0.04; 
BG_contrast = 0.3; % Keeps the background at your desired 0.3 dim level
basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';

%%---------------------------LOAD BASE BACKGROUND 2--------------------
% Loading only BG2 directly
imgBG = imread(fullfile(basePath, 'BG2.jpg'));
if isa(imgBG, 'uint8'), imgBG = double(imgBG)/255; elseif isa(imgBG, 'logical'), imgBG = double(imgBG); end
if size(imgBG, 3) > 1, imgBG = rgb2gray(imgBG); end 

% CHANGED: Swapped the midpoint anchor from 0.5 to 0.75 to make the gray a much lighter shade
texCorridor = 0.60 + (imgBG - 0.5) * BG_contrast;

finalBGlength = size(texCorridor, 2);
finalBGheight = size(texCorridor, 1);

%%---------------------------LOAD & SCALE LANDMARKS----------------------
% Vertical Grating
imgG = imread(fullfile(basePath, 'grating_vertical.jpg'));
if isa(imgG, 'uint8'), imgG = double(imgG)/255; elseif isa(imgG, 'logical'), imgG = double(imgG); end
if size(imgG, 3) > 1, imgG = rgb2gray(imgG); end
% Maximize landmark contrast using histogram equalization
imgG = histeq(imgG); 

% Plaid
imgP = imread(fullfile(basePath, 'plaid.jpg'));
if isa(imgP, 'uint8'), imgP = double(imgP)/255; elseif isa(imgP, 'logical'), imgP = double(imgP); end
if size(imgP, 3) > 1, imgP = rgb2gray(imgP); end
% Maximize landmark contrast using histogram equalization
imgP = histeq(imgP);

% Resize landmarks to match the expected corridor slice dimensions
final_width = round(finalBGlength * texwidth); 
[u, v] = meshgrid(1:size(imgP,2), 1:size(imgP,1));
[uq, vq] = meshgrid(linspace(1, size(imgP,2), final_width), linspace(1, size(imgP,1), finalBGheight));
grating_res = interp2(u, v, imgG, uq, vq);
plaid_res   = interp2(u, v, imgP, uq, vq);

%%---------------------------SUPERIMPOSE LANDMARKS----------------------
centers_px = [0.20, 0.40, 0.60, 0.80] * finalBGlength;
for idx = 1:4
    c = round(centers_px(idx));
    start_col = c - round(final_width/2) + 1;
    
    % Alternate Grating and Plaid onto the dim background frame
    if mod(idx,2) == 1
        texCorridor(:, start_col:start_col+final_width-1) = grating_res;
    else
        texCorridor(:, start_col:start_col+final_width-1) = plaid_res;
    end
end

%%---------------------------PLOT AND EXPORT PANEL----------------------
figSingle = figure('Name', 'Max Contrast Landmarks Panel', 'Position', [100, 100, 1000, 250]);
imagesc([0, 200], [0, 12], texCorridor, [0 1]);
colormap gray; 
axis tight; 
box off;

% Format clean scientific labels
xlabel('Position (cm)', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'XTick', [40, 80, 120, 160, 200], ...
         'XTickLabel', {'40', '80', '120', '160', '200'}, ...
         'YColor', 'none', 'FontSize', 11, 'FontName', 'Arial', 'TickDir', 'out');

saveFigureFormats(figSingle, 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\VRTexture_Cartoons\corridor_cartoon');