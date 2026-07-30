%%---------------------------PARAMETERS---------------------------------
corridorL = 200; 
corridorH = 12; 
texwidth = 0.04; 
BG_contrast = 0.7; % Background dimming coefficient # previous version was 0.3
basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';
saveDir  = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\VRCartoonTextures-conditions_higherBGContrast\';

%%---------------------------LOAD BASE BACKGROUND 2--------------------
imgBG = imread(fullfile(basePath, 'BG2.jpg'));
if isa(imgBG, 'uint8'), imgBG = double(imgBG)/255; elseif isa(imgBG, 'logical'), imgBG = double(imgBG); end
if size(imgBG, 3) > 1, imgBG = rgb2gray(imgBG); end 
% Midpoint anchor adjusted to 0.75 for a much lighter background gray shade
cleanBG = 0.60 + (imgBG - 0.5) * BG_contrast;
finalBGlength = size(cleanBG, 2);
finalBGheight = size(cleanBG, 1);

%%---------------------------LOAD & SCALE LANDMARKS----------------------
% Vertical Grating
imgG = imread(fullfile(basePath, 'grating_vertical.jpg'));
if isa(imgG, 'uint8'), imgG = double(imgG)/255; elseif isa(imgG, 'logical'), imgG = double(imgG); end
if size(imgG, 3) > 1, imgG = rgb2gray(imgG); end
imgG = histeq(imgG); 

% Plaid
imgP = imread(fullfile(basePath, 'plaid.jpg'));
if isa(imgP, 'uint8'), imgP = double(imgP)/255; elseif isa(imgP, 'logical'), imgP = double(imgP); end
if size(imgP, 3) > 1, imgP = rgb2gray(imgP); end
imgP = histeq(imgP); 

% Resize landmarks to match the expected corridor slice dimensions
final_width = round(finalBGlength * texwidth); 
[u, v] = meshgrid(1:size(imgP,2), 1:size(imgP,1));
[uq, vq] = meshgrid(linspace(1, size(imgP,2), final_width), linspace(1, size(imgP,1), finalBGheight));
grating_res = interp2(u, v, imgG, uq, vq);
plaid_res   = interp2(u, v, imgP, uq, vq);

%%---------------------------DEFINE ALL 6 EXACT LAYOUTS------------------
centers_px = [0.20, 0.40, 0.60, 0.80] * finalBGlength;
condList = struct('name', {}, 'type', {}, 'activeIndices', {}, 'textures', {}, 'filename', {});

% 1. Baseline Layout (G - P - G - P)
condList(1).name          = 'Baseline';
condList(1).type          = 'loop';
condList(1).activeIndices = 1:4;
condList(1).textures      = {grating_res, plaid_res, grating_res, plaid_res};
condList(1).filename      = 'corridor_cartoon_baseline';

% 2. Swap 2 and 3 (G - G - P - P)
condList(2).name          = 'Swap 2 3';
condList(2).type          = 'loop';
condList(2).activeIndices = 1:4;
condList(2).textures      = {grating_res, grating_res, plaid_res, plaid_res};
condList(2).filename      = 'corridor_cartoon_swap23';

% 3. Swap 3 and 4 (G - P - P - G)
condList(3).name          = 'Swap 3 4';
condList(3).type          = 'loop';
condList(3).activeIndices = 1:4;
condList(3).textures      = {grating_res, plaid_res, plaid_res, grating_res};
condList(3).filename      = 'corridor_cartoon_swap34';

% 4. Omit Landmark 2 (G - X - G - P)
condList(4).name          = 'Omit 2';
condList(4).type          = 'loop';
condList(4).activeIndices = [1, 3, 4];
condList(4).textures      = {grating_res, [], grating_res, plaid_res}; % Index 2 left empty
condList(4).filename      = 'corridor_cartoon_omit2';

% 5. Omit Landmark 3 (G - P - X - P)
condList(5).name          = 'Omit 3';
condList(5).type          = 'loop';
condList(5).activeIndices = [1, 2, 4];
condList(5).textures      = {grating_res, plaid_res, [], plaid_res};    % Index 3 left empty
condList(5).filename      = 'corridor_cartoon_omit3';

% 6. Omit Landmark 4 (G - P - G - X)
condList(6).name          = 'Omit 4';
condList(6).type          = 'loop';
condList(6).activeIndices = [1, 2, 3];
condList(6).textures      = {grating_res, plaid_res, grating_res, []};    % Index 4 left empty
condList(6).filename      = 'corridor_cartoon_omit4';

%%---------------------------LOOP, BUILD, AND SAVE PANELS----------------
for cIdx = 1:length(condList)
    % Assign fresh baseline context to completely isolate conditions
    texCorridor = cleanBG;
    indices = condList(cIdx).activeIndices;
    
    for idx = indices
        c = round(centers_px(idx));
        start_col = c - round(final_width/2) + 1;
        
        % Populate using your exact structural mappings
        texCorridor(:, start_col:start_col+final_width-1) = condList(cIdx).textures{idx};
    end
    
    % Initialize figure panel with properties matching your thesis standard
    figSingle = figure('Name', condList(cIdx).name, 'Position', [100, 100, 1000, 250], 'Color', 'w');
    imagesc([0, 200], [0, 12], texCorridor, [0 1]);
    colormap gray; 
    axis tight; 
    box off;
    
    % Format scientific ticks and axes labels
    xlabel('Position (cm)', 'FontSize', 12, 'FontName', 'Arial');
    title(condList(cIdx).name, 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
    
    set(gca, 'XTick', [40, 80, 120, 160], ...
             'XTickLabel', {'40', '80', '120', '160'}, ...
             'YColor', 'none', 'FontSize', 11, 'FontName', 'Arial', 'TickDir', 'out');
         
    % Export to target path
    saveFigureFormats(figSingle, fullfile(saveDir, condList(cIdx).filename));
    close(figSingle);
end

disp('All 6 corrected corridor cartoons exported successfully.');