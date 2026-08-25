%% Background texture (BG2), one repeat cut into 6 chunks
% Standalone script - loads BG2.jpg directly and splits one period into
% the 6 chunks defined by BGseg_pct in your GLM code.

%% Parameters
basePath = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\VRCorridorFinal-TextureCopyForSonali\';
saveSubfolder = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter1-VISpSomas\Fig2_7_Section3';
nameBase = 'BG2_6chunks';

corridorL = 200;      % corridor length (cm)
period_pct = 12;      % one BG repeat, % of corridor length
n_chunks = 6;          % chunks per repeat

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

%% Cut one period into 6 chunks
chunk_edges_pct = linspace(0, period_pct, n_chunks + 1);
chunks = cell(1, n_chunks);
for i = 1:n_chunks
    s_px = round(chunk_edges_pct(i)/100 * corridorL * px_per_cm) + 1;
    e_px = round(chunk_edges_pct(i+1)/100 * corridorL * px_per_cm);
    chunks{i} = img(:, s_px:e_px);
end

%% Colours (red -> yellow gradient, matching your existing BG convention)
col_start = [247 33 10]/255;
col_end   = [247 242 143]/255;
bg_cols = [linspace(col_start(1), col_end(1), n_chunks)', ...
           linspace(col_start(2), col_end(2), n_chunks)', ...
           linspace(col_start(3), col_end(3), n_chunks)'];

%% Plot
figB = figure('Name', 'BG2 - 6 chunks', 'Position', [100 300 1200 250]);
for i = 1:n_chunks
    subplot(1, n_chunks, i);
    imagesc(chunks{i}, [0 1]);
    colormap gray; axis image off;
    title(sprintf('BG%d', i), 'FontSize', 11, 'FontWeight', 'bold', 'Color', bg_cols(i,:));
end
sgtitle('Background texture (BG2), one repeat cut into 6 chunks', 'FontSize', 13);

%% Save
saveFigureFormats(figB, fullfile(saveSubfolder, nameBase));