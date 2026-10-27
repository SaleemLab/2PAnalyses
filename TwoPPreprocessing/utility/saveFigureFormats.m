% function saveFigureFormats(figHandle, savePath)
% % exporter that saves a PNG, PDF, SVG, and FIG file.
% 
%     if isempty(figHandle) || ~isvalid(figHandle)
%         error('Invalid or empty figure handle provided.');
%     end
%     
%     [folder, name, ~] = fileparts(savePath);
%     
%     % Ensure target directory exists before saving
%     if ~isempty(folder) && ~exist(folder, 'dir')
%         mkdir(folder);
%     end
%     
%     basePath = fullfile(folder, name);
%     fprintf('Exporting figures to multi-formats...\n');
%     
%     % Check if figure contains 3D surf/texturemap objects
%     hasTextures = ~isempty(findobj(figHandle, 'Type', 'surface'));
%     
%     if hasTextures
%         % OpenGL is required for 3D surfaces and texturemaps
%         set(figHandle, 'Renderer', 'opengl');
%         
%         % Export PNG (Raster)
%         exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
%         
%         % PDF/SVG: Capture textured surface cleanly
%         exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'image', 'Resolution', 300);
%         print(figHandle, [basePath, '.svg'], '-dsvg');
%     else
%         % Standard 2D vector plots
%         set(figHandle, 'Renderer', 'painters');
%         exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
%         exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'vector');
%         print(figHandle, [basePath, '.svg'], '-dsvg');
%     end
%     
%     % Save as a MATLAB .fig file
%     savefig(figHandle, [basePath, '.fig']);
%     
%     fprintf('Export complete! Formats saved successfully to:\n   %s[.png, .pdf, .svg, .fig]\n', basePath);
% end


function saveFigureFormats(figHandle, savePath)
% exporter that saves a PNG, PDF, SVG, and FIG file.
if isempty(figHandle) || ~isvalid(figHandle)
        error('Invalid or empty figure handle provided.');
end
    [folder, name, ~] = fileparts(savePath);
% Ensure target directory exists before saving
if ~isempty(folder) && ~exist(folder, 'dir')
        mkdir(folder);
end
    basePath = fullfile(folder, name);
    fprintf('Exporting figures to multi-formats...\n');
% Distinguish a genuinely 3D surface (needs OpenGL/raster — a flat vector
% export can't represent real lighting/perspective) from a FLAT surface
% object viewed top-down, like pcolor heatmaps (view elevation ~90°).
% pcolor creates Type='surface' even though it's a 2D heatmap, so object
% type alone isn't enough — check the axes' actual view angle. Flat
% surfaces get exported as vector, so each cell stays an individually
% editable/movable/recolourable shape in Illustrator/Inkscape, matching
% how a normal 2D vector plot is handled below. A true 2D image object
% (imagesc) is always raster regardless, nothing to preserve as vector.
    axesList = findall(figHandle, 'Type', 'axes');
    has3DSurface = false;
    for iAx = 1:numel(axesList)
        surfObjs = findobj(axesList(iAx), 'Type', 'surface');
        if ~isempty(surfObjs)
            viewAngles = get(axesList(iAx), 'View'); % [azimuth, elevation]
            if abs(viewAngles(2)) < 89 % elevation away from top-down => truly 3D
                has3DSurface = true;
            end
        end
    end
    hasImage = ~isempty(findobj(figHandle, 'Type', 'image'));

if has3DSurface
% OpenGL is required for genuine 3D surfaces and texturemaps
        set(figHandle, 'Renderer', 'opengl');
        exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
        exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'image', 'Resolution', 300);
        print(figHandle, [basePath, '.svg'], '-dsvg');
elseif hasImage
% Flat 2D image content (e.g. imagesc) — always raster, painters
% renderer avoids the OpenGL tiling artifact seen at high resolution.
        set(figHandle, 'Renderer', 'painters');
        exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
        exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'image', 'Resolution', 300);
        print(figHandle, [basePath, '.svg'], '-dsvg');
else
% Standard 2D vector plots — includes flat top-down pcolor heatmaps, so
% each cell remains an individually editable shape in the saved PDF/SVG.
        set(figHandle, 'Renderer', 'painters');
        exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
        exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'vector');
        print(figHandle, [basePath, '.svg'], '-dsvg');
end
% Save as a MATLAB .fig file
    savefig(figHandle, [basePath, '.fig']);
    fprintf('Export complete! Formats saved successfully to:\n   %s[.png, .pdf, .svg, .fig]\n', basePath);
end