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
    
    % Check if figure contains 3D surf/texturemap objects
    hasTextures = ~isempty(findobj(figHandle, 'Type', 'surface'));
    
    if hasTextures
        % OpenGL is required for 3D surfaces and texturemaps
        set(figHandle, 'Renderer', 'opengl');
        
        % Export PNG (Raster)
        exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
        
        % PDF/SVG: Capture textured surface cleanly
        exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'image', 'Resolution', 300);
        print(figHandle, [basePath, '.svg'], '-dsvg');
    else
        % Standard 2D vector plots
        set(figHandle, 'Renderer', 'painters');
        exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 300);
        exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'vector');
        print(figHandle, [basePath, '.svg'], '-dsvg');
    end
    
    % Save as a MATLAB .fig file
    savefig(figHandle, [basePath, '.fig']);
    
    fprintf('Export complete! Formats saved successfully to:\n   %s[.png, .pdf, .svg, .fig]\n', basePath);
end
% function saveFigureFormats(figHandle, savePath)
% % exporter that saves a  PNG,
% % a vector-safe PDF, and a classic-traced layered SVG using the print engine.
% % gemini 
% 
%     if isempty(figHandle) || ~isvalid(figHandle)
%         error('Invalid or empty figure handle provided.');
%     end
%     
%     [folder, name, ~] = fileparts(savePath);
%     basePath = fullfile(folder, name);
% 
%     fprintf('Exporting figures to multi-formats...\n');
% 
%     % Ensure figure has vector rendering parameters enabled
%     set(figHandle, 'Renderer', 'painters');
% 
% 
%     exportgraphics(figHandle, [basePath, '.png']);
%     
%     exportgraphics(figHandle, [basePath, '.pdf'],'ContentType', 'vector');
%     
%     % Gemini suggested fix Layered SVG 
%     % Uses the native print driver to map every text, line, and matrix component 
%     % into fully editable vector group hierarchies. 
%     print(figHandle, [basePath, '.svg'], '-dsvg');
%     
%     % Save as a MATLAB .fig file
%     savefig(figHandle, [basePath, '.fig']);
% 
% 
%     
%     fprintf('Export complete! Formats saved successfully to:\n   %s[.png, .pdf, .svg, .fig]\n', basePath);
% end