function saveFigureFormats(figHandle, savePath)
% saveFigureFormats: High-quality exporter that saves a clean raster PNG,
% a vector-safe PDF, and a classic-traced layered SVG using the print engine.
% gemini 

    if isempty(figHandle) || ~isvalid(figHandle)
        error('Invalid or empty figure handle provided.');
    end
    
    [folder, name, ~] = fileparts(savePath);
    basePath = fullfile(folder, name);

    fprintf('Exporting figures to multi-formats...\n');

    % Ensure figure has vector rendering parameters enabled
    set(figHandle, 'Renderer', 'painters');


    exportgraphics(figHandle, [basePath, '.png'], 'Resolution', 600);
    
    exportgraphics(figHandle, [basePath, '.pdf'], 'ContentType', 'vector');
    
    % 3. FIXED Layered SVG (Bypasses exportgraphics restrictions entirely)
    % Uses the native print driver to map every text, line, and matrix component 
    % into fully editable vector group hierarchies.
    print(figHandle, [basePath, '.svg'], '-dsvg', '-vector');
    
    fprintf('Export complete! Formats saved successfully to:\n   %s[.png, .pdf, .svg]\n', basePath);
end