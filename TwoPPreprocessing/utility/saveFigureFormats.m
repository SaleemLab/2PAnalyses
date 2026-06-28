function saveFigureFormats(figHandle, savePath)
% exporter that saves a  PNG,
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


    exportgraphics(figHandle, [basePath, '.png']);
    
    exportgraphics(figHandle, [basePath, '.pdf'],'ContentType', 'vector');
    
    % Gemini suggested fix Layered SVG 
    % Uses the native print driver to map every text, line, and matrix component 
    % into fully editable vector group hierarchies. 
    print(figHandle, [basePath, '.svg'], '-dsvg');
    
    % Save as a MATLAB .fig file
    savefig(figHandle, [basePath, '.fig']);


    
    fprintf('Export complete! Formats saved successfully to:\n   %s[.png, .pdf, .svg, .fig]\n', basePath);
end