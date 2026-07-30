function plotOrientedTextures(filename, nY, nX)
% PLOTORIENTEDTEXTURES Reads a horizontal texture image, reshapes it inline,
% and plots it vertically with correct spatial orientation.
%
% Usage:
%   plotOrientedTextures('BG1.jpg', 16, 16)

    % 1. Read image and convert to 2D double matrix inline
    img = double(imread(filename));
    if ndims(img) == 3, img = img(:,:,1); end % Keep single channel if RGB

    % 2. Reshape into 3D volume [nY x nX x nUnits] and flatten vertically inside imagesc
    figure('Name', filename);
    imagesc(reshape(reshape(img, [nY, nX, []]), [], nX));
    
    % 3. Format visual properties
    colormap gray; 
    axis image; % Preserves 1:1 pixel aspect ratio
    xlabel('Spatial X');
    ylabel('Unit Index \times Spatial Y');
    title(filename, 'Interpreter', 'none');

end