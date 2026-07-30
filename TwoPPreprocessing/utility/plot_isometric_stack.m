function fig = plot_isometric_stack(images, z_spacing, rot_deg)
% PLOT_ISOMETRIC_STACK Visualizes multi-plane imaging FOVs as horizontal 3D slabs.
%
% images    : cell array of 2D image matrices
% z_spacing : distance between planes along vertical axis (default: 0.3)
% rot_deg   : rotation in degrees: 0, 90, 180, or 270 (default: 180)

if nargin < 2, z_spacing = 0.3; end
if nargin < 3, rot_deg = 180; end % Default to 180 to flip back-to-front

num_planes = numel(images);
fig = figure('Color', 'w', 'Position', [100, 100, 700, 600]);
ax = axes('Parent', fig); 
hold(ax, 'on');

% Apply rotation to the first image to establish spatial dimensions
k = round(rot_deg / 90); % Convert degrees (90, 180, 270) to rot90 steps
img1_rot = rot90(im2double(images{1}), k);

[h, w] = size(img1_rot);
[X, Y] = meshgrid(1:w, 1:h);

for p = 1:num_planes
    img = im2double(images{p});
    
    % Rotate the image matrix before rendering
    img = rot90(img, k);
    
    % Z height for this plane
    Z = ones(h, w) * (-p * z_spacing * h);
    
    % Draw image as a 3D surface
    surf(ax, X, Y, Z, img, ...
        'FaceColor', 'texturemap', ...
        'EdgeColor', 'none', ...
        'CDataMapping', 'scaled');
end

% Set gray colormap and standard camera view
colormap(ax, 'gray');
axis(ax, 'equal', 'off');

% Pitch the view down so planes look like tilted slabs
view(ax, [-37.5, 30]); 

hold(ax, 'off');
end