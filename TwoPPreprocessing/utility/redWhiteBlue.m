function cmap = redWhiteBlue(n)
% Returns an n-point red-white-blue diverging colormap
% Low = blue, middle = white, high = red

if nargin < 1
    n = 256;
end

half = floor(n/2);
r = [linspace(0,1,half)'; ones(n-half,1)];
g = [linspace(0,1,half)'; linspace(1,0,n-half)'];
b = [ones(half,1); linspace(1,0,n-half)'];

cmap = [r g b];
end
