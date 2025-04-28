function cmap = whiteBlueRed(N)
% WHITEBLUERED  Generate a white→blue→red colormap
%
%   CMAP = WHITEBLUERED(N) returns an N-by-3 array of RGB values
%   that smoothly transition from white to blue in the first half
%   and from blue to red in the second half. If N is omitted, it
%   defaults to 256.
%
% Example:
%   colormap(whiteBlueRed(128));
%   colorbar;

if nargin<1 || isempty(N)
    N = 256;
end
half = floor(N/2);
cmap = zeros(N,3);

% white -> blue
cmap(1:half,1) = linspace(1, 0, half)';   % R: 1->0
cmap(1:half,2) = linspace(1, 0, half)';   % G: 1->0
cmap(1:half,3) = 1;                       % B: 1

% blue -> red
rem = N - half;
cmap(half+1:end,1) = linspace(0, 1, rem)';  % R: 0->1
cmap(half+1:end,2) = 0;                     % G: 0
cmap(half+1:end,3) = linspace(1, 0, rem)';  % B: 1->0
end
