%% Local function
function xt = compressAxis(x, breakpt, compressFactor)
% compressAxis  Warp x values above breakpt into a compressed range.
% Values at or below breakpt keep a 1:1 mapping.
%   x               values to transform
%   breakpt         x value where compression starts (e.g. 80)
%   compressFactor  how much to squeeze the region above breakpt.
%                   Smaller means more squeeze. 0.25 squeezes that
%                   region to a quarter of its natural width.
    xt = x;
    above = x > breakpt;
    xt(above) = breakpt + (x(above) - breakpt) * compressFactor;
end