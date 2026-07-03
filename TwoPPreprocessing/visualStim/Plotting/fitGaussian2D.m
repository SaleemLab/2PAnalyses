function [fitFun, fwhmX, fwhmY] = fitGaussian2D(x, y, Z)
% fitGaussian2D: Fits
%   Z = a*exp(-((X-x0)^2/(2*sx^2) + (Y-y0)^2/(2*sy^2))) + d
% via lsqnonlin and returns a handle to the fitted surface plus FWHM in x/y.
%
% Requires Optimization Toolbox (lsqnonlin).
%
% INPUTS:
%   x - vector of x-axis grid values (e.g. azimuth offset, degrees)
%   y - vector of y-axis grid values (e.g. elevation offset, degrees)
%   Z - 2D response matrix, size [length(y), length(x)]
%
% OUTPUTS:
%   fitFun - function handle, fitFun(xq, yq) evaluates the fitted surface
%            at meshgrid-style query points (xq, yq must be same size)
%   fwhmX  - full width at half maximum along x
%   fwhmY  - full width at half maximum along y
%
% USAGE:
%   [fit2D, fwhmX, fwhmY] = fitGaussian2D(azOffset, elOffset, meanGrid2D);

    [X, Y] = meshgrid(x, y);
    validMask = ~isnan(Z);

    gauss2DEq = @(p, xy) p(1) * exp(-(((xy(:,:,1) - p(2)).^2)/(2*p(4)^2) + ...
                                       ((xy(:,:,2) - p(3)).^2)/(2*p(5)^2))) + p(6);

    xy = cat(3, X, Y);
    [~, peakIdx] = max(Z(:), [], 'omitnan');
    [rPk, cPk] = ind2sub(size(Z), peakIdx);
    p0 = [max(Z(:), [], 'omitnan') - min(Z(:), [], 'omitnan'), ...
          x(cPk), y(rPk), range(x)/4, range(y)/4, min(Z(:), [], 'omitnan')];
    lb = [0, min(x), min(y), 0.1, 0.1, -Inf];
    ub = [Inf, max(x), max(y), range(x)*2, range(y)*2, Inf];

    ZFit = Z; ZFit(~validMask) = 0;  % lsqnonlin can't handle NaNs directly

    opts = optimoptions('lsqnonlin', 'Display', 'off');
    weights = double(validMask);
    costFun = @(p) weights .* (gauss2DEq(p, xy) - ZFit);
    pFit = lsqnonlin(costFun, p0, lb, ub, opts);

    fitFun = @(xq, yq) gauss2DEq(pFit, cat(3, xq, yq));
    fwhmX  = 2.3548 * abs(pFit(4));
    fwhmY  = 2.3548 * abs(pFit(5));
end
