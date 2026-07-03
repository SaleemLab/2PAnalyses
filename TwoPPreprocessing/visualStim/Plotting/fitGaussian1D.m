function [fitFun, fwhm] = fitGaussian1D(x, y)
% fitGaussian1D: Fits y = a*exp(-(x-b)^2/(2*c^2)) + d via lsqcurvefit and
% returns a handle to the fitted function plus its FWHM.
%
% Requires Optimization Toolbox (lsqcurvefit).
%
% INPUTS:
%   x - independent variable (e.g. azimuth or elevation offset, degrees)
%   y - response values (same size as x)
%
% OUTPUTS:
%   fitFun - function handle, fitFun(xq) evaluates the fitted Gaussian at xq
%   fwhm   - full width at half maximum of the fitted Gaussian
%
% USAGE:
%   [azFit, azFWHM] = fitGaussian1D(azOffset, meanAzProfile);
%   xFine = linspace(min(azOffset), max(azOffset), 200);
%   plot(xFine, azFit(xFine));

    x = x(:); y = y(:);
    valid = ~isnan(y);
    x = x(valid); y = y(valid);

    gaussEq = @(p, x) p(1) * exp(-((x - p(2)).^2) / (2*p(3)^2)) + p(4);

    [~, peakIdx] = max(y);
    p0 = [max(y)-min(y), x(peakIdx), range(x)/4, min(y)];  % [amp, mu, sigma, offset]
    lb = [0, min(x), 0.1, -Inf];
    ub = [Inf, max(x), range(x)*2, Inf];

    opts = optimoptions('lsqcurvefit', 'Display', 'off');
    pFit = lsqcurvefit(gaussEq, p0, x, y, lb, ub, opts);

    fitFun = @(xq) gaussEq(pFit, xq);
    fwhm   = 2.3548 * abs(pFit(3));  % 2*sqrt(2*ln2)*sigma
end
