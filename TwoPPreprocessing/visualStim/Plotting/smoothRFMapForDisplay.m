function smoothedMap = smoothRFMapForDisplay(rawMap, scaleFactor, sigma)
% smoothRFMapForDisplay: Upsamples then Gaussian-blurs a spatial map for a
% nicer continuous-looking display, matching the imresize+imgaussfilt
% pattern used for the spatial maps in sparseNoiseAnalysis.m. This looks
% much smoother than applying imgaussfilt directly to a small (e.g. 5x3)
% raw grid, since there's more resolution for the blur to work with.
%
% INPUTS:
%   rawMap      - nEl x nAz raw spatial map (e.g. meanGridResponse, or
%                 the output of computeSVDMap)
%   scaleFactor - (optional) upsampling factor passed to imresize. Default: 10.
%   sigma       - (optional) Gaussian blur sigma (in upsampled pixels)
%                 passed to imgaussfilt. Default: 3.
%
% OUTPUT:
%   smoothedMap - upsampled, blurred map. Pass straight into imagesc with
%                 the ORIGINAL uAz/uEl_plot axis limits — imagesc stretches
%                 the image to fit the given X/Y range regardless of the
%                 map's pixel resolution, so no coordinate adjustment needed.
%
% USAGE:
%   smoothedMap = smoothRFMapForDisplay(boutonData.meanGridResponse);
%   imagesc(uAz, uEl_plot, smoothedMap);

    if nargin < 2 || isempty(scaleFactor), scaleFactor = 10; end
    if nargin < 3 || isempty(sigma), sigma = 3; end

    smoothedMap = imresize(rawMap, scaleFactor);
    smoothedMap = imgaussfilt(smoothedMap, sigma);
end
