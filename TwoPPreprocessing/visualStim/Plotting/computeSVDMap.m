function svdMap = computeSVDMap(tempStack, timeVector, respWin)
% computeSVDMap: Extracts a denoised spatial tuning map from a bouton's
% per-position temporal responses via SVD, analogous to the
% getResponseSVD/makeRFmap approach in sparseNoiseAnalysis.m — instead of
% just averaging each position's trace within the response window, this
% finds the dominant spatiotemporal component across ALL positions at once
% and uses its spatial loadings as the map. This is more robust to noisy
% individual-position traces than simple windowed averaging.
%
% INPUTS:
%   tempStack   - nTpts x nEl x nAz, a bouton's meanTemporalResponse
%   timeVector  - nTpts x 1, shared time vector
%   respWin     - (optional) [start end] response window (seconds) to
%                 restrict the SVD to; if omitted, uses the full trace
%
% OUTPUT:
%   svdMap      - nEl x nAz spatial map (first SVD spatial component,
%                 sign-corrected so it points in the direction of a
%                 positive average response)
%
% USAGE:
%   svdMap = computeSVDMap(boutonData.meanTemporalResponse, timeVector, respWin);
%   imagesc(uAz, uEl_plot, svdMap);

    [nTpts, nEl, nAz] = size(tempStack);

    dataMat = reshape(tempStack, nTpts, nEl * nAz);  % nTpts x nPositions

    if nargin >= 3 && ~isempty(respWin)
        respIdx = timeVector >= respWin(1) & timeVector <= respWin(2);
        dataMat = dataMat(respIdx, :);
    end

    % SVD can't handle NaNs — zero-fill any missing samples
    dataMat(isnan(dataMat)) = 0;

    if all(dataMat(:) == 0)
        svdMap = zeros(nEl, nAz);
        return;
    end

    [~, ~, V] = svd(dataMat, 'econ');
    spatialWeights = V(:, 1);

    % sign-correct: SVD components are sign-ambiguous, so flip if the
    % average weight is negative (we want positive = "responsive")
    if mean(spatialWeights) < 0
        spatialWeights = -spatialWeights;
    end

    svdMap = reshape(spatialWeights, nEl, nAz);
end
