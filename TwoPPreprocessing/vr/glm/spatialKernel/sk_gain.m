function [g, e, gainPred] = sk_gain(Db, Vb)
% Multiplicative gain of data over the vision+speed prediction, estimated from
% base-condition means. Db, Vb: [nCells x nBins] (Hz), or [nCells x 1 x nBins].
% The Gaussian (identity-link) GLM can predict negative rates, so V is clipped
% at 0; e (20% of the cell's mean rate) regularises the ratio, and g is capped
% at 10 so near-silent bins cannot dominate.
%   gainPred(V) = (max(V,0) + e) .* g - e     prediction for any V of matching size
d = ndims(Db);
e = 0.2*mean(Db, d, 'omitnan') + 0.01;
g = (Db + e) ./ (max(Vb,0) + e);
g = min(max(g, 0), 10);
gainPred = @(V) (max(V,0) + e) .* g - e;
end
