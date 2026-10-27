function adj_p = fdr_bh(p)
% FDR_BH  Benjamini-Hochberg FDR-adjusted p-values.
%   adj_p = fdr_bh(p) returns adjusted p-values the same size as p.
%   Compare adj_p < alpha (e.g. 0.05) instead of raw p < alpha.
%
%   NaNs in p are preserved as NaN in the output and excluded from the
%   correction (use this for ROIs you've already flagged unreliable /
%   excluded from testing).

sz = size(p);
p  = p(:);

validIdx = ~isnan(p);
pValid   = p(validIdx);
n        = length(pValid);

adj_p = nan(size(p));

if n == 0
    adj_p = reshape(adj_p, sz);
    return
end

[sortedP, sortOrder] = sort(pValid);
ranks = (1:n)';

% BH adjusted p-values, with the standard monotonicity enforcement
% (cumulative minimum from the largest rank down)
adjSorted = sortedP .* n ./ ranks;
adjSorted = cummin(adjSorted(end:-1:1));
adjSorted = adjSorted(end:-1:1);
adjSorted = min(adjSorted, 1); % clip at 1

adjValid(sortOrder) = adjSorted;
adj_p(validIdx) = adjValid;

adj_p = reshape(adj_p, sz);
end
