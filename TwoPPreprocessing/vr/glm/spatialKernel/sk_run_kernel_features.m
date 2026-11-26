function S = sk_run_kernel_features(sets)
% Feature tables for several kernel sets, side by side. Writes
% <outRoot>/kernelFeatures/<set>.csv and kernel_features_summary.txt (medians over
% significant cells), so a new fit is compared by adding its name to SETS.
if nargin < 1, sets = {'glm1', 'refit1_vsp', 'refit1_vsip', 'refit2_en_vsip'}; end
C = sk_config(); od = fullfile(C.outRoot, 'kernelFeatures'); if ~exist(od, 'dir'), mkdir(od); end
S = struct();
for i = 1:numel(sets)
    ks = sk_kernelset(sets{i}); F = sk_kernel_features(ks);
    writetable(F, fullfile(od, [sets{i} '.csv'])); S.(sets{i}) = F;
end
cols = {'depth','fwhm','asym','nPeaks','nTroughs','rough','onset','endEff','pow24','pow40','pkTexDist', ...
    'corrObs','corrVis','corrPos','lagVis','corrSpeedMap','speedCorr','posShare','relMap'};
fid = fopen(fullfile(od, 'kernel_features_summary.txt'), 'w');
out = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
out('Median over cells with position significant in that fit (n); pkTexDist over peaks within 20 cm\n');
out('%-14s', 'feature'); for i = 1:numel(sets), out('%16s', sets{i}); end; out('\n');
out('%-14s', 'n sig / all'); for i = 1:numel(sets), F = S.(sets{i}); out('%16s', sprintf('%d/%d', nnz(F.sig), height(F))); end; out('\n');
for c = cols
    out('%-14s', c{1}); for i = 1:numel(sets), F = S.(sets{i}); out('%16.3f', median(F.(c{1})(F.sig), 'omitnan')); end; out('\n');
end
out('%-14s', 'pre:post pk'); for i = 1:numel(sets), F = S.(sets{i}); out('%16s', sprintf('%d:%d', sum(F.nPkPre(F.sig)), sum(F.nPkPost(F.sig)))); end; out('\n');
fclose(fid);
end
