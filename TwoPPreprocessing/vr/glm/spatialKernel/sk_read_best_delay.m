function d = sk_read_best_delay(f)
% Per-cell best visual delay (ms) chosen by the original GLM1 fit,
% GLMs{1}.Perf.bestVisdelay, cached next to the exported data.
C = sk_config();
[~, n] = fileparts(f);
cache = fullfile(C.outRoot, [n '_bestVisdelay.mat']);
if exist(cache, 'file'), S = load(cache); d = S.d; return; end
S = load(f, 'GLMs');
d = double(S.GLMs{1}.Perf.bestVisdelay(:));
save(cache, 'd');
end
