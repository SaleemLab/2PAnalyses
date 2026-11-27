function F = sk_kernel_features(ks)
% One row of shape features per cell for a kernel set from sk_kernelset. All shape
% features use ks.G (fractional change in activity, comparable across fits).
%
% Shape
%   depth, peakG, troughG   max - min, max and min of G
%   mainPeak, mainTrough    position (cm) of the largest peak / deepest trough
%   fwhm, riseW, fallW      width of the main peak at half height (re median), and its
%                           left (rising) and right (falling) parts, cm
%   asym                    (riseW - fallW) / (riseW + fallW); > 0 slow rise, sharp fall
%   nPeaks, nTroughs        peaks / troughs with prominence > 50% of the range
%   peaks                   all peak positions (text, cm)
%   rough                   mean squared 2nd difference / variance (0 smooth, ~6 white noise)
%   onset, endEff           mean G over 0-10 cm and 190-200 cm minus mean over 10-190 cm
%   pow24, pow40            fraction of variance (10-200 cm) in a sinusoid of period 24 / 40 cm
% Landmarks (textures pass the animal at 40, 80, 120, 160 cm)
%   pkTexDist               signed distance of the main peak to the nearest texture (cm,
%                           < 0 before it), NaN if > 20 cm away
%   nPkPre, nPkPost         peaks 0-8 cm before / after a texture
% Relation to the rest of the cell's activity (GLM1 reference, base trials)
%   corrObs, corrVis, corrPos   correlation of G with the measured map, GLM1 vision+speed
%                               prediction, and GLM1 position prediction
%   lagVis                  shift (cm) maximising corr of G(x) with vision(x + lag), |lag| <= 20;
%                           > 0: the kernel feature comes before the matching visual feature
%   corrSpeedMap            correlation of G with the session's mean speed profile
%   corrSpeedMap10          the same beyond 10 cm (the animal is slowest at trial onset)
%   speedCorr               correlation of activity with running speed (trial x bin samples)
%   visDepth, posShare      range of GLM1 vision prediction / its mean; share of GLM1 track
%                           modulation from position: var(pos) / (var(pos) + var(vision))
%   relMap                  odd/even trial reliability of the measured map
%   snr                     var(K) / mean(SE^2) where SEs exist
x = ks.x(:)'; n = size(ks.G, 1); tex = [40 80 120 160]; mid = x > 10 & x < 190; late = x > 10;
nanv = nan(n, 1); F = table(nanv, nanv, nanv, nanv, nanv, nanv, nanv, nanv, nanv, nanv, nanv, strings(n,1), nanv, nanv, nanv, nanv, nanv, ...
    nanv, nanv, nanv, 'VariableNames', {'depth','peakG','troughG','mainPeak','mainTrough','fwhm','riseW','fallW','asym','nPeaks','nTroughs','peaks', ...
    'rough','onset','endEff','pow24','pow40','pkTexDist','nPkPre','nPkPost'});
pw = @(g, P) 2*abs(sum(g .* exp(-2i*pi*x(late)/P)))^2 / (nnz(late) * sum(g.^2));
for c = 1:n
    g = ks.G(c,:); rg = max(g) - min(g);
    [mx, mp] = max(g); [mn, tp] = min(g); md = median(g); h = md + (mx - md)/2;
    l = mp; while l > 1 && g(l-1) > h, l = l - 1; end
    r = mp; while r < numel(g) && g(r+1) > h, r = r + 1; end
    F.depth(c) = rg; F.peakG(c) = mx; F.troughG(c) = mn; F.mainPeak(c) = x(mp); F.mainTrough(c) = x(tp);
    F.riseW(c) = 2*(mp - l) + 1; F.fallW(c) = 2*(r - mp) + 1; F.fwhm(c) = 2*(r - l + 1);
    F.asym(c) = (F.riseW(c) - F.fallW(c)) / (F.riseW(c) + F.fallW(c));
    if rg > 0
        [~, lp] = findpeaks([min(g) g min(g)], 'MinPeakProminence', 0.5*rg); pk = x(lp - 1);
        [~, lt] = findpeaks(-[max(g) g max(g)], 'MinPeakProminence', 0.5*rg);
    else, pk = []; lt = []; end
    F.nPeaks(c) = numel(pk); F.nTroughs(c) = numel(lt); F.peaks(c) = strjoin(string(pk), ' ');
    d = pk(:) - tex; [~, j] = min(abs(d), [], 2); d = d(sub2ind(size(d), (1:numel(pk))', j));
    F.nPkPre(c) = nnz(d >= -8 & d < 0); F.nPkPost(c) = nnz(d > 0 & d <= 8);
    dm = x(mp) - tex; [~, j] = min(abs(dm)); F.pkTexDist(c) = dm(j); if abs(dm(j)) > 20, F.pkTexDist(c) = NaN; end
    F.rough(c) = mean(diff(g, 2).^2) / var(g);
    F.onset(c) = mean(g(x < 10)) - mean(g(mid)); F.endEff(c) = mean(g(x > 190)) - mean(g(mid));
    gl = g(late) - mean(g(late)); F.pow24(c) = pw(gl, 24); F.pow40(c) = pw(gl, 40);
end

% relation to the other parts of the cell's activity
cr = @(A, B) arrayfun(@(i) corr(A(i,:)', B(i,:)', 'rows', 'complete'), (1:n)');
V = ks.comp.mapVS; Pp = ks.comp.mapPos;
F.corrObs = cr(ks.G, ks.comp.mapObs); F.corrVis = cr(ks.G, V); F.corrPos = cr(ks.G, Pp);
F.corrSpeedMap = cr(ks.G, ks.comp.speedMap);
F.corrSpeedMap10 = cr(ks.G(:, late), ks.comp.speedMap(:, late));    % without the onset bins
lags = -10:10; F.lagVis = nan(n,1);
for c = 1:n
    rr = nan(size(lags));
    for k = 1:numel(lags)
        s = lags(k); i = max(1, 1 - s):min(numel(x), numel(x) - s);
        rr(k) = corr(ks.G(c, i)', V(c, i + s)');
    end
    [~, k] = max(rr); F.lagVis(c) = 2*lags(k);
end
F.speedCorr = ks.comp.speedCorr; F.relMap = ks.comp.relMap;
F.visDepth = range(V, 2) ./ mean(V, 2); F.posShare = var(Pp, 0, 2) ./ (var(Pp, 0, 2) + var(V, 0, 2));
F.snr = var(ks.K, 0, 2) ./ mean(ks.SE.^2, 2);
F = [ks.cells(:, {'CellIdx','mouseID','sessionid','cellNum'}), table(ks.sig, ks.sigLR, 'VariableNames', {'sig','sigLR'}), F];
end
