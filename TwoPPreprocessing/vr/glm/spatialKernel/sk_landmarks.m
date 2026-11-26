function sk_landmarks()
% Landmark positions: on/off edges of angle bins 1, 2 and 19 of each Tex/BG/End element, base condition
% Every session has the same layout: textures pass the animal at 37-43, 77-83,
% 117-123 and 157-163 cm; the six background elements repeat every 24 cm.
C = sk_config(); od = C.outRoot; x = 1:2:199;
for tag = arrayfun(@(i) sprintf('%s_%d', C.sessions{i,1}, C.sessions{i,2}), 1:size(C.sessions,1), 'UniformOutput', false)
    S = load(fullfile(od, [tag{1} '_samples.mat']), 'fit');
    B = load(fullfile(od, [tag{1} '_binned.mat']), 'condDesign', 'meta');
    bi = S.fit.blockIdx(B.meta.visCols); D = B.condDesign.base;
    fprintf('%s\n', tag{1});
    for k = [3 4 5]
        cols = find(bi == k); nEl = max(1, floor(numel(cols) / 19));
        fprintf(' block %d: %d columns\n', k, numel(cols));
        for e = 1:nEl
            for a = [1 2 19]
                j = cols(min((e-1)*19 + a, numel(cols)));
                d = diff([0; D(:, j) > 0.5; 0]); st = x(d(1:end-1) == 1); en = x(d(2:end) == -1);
                fprintf('  el %d angle %2d: %s\n', e, a, strjoin(arrayfun(@(s,t) sprintf('%d-%d', s, t), st, en, 'UniformOutput', false), ', '));
            end
        end
    end
end
end
