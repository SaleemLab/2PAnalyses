function ks = sk_kernelset(source)
% Load a set of position kernels in one standard form, so that every shape analysis
% runs unchanged on GLM1, the refits, and later fits (add a case below for a new model).
%
%   ks.name    source name;  ks.link 'identity' (additive kernel) or 'log' (gain kernel)
%   ks.x       bin centres (cm), 1 x 100
%   ks.K       kernels as fitted, nCells x 100 (gaps filled linearly)
%   ks.SE      kernel standard errors where the fit gives them, else NaN
%   ks.G       kernels on a common scale: fractional change in activity relative to the
%              track mean. log link: exp(K - mean K) - 1. identity link: (K - mean K) /
%              mean GLM1 full prediction on the track (base trials).
%   ks.cells   table of cell information (CellIdx, CellName, mouseID, sessionid, date, cellNum, spatialcells)
%   ks.sig     position significant for this fit (fold t-test, p < 0.05 and LLHrel > 1;
%              GLM1: spatialcells)
%   ks.sigLR   dispersion-corrected likelihood-ratio test where available, else NaN
%   ks.comp    rows of sk_cell_components for these cells
%
% Sources: 'glm1' (GLM1 spatial cells), 'glm1_all', 'refit1_vsp', 'refit1_vsip',
% 'refit1_gauss', 'refit2_en_vsip', 'refit2_ridge_vsip', 'refit2_en_vsp', 'refit2_ridge_vsp',
% and simulated cells without a position term (sk_run_sim_refit): 'sim0_vsip',
% 'simlag6_vsip', 'simlead6_vsp', ... For simulated sets, ks.comp is still the recorded
% cell's (its mapVS is the simulation's generator, its mapObs is the real data).
% Modifiers: append ':all' to treat every cell as significant and/or ':s33' to keep one
% session, e.g. 'refit1_vsip:s33:all'.
C = sk_config(); od = C.outRoot;
T = readtable(C.tablePath);
ks.name = source; ks.x = C.binCentres; SE = [];
parts = strsplit(source, ':'); source = parts{1}; mods = parts(2:end);
switch source
    case {'glm1', 'glm1_all'}
        g = 'D:\DATA\Sonali\GLM\ForAS\V1_GLMData.mat';
        K = h5read(g, '/Kernels/PosKernel'); SE = h5read(g, '/Kernels/PosKernelSE');
        idx = T.CellIdx; if strcmp(source, 'glm1'), idx = T.CellIdx(T.spatialcells == 1); end
        % GLM1 kernels are in activity per second; the binned maps are per 60 Hz sample
        K = K(idx,:) / 60; SE = SE(idx,:) / 60; ks.link = 'identity';
        sig = T.spatialcells(idx) == 1; sigLR = nan(numel(idx), 1);
    case {'refit1_vsp', 'refit1_vsip', 'refit1_gauss'}
        R = readtable(fullfile(od, 'refit_spatialcells.csv')); Kf = load(fullfile(od, 'refit_kernels.mat'));
        assert(isequal(Kf.CellIdx, R.CellIdx)); idx = R.CellIdx;
        switch source
            case 'refit1_vsp',   K = Kf.pois_VSP;  ks.link = 'log';      sig = R.p_pois_pos < 0.05 & R.LLHrel_pois > 1;
            case 'refit1_vsip',  K = Kf.pois_VSIP; ks.link = 'log';      sig = R.p_pois_posI < 0.05 & R.LLHrelI_pois > 1;
            case 'refit1_gauss', K = Kf.gauss_VSP; ks.link = 'identity'; sig = R.p_gauss_pos < 0.05 & R.LLHrel_gauss > 1;
        end
        sigLR = nan(numel(idx), 1);
    case {'refit2_en_vsip', 'refit2_ridge_vsip', 'refit2_en_vsp', 'refit2_ridge_vsp'}
        p = regexp(source, 'refit2_(\w+)_(vsi?p)', 'tokens', 'once'); pen = p{1}; mdl = upper(p{2});
        sfx = ''; if strcmp(mdl, 'VSIP'), sfx = 'I'; end
        K = []; idx = []; sig = []; sigLR = []; ks.link = 'log';
        for f = dir(fullfile(od, 'refit2', 'refit2_M*.mat'))'
            E = load(fullfile(f.folder, f.name)); E = E.Ex;
            N = readtable(fullfile(f.folder, strrep(f.name, '.mat', '.csv')));
            assert(isequal(N.cellNum, E.cellNum));
            [tf, r] = ismember([repmat(E.sessionid, numel(E.cellNum), 1) E.cellNum], [T.sessionid T.cellNum], 'rows'); assert(all(tf));
            K = [K; E.kernels.([pen '_' mdl])]; idx = [idx; T.CellIdx(r)]; %#ok<AGROW>
            sig = [sig; N.(['p_' pen '_pos' sfx]) < 0.05 & N.(['LLHrel' sfx '_' pen]) > 1]; %#ok<AGROW>
            sigLR = [sigLR; N.(['pLRq_' pen '_pos' sfx]) < 0.05]; %#ok<AGROW>
        end
    otherwise
        tk = regexp(source, '^sim(0|lag\d+|lead\d+)_(vsp|vsip|gauss)$', 'tokens', 'once');
        if isempty(tk), error('unknown kernel set %s', source); end
        K = []; idx = []; sig = []; sigLR = [];
        for f = dir(fullfile(od, 'simrefit', sprintf('sim%s_M*.csv', tk{1})))'
            N = readtable(fullfile(f.folder, f.name)); E = load(strrep(fullfile(f.folder, f.name), '.csv', '.mat')); E = E.Ex;
            [tf, r] = ismember([N.sessionid N.cellNum], [T.sessionid T.cellNum], 'rows'); assert(all(tf));
            switch tk{2}
                case 'vsp',   k = E.kernels.pois_VSP;  s = N.p_pois_pos < 0.05 & N.LLHrel_pois > 1;
                case 'vsip',  k = E.kernels.pois_VSIP; s = N.p_pois_posI < 0.05 & N.LLHrelI_pois > 1;
                case 'gauss', k = E.kernels.gauss_VSP; s = N.p_gauss_pos < 0.05 & N.LLHrel_gauss > 1;
            end
            K = [K; k]; idx = [idx; T.CellIdx(r)]; sig = [sig; s]; sigLR = [sigLR; nan(numel(r), 1)]; %#ok<AGROW>
        end
        if isempty(K), error('no simulated refits for %s yet', source); end
        if strcmp(tk{2}, 'gauss'), ks.link = 'identity'; else, ks.link = 'log'; end
end
sig = double(sig); sigLR = double(sigLR);
% modifiers: ':all' treats every cell as selected, ':sNN' keeps session NN only
for md = mods
    if strcmp(md{1}, 'all'), sig = true(size(sig));
    elseif md{1}(1) == 's'
        sid = str2double(md{1}(2:end)); [~, r] = ismember(idx, T.CellIdx); keep = T.sessionid(r) == sid;
        K = K(keep,:); idx = idx(keep); sig = sig(keep); sigLR = sigLR(keep); if ~isempty(SE), SE = SE(keep,:); end
    else, error('unknown modifier %s', md{1}); end
end
ks.K = fillmissing(double(K), 'linear', 2, 'EndValues', 'nearest');
if isempty(SE), ks.SE = nan(size(ks.K)); else, ks.SE = double(SE); end
[~, r] = ismember(idx, T.CellIdx);
ks.cells = T(r, {'CellIdx', 'CellName', 'mouseID', 'sessionid', 'date', 'cellNum', 'spatialcells'});
ks.sig = logical(sig(:)); ks.sigLR = double(sigLR(:));
P = sk_cell_components(); [~, pr] = ismember(idx, P.CellIdx);
for fn = fieldnames(P)'
    v = P.(fn{1}); if size(v, 1) == numel(P.CellIdx), ks.comp.(fn{1}) = double(v(pr,:)); end
end
ks.comp.x = P.x;
kc = ks.K - mean(ks.K, 2);
if strcmp(ks.link, 'log'), ks.G = exp(kc) - 1;
else, ks.G = kc ./ mean(ks.comp.mapFull, 2, 'omitnan'); end
end
