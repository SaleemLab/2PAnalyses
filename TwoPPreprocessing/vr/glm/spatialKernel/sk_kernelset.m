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
% 'refit1_gauss', 'refit2_en_vsip', 'refit2_ridge_vsip', 'refit2_en_vsp', 'refit2_ridge_vsp'.
C = sk_config(); od = C.outRoot;
T = readtable(C.tablePath);
ks.name = source; ks.x = C.binCentres; SE = [];
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
        error('unknown kernel set %s', source);
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
