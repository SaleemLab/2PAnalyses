function [Y, phi] = sk_simulate_null(iSess, lagSamples, seed)
% Simulated activity for one session with NO position term, for checking which kernel
% shapes the fitting pipeline produces on its own.
%   rate  = GLM1's vision + speed prediction (/Spk/GLMprediction_noSpace1), per sample,
%           shifted by LAGSAMPLES (> 0: the simulated cell responds later than the
%           model's visual timing, < 0: earlier; 0: the model's own timing)
%   Y     = phi * Poisson(rate / phi): overdispersed counts with each cell's own
%           dispersion phi (Pearson, from the recorded activity against GLM1's full
%           prediction), so the noise level matches the data
% Returns Y with the same size as the recorded Y in <session>_samples.mat.
if nargin < 2, lagSamples = 0; end
if nargin < 3, seed = 1; end
C = sk_config(); rng(seed);
mouse = C.sessions{iSess,1}; date = C.sessions{iSess,2};
d = dir(fullfile(C.rawRoot, mouse, num2str(date), 'SpkLin_EXP_*.mat')); f = fullfile(d(1).folder, d(1).name);
mu = double(h5read(f, '/Spk/GLMprediction_noSpace1'));
full = double(h5read(f, '/Spk/GLMprediction1'));
L = load(fullfile(C.outRoot, sprintf('%s_%d_samples.mat', mouse, date)), 'Y', 'fit');
assert(isequal(size(mu), size(L.Y)), 'prediction and activity sizes differ');
in = L.fit.inGLM;
% dispersion per cell, from the recorded activity against GLM1's full prediction
phi = zeros(1, size(mu, 2));
for c = 1:size(mu, 2)
    m = max(full(in, c), 1e-3 * mean(full(in, c))); y = double(L.Y(in, c));
    phi(c) = max(sum((y - m).^2 ./ m) / (nnz(in) - 1), 1);
end
% rate: model timing shifted by the lag, floored at 1% of the cell's mean. The prediction
% is NaN between trials, so where the shift brings a gap into a predicted sample (the
% first or last LAGSAMPLES of each stretch) the unshifted value is kept.
mu0 = mu;
if lagSamples > 0, mu = [repmat(mu(1,:), lagSamples, 1); mu(1:end - lagSamples, :)];
elseif lagSamples < 0, mu = [mu(1 - lagSamples:end, :); repmat(mu(end,:), -lagSamples, 1)]; end
gap = isnan(mu) & ~isnan(mu0); mu(gap) = mu0(gap);
mu = max(mu, 0.01 * mean(mu(in,:), 1));
Y = single(phi .* poissrnd(mu ./ phi));
end
