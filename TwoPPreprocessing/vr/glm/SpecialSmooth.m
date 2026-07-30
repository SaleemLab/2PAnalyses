function output = SpecialSmooth(input, win, nGrid, Fcircular)
%output = SpecialSmooth(input, win, nGrid, Fcircular) 
%SpecialSmooth performs nD gaussian smoothing.
%INPUTS:
%- input: nD array to smooth
%- win: 1 x n vector with smoothing window along each dimension, in inverse
%relative units.
%- nGrid: number of bins along each dimension
%- Fcircular: 1 x n vector of boolean. True if the corresponding dimension
%should be smoothed circularly.
%
%OUTPUT:
%- output: smoothed nD array.

sz0 = size(input);
Ndim = ndims(input);
if nargin < 4
    Fcircular = false(1,Ndim);
end

if numel(win) > Ndim
    win = win(1:Ndim);
end
if numel(nGrid) > Ndim
    nGrid = nGrid(1:Ndim);
end
if numel(Fcircular) > Ndim
    Fcircular = Fcircular(1:Ndim);
end

output = input;

for k = 1:Ndim
    if Fcircular(k)
        vrep = ones(1,Ndim);
        vrep(k) = 3;
        input = repmat(input, vrep);
    end
end

if sum(win > 0) > 0
    input0 = input;
    input(isnan(input)) = 0;
    
    output = specialsmooth(input,win,nGrid);
    
    if sum(isnan(input0(:)))>0
        input = double(~isnan(input0));
        flat = specialsmooth(input,win,nGrid);
    else
        input = ones(size(input0));
        flat = specialsmooth(input,win,nGrid); %convn of a uniform array is
        %     extremely slow in Matlab. Must have to do with the FFT transform.
    end
    output = output./flat;
    
    output(isnan(input0)) = NaN;
    
    Fkeep = true(size(output));
    for k = 1:Ndim
        if Fcircular(k)
            Fkeepk = false(size(output,k),1);
            Fkeepk((size(output,k)/3+1):(2*size(output,k)/3)) = true;
            vperm = 1:Ndim;
            vperm([1 k]) = vperm([k 1]);
            Fkeepk = permute(Fkeepk,vperm);
            sz = size(output);
            sz(k) = 1;
            Fkeep = Fkeep .* repmat(Fkeepk, sz);
        end
    end
    output = reshape(output(Fkeep>0),sz0);
end

end

function output = specialsmooth(input,win,nGrid)
Ndim = ndims(input);
npts = 5;
win(isnan(win) | isinf(win)) = 0;
Smoother = ones(2*(npts * round(max(win, 1./nGrid) .* nGrid)) + 1);
for k = 1:Ndim
%     x = (-nGrid(k):nGrid(k))/nGrid(k);
    x = (-(npts * max(win(k), 1./nGrid(k)) * nGrid(k)):(npts * max(win(k), 1./nGrid(k)) * nGrid(k)));
    if win(k) ~= 0 && ~isnan(win(k)) && ~isinf(win(k))
        Smoother_1D = exp(-x.^2/(win(k) * nGrid(k))^2/2);
    else
        Smoother_1D = double(x == 0);
    end
    Smoother_1D = Smoother_1D(:);
    vperm = 1:Ndim;
    vperm([1 k]) = vperm([k 1]);
    Smoother_1D = permute(Smoother_1D,vperm);
    sz = size(Smoother);
    sz(k) = 1;
    Smoother = Smoother .* repmat(Smoother_1D, sz);
end
Smoother = Smoother / sum(Smoother(:));
output = convn(input,Smoother,'same');

% Ndim = ndims(input);
% Smoother = ones(2*nGrid+1);
% for k = 1:Ndim
%     x = (-nGrid(k):nGrid(k))/nGrid(k);
%     if win(k) ~= 0 && ~isnan(win(k)) && ~isinf(win(k))
%         Smoother_1D = exp(-x.^2/win(k)^2/2);
%     else
%         Smoother_1D = double(x == 0);
%     end
%     Smoother_1D = Smoother_1D(:);
%     vperm = 1:Ndim;
%     vperm([1 k]) = vperm([k 1]);
%     Smoother_1D = permute(Smoother_1D,vperm);
%     sz = size(Smoother);
%     sz(k) = 1;
%     Smoother = Smoother .* repmat(Smoother_1D, sz);
% end
% Smoother = Smoother / sum(Smoother(:));
% output = convn(input,Smoother,'same');

end

