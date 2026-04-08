function [signalTrace,neuropCorrPars]=estimateNeuropil(cellRoiTrace,neuropRoiTrace,opt)
%
% estimateNeuropil estimates the correction factor r to be applied to the
% neuropil subtraction assuming that being C the proper signal of the cell,
% N the contamination from neuropil, and S the measured signal of the cell
% then C = S - rN
%
% INPUTS:
% cellRoiTrace: (S) time traces of fluorescence in the cell ROI
% neuropRoiTrace: (N) time traces of fluorescence surrounding the cell ROI
% opt.numN: is the number of windows in whoch compute the discretized
% values of Neuropil (default is 20)
% opt.minNp: is the lowest percentile of N from which
% compute the windows (default is 10)
% opt.maxNp: is the highest percentile of N from which
% compute the windows (default is 90)
% opt.pCell: is the lowest percentile of S to be applied to each discrete
% window of N (default is 5)
% opt.noNeg: set to 1 if negative correlations between C and N are
% prohibited (will be set to 0), otherwise set to 0
%
% OUTPUTS:
% cellRoiTrace: (C) proper signal of the cell computed by C=S-rN
% neuropCorrPars.fitNeuro: are the discretized values of neuropil used to
% compute the correction factor r
% neuropCorrPars.corrFactor: nCells x 2 matrix, second column (r) is the
% correction factor. It is obtained by a linear fit of fitNeuro and lowCell
% neuropCorrPars.lowCell: is the lowest percentile of S for each discrete
% window of N (default is 5)
%
% 2015.06.09 Mario Dipoppa - Created
% 2015.11.08 Sylvia - added opt.noNeg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<3 || ~isfield(opt, 'numN')
    numN=20; 
else
    numN=opt.numN;
end
if nargin<3 || ~isfield(opt, 'minNp')
    minNp=10;
else
    minNp=opt.minNp;
end
if nargin<3 || ~isfield(opt, 'maxNp')
    maxNp=90;
else
    maxNp=opt.maxNp;
end
if nargin<3 || ~isfield(opt, 'pN')
    pCell=5;
else
    pCell=opt.pCell;
end
if nargin<3 || ~isfield(opt, 'noNeg')
    noNeg = 1;
else
    noNeg = opt.noNeg;
end
if nargin<3 || ~isfield(opt, 'constrainedFit')
    constrainedFit = 0;
else
    constrainedFit = opt.constrainedFit;
end
if nargin<3 || ~isfield(opt, 'window')
    window = Inf;
else
    window = opt.window;
end
if nargin<3 || ~isfield(opt, 'stepSize')
    stepSize = round(window / 20);
else
    stepSize = opt.stepSize;
end
if nargin<3 || ~isfield(opt, 'verbose')
    verbose = 1;
else
    verbose = opt.verbose;
end

[nCells,nT]=size(neuropRoiTrace);

fitNeuro=nan(nCells,numN);
lowCell=nan(nCells,numN);
corrFactor=nan(nCells,2);
signalTrace=nan(nCells,nT);
if window < nT
    numPad = stepSize - mod(nT - window, stepSize);
    cols = ceil((nT - window) / stepSize) + 1;
    
    lowCell = repmat(lowCell, 1, 1, cols);
    fitNeuro = repmat(fitNeuro, 1, 1, cols);
    corrFactor = repmat(corrFactor, 1, 1, nT);
end
if verbose ~= 1
    options = optimoptions(@fmincon, 'Display', 'none');
else
    options = optimoptions(@fmincon, 'Display', 'final');
end

for iCell=1:nCells
    if verbose == 1
        if mod(iCell,5)>0
            fprintf([num2str(iCell) ' '])
        else
            fprintf([num2str(iCell) '\n'])
        end
    end
    
    traceN = neuropRoiTrace(iCell,:)';
    traceC = cellRoiTrace(iCell,:)';
    % in case, neuropil and neural traces are fit in moving windows,
    % extract traces of these windows in matrix [windowLength x numWindows]
    if window < length(traceN)
        traceN = [traceN; NaN(numPad,1)];
        traceC = [traceC; NaN(numPad,1)];
        ind = bsxfun(@plus, (1:window)', (0:cols-1).*stepSize);
%         ind = spdiags(ones(cols, window), -window+1:0, length(traceN), cols);
        traceN = traceN(ind);
        traceC = traceC(ind);
    end
    % get low and high percentile of neuropil trace in each window
    prcNeurop = prctile(traceN, [minNp maxNp], 1);
    % for each window, divide neuropil values into numN groups
    binSize = (prcNeurop(2,:) - prcNeurop(1,:)) ./ numN;
    discrNeuro=round((numN-1) .* bsxfun(@rdivide, ...
        bsxfun(@minus, traceN, prcNeurop(1,:)+0.5*binSize), ...
        binSize .* (numN-1))) + 1;
    
    %discrNeuro are the discretized values of neuropil between minN and
    % maxN, with numN elements
    
    % for each window and each neuropil value group, find the according
    % value of neural trace (only regard low values)
    for iN=1:numN
        tmp = NaN(size(traceC));
        tmp(discrNeuro == iN) = traceC(discrNeuro == iN);
        lowCell(iCell, iN, :) = prctile(tmp, pCell, 1);
    end
    
    % values of neuropil that will be fit to neural trace values
    fitNeuro(iCell,:,:)= permute(bsxfun(@plus, prcNeurop(1,:), ...
        bsxfun(@times, (1:numN)'./ numN, diff(prcNeurop, 1, 1))), [3 1 2]);
    
    if constrainedFit == 1
        if window < numel(traceN)
            roi = squeeze(lowCell(iCell,:,:));
            npil = squeeze(fitNeuro(iCell,:,:));
            % minimize: roi - (alpha*npil + beta) under constraints 
            % 0 <= alpha <= 2 and 0 <= beta
            % starting condition: alpha = 0, beta = mean(roi)
%             corrFTmp = fmincon(@(pars) myRegress(pars, npil, roi), ...
%                 [mean(roi,1); zeros(1, size(roi,2))], [], [], [], [], ...
%                 zeros(2, size(roi,2)), repmat([Inf 2]', 1, size(roi,2)));
            corrFTmp = NaN(2, cols);
            for win = 1:size(roi,2)
                corrFTmp(:,win) = fmincon(@(pars) myRegress(pars, ...
                    npil(:,win), roi(:,win)), [mean(roi(:,win)); 0], ...
                    [], [], [], [], [0; 0], [Inf; 2], [], options);
            end
        else
            roi = lowCell(iCell,:);
            npil = fitNeuro(iCell,:);
            corrFactor(iCell,:) = fmincon(@(pars) myRegress(pars, npil, roi), ...
                [mean(roi); 0], [], [], [], [], [0; 0], [Inf; 2], [], options);
        end
    else

        if window < numel(traceN)
            corrFTmp = NaN(2,cols);
            for win = 1:size(lowCell,3)
                corrFTmp(:,win) = robustfit(fitNeuro(iCell,:,win), ...
                    lowCell(iCell,:,win));
            end
        else
           corrFactor(iCell,:) = robustfit(fitNeuro(iCell,:), lowCell(iCell,:), [], [], ops);
        end
        % if negative scaling factors are to be avoided, set all negative
        % scalings to zero and set the intercept to the mean of neural
        % response data points
        if noNeg == 1
            if window < numel(traceN)
                ind = corrFTmp(2,:) < 0;
                corrFTmp(2,ind) = 0;
                corrFTmp(1,ind) = mean(lowCell(iCell,:,ind),2);
            else
                ind = corrFactor(2,:) < 0;
                corrFactor(2,ind) = 0;
                corrFactor(iCell,1,ind) = mean(lowCell(iCell,:,ind),2);
            end
        end
        %fit between discretized Neuropil and lowest percentile of signal in
        %the cell ROI
    end
    
    if size(lowCell,3) > 1
        t = [1 round(window/2)+(0:size(lowCell,3)-1).*stepSize nT];
        alphas = interp1(t, corrFTmp(2,[1 1:end end]), 1:nT, 'pchip');
        betas = interp1(t, corrFTmp(1,[1 1:end end]), 1:nT, 'pchip');
        corrFactor(iCell,:,:) = reshape([betas; alphas], 1, 2, nT);
    end
    
    signalTrace(iCell,:)=cellRoiTrace(iCell,:) - ...
        reshape(corrFactor(iCell,2,:),1,[]) .* neuropRoiTrace(iCell,:);
    
end

neuropCorrPars.fitNeuro=fitNeuro;
neuropCorrPars.corrFactor=corrFactor;
neuropCorrPars.lowCell=lowCell;


function error = myRegress(pars, x, y)
error = y - (bsxfun(@plus, bsxfun(@times, pars(2,:), x), pars(1,:)));
error = sqrt(sum(error(:).^2));