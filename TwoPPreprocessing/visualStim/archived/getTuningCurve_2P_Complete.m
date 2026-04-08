function [tuningCurveMean, tuningCurveSEM, varAxis, pVal, pVal2, binnedSigArray, peakVal, model, cvmodel, vals, options] = ...
    getTuningCurve_2P_Complete(sig, varTimes, varValues, options)

%% 1. Initialization (Exact Original Struct)
tuningCurveMean = []; tuningCurveSEM = []; varAxis = []; 
pVal = []; pVal2 = []; peakVal = []; model = []; cvmodel = []; vals = [];

if ~exist('options','var'), options=struct; end
if ~isfield(options,'binSize'),      options.binSize= 0.1;           end
if ~isfield(options,'startTime'),    options.startTime= varTimes(1); end
if ~isfield(options,'stopTime'),     options.stopTime= varTimes(end);end
if ~isfield(options,'varLims'),      options.varLims= [-inf inf];    end
if ~isfield(options,'equalBins'),    options.equalBins= true;        end
if ~isfield(options,'nBins'),        options.nBins= 20;              end
if ~isfield(options,'catBelow'),     options.catBelow= -inf;         end
if ~isfield(options,'getSig'),       options.getSig=true;            end
if ~isfield(options,'nPerms'),       options.nPerms=1000;            end
if ~isfield(options,'fitCustomFun'), options.fitCustomFun=true;      end 
if ~isfield(options,'CVcustomFun'),  options.CVcustomFun=true;       end
if ~isfield(options,'kfold'),        options.kfold=5;                end

%% 2. Timing and Truncation
tIdx = varTimes >= options.startTime & varTimes <= options.stopTime;
v = varValues(tIdx);
s = sig(tIdx);

% remove bins where vars outside of var lims
toRemove = (options.varLims(2) < v | v < options.varLims(1));
v(toRemove) = [];
s(toRemove) = [];

%% 3. Generate Bins (Original Quantile Logic)
if options.equalBins
    nBins = options.nBins;
    v_moving = v(v >= options.catBelow);
    if isempty(v_moving), binedges = options.catBelow; else
        binedges = quantile(v_moving, linspace(0, 1, nBins-1));
    end
    
    [~, binnedVar] = histc(v, [-inf; binedges(:); inf]);
    makeone = find(binedges < options.catBelow);
    binnedVar(ismember(binnedVar, makeone)) = 1;
else
    binnedVar = round(v ./ options.binSize) .* options.binSize;
end

uniqueVarVals = unique(binnedVar);
nVarVals = numel(uniqueVarVals);

for ival = 1:nVarVals
    idx = find(binnedVar == uniqueVarVals(ival));
    vals(ival).idx = idx;
    vals(ival).val = nanmean(v(idx));
    vals(ival).nTimes = numel(idx);
    vals(ival).binnedSpikes = s(idx)'; 
    vals(ival).meanRate = mean(s(idx), 'omitnan');
    vals(ival).semRate = std(s(idx), 0, 'omitnan') / sqrt(vals(ival).nTimes);
end
binnedSigArray = {vals.binnedSpikes};
tuningCurveMean = [vals.meanRate];
tuningCurveSEM = [vals.semRate];

%% 4. Significance (Edward's Random Shuffling)
tuningCurveVar = nanvar(tuningCurveMean);
tuningCurveVarNotZero = nanvar(tuningCurveMean(2:end));

if options.getSig
    permTuneVars = nan*ones(options.nPerms,1);
    permTuneVarsNotZero = nan*ones(options.nPerms,1);
    
    for iperm = 1:options.nPerms
        % --- RESTORED: Edward's randperm logic ---
        sPerm = s(randperm(length(s)));
        
        tMeans = nan(1, nVarVals);
        for ival = 1:nVarVals
            tMeans(ival) = mean(sPerm(vals(ival).idx), 'omitnan');
        end
        
        permTuneVars(iperm) = nanvar(tMeans);
        permTuneVarsNotZero(iperm) = nanvar(tMeans(2:end));
    end
    
    pVal = sum(permTuneVars > tuningCurveVar) * (1/options.nPerms);
    pVal2 = sum(permTuneVarsNotZero > tuningCurveVarNotZero) * (1/options.nPerms);
end

[~, peakIdx] = max(tuningCurveMean);
varAxis = [vals.val];
peakVal = varAxis(peakIdx);

%% 5. Model Fitting (Edward's Global Fit)
if options.fitCustomFun
    % Original Exponential Function
    func = @(params, x) params(1) * exp(params(2) * x) + params(3);
    
    xvals = varAxis(2:end)'; 
    yvals = tuningCurveMean(2:end)';
    
    valid = ~isnan(xvals) & ~isnan(yvals);
    if sum(valid) > 3
        % Robust Initial Guesses: [Amplitude, Growth, Offset]
        params0 = [range(yvals), 0.01, min(yvals)];
        
        % Use Edward's MultiStart approach
        problem = createOptimProblem('lsqcurvefit','x0',params0,'objective',func,...
            'lb',[0 -inf -inf],'ub',[inf inf inf],'xdata',xvals(valid),'ydata',yvals(valid));
        ms = MultiStart('Display', 'off');
        
        try
            [params, ~] = run(ms, problem, 50);
            model.params = params;
            model.fun = func;
        catch
            model.params = [];
        end
    end
end

%% 6. Cross-Validation (Edward's Global-to-Fold Logic)
if options.CVcustomFun && isfield(model, 'params') && ~isempty(model.params)
    k = options.kfold;
    for ik = 1:k
        % Pre-allocate so we don't hit "Index Exceeds Bounds"
        cvmodel(ik).meanTest = nan(1, nVarVals); 
        
        for ival = 1:nVarVals
            v_idx = vals(ival).idx;
            foldSize = floor(length(v_idx)/k);
            if foldSize > 0
                testIdx = v_idx( (ik-1)*foldSize+1 : ik*foldSize );
                cvmodel(ik).meanTest(ival) = mean(s(testIdx), 'omitnan');
            end
        end
        
        % --- THE CRITICAL FIX FOR NEGATIVE R2 ---
        y_obs = cvmodel(ik).meanTest(2:end);  % This is a Row vector [1 x N]
        x_moving = varAxis(2:end);            % This is a Row vector [1 x N]
        
        % Get prediction and ensure it's a Row vector to match y_obs
        y_pred = model.fun(model.params, x_moving); 
        if size(y_pred, 1) > size(y_pred, 2), y_pred = y_pred'; end
        
        % Use element-wise subtraction and ignore NaNs
        residuals = (y_obs - y_pred);
        SSres = sum(residuals.^2, 'omitnan');
        
        total_var = (y_obs - mean(y_obs, 'omitnan'));
        SStot = sum(total_var.^2, 'omitnan');
        
        if SStot > 0 && ~isnan(SSres)
            cvmodel(ik).R2 = 1 - (SSres / SStot);
        else
            cvmodel(ik).R2 = NaN;
        end
    end
    % Average R2 across folds
    model.cv_R2_avg = mean([cvmodel.R2], 'omitnan');
end
end

