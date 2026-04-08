function [tuningCurveMean, tuningCurveSEM, varAxis, pVal, pVal2, binnedSpikeArray, peakVal, model, cvmodel, vals, options]...
    = getTuningCurve_continous(spikeTimes, varTimes, varValues, options)
% inputs: spikeTimes (vector of spike times)
% varTimes, varValues (vector of var times and corresponding values)
% options struct
% Edward Horrocks
% Adapted for 2p data 

%% To do:
% smoothing options?

%%
tuningCurveMean = [];
tuningCurveSEM = []; 
varAxis = []; 
pVal = []; 
pVal2 = [];
peakVal = []; 
model = []; 
cvmodel = [];
vals = [];

%% options
if ~exist('options','var'),                 options=struct;                 end

% generic options
if ~isfield(options,'binSize'),             options.binSize= 0.1;           end
if ~isfield(options,'varRoundVal'),         options.varRoundVal= 0.1;       end
if ~isfield(options,'startTime'),           options.startTime= [];          end
if ~isfield(options,'stopTime'),            options.stopTime= [];           end
if ~isfield(options,'varLims'),             options.varLims= [-inf inf];    end
if ~isfield(options,'plot'),                options.plot=false;             end

% equal bins tuning curve options
if ~isfield(options,'equalBins'),           options.equalBins= true;        end
if ~isfield(options,'nBins'),               options.nBins= [20];            end
if ~isfield(options,'catBelow'),            options.catBelow= -inf;         end
if ~isfield(options,'catAbove'),            options.catAbove = inf;         end

% significance testing options
if ~isfield(options,'getSig'),              options.getSig=false;           end
if ~isfield(options,'sigType'),             options.sigType='perm';         end % or 'KW'
if ~isfield(options,'nPerms'),              options.nPerms=100;             end % only variance atm.

% descriptive function options
if ~isfield(options,'fitCustomFun'),        options.fitCustomFun=false;     end 
if ~isfield(options,'customFun'),           options.customFun=[];           end 
if ~isfield(options,'customFunIni'),        options.customFunIni=[];        end 
if ~isfield(options,'customFunStartVal'),   options.customFunStartVal=1;    end
if ~isfield(options,'CVcustomFun'),         options.CVcustomFun=false;      end
if ~isfield(options,'kfold'),               options.kfold=5;                end

%% Deal with timing information
binSize = options.binSize;
varRoundVal = options.varRoundVal;

% find timing over which to compute tuning curve
if isempty(options.startTime)
    varStart = varTimes(1);
    spikeStart = spikeTimes(1);
    startTime = max([spikeStart varStart]);
else
    startTime = options.startTime;
end

if isempty(options.stopTime)
    varEnd = varTimes(end);
    spikesEnd = spikeTimes(end);
    stopTime = min([varEnd, spikesEnd]);
else
    stopTime = options.stopTime;
end

% bin spikes into binSize bins and get midpoints of bins
binedges = startTime:binSize:stopTime;
binMidPoints = binedges+(binSize/2);
binMidPoints(end) = [];

%% Bin spikes, interpolate, round and truncate input var
binnedSpikes = histcounts(spikeTimes, binedges);

% interpolate varValues using new bin midpoints
var_interped = interp1(varTimes,varValues,binMidPoints);

% remove bins where vars outside of var lims
toRemove = (options.varLims(2)<var_interped | var_interped<options.varLims(1));
var_interped(toRemove)=[];
binnedSpikes(toRemove)=[];

%% generate struct for each unique var and compute firing rate for it

if options.equalBins
    nBins = options.nBins;
    binedges = quantile(var_interped,nBins-1);
    if binedges(1) > options.catBelow
        binedges = quantile(var_interped,nBins-2);
        binedges = [options.catBelow, binedges];
    end
    makeone = find(binedges<options.catBelow);
    makemax = find(binedges>options.catAbove);
    [~, binnedVar] = histc(var_interped,[-inf;binedges(:);inf]);
    binnedVar(ismember(binnedVar,makeone)) = 1;
    binnedVar(ismember(binnedVar,makemax)) = nBins;
    
else
    if varRoundVal ~= 0
        var_interped = round(var_interped./varRoundVal).*varRoundVal;
    end
    binnedVar = var_interped;
end


% unique var vals are defined as input space
% pre-allocate struct array for each unique value of input variable
uniqueVarVals = unique(binnedVar);
nVarVals = numel(uniqueVarVals);
vals(nVarVals).idx = [];
vals(nVarVals).nTimes = [];
vals(nVarVals).meanRate = [];
vals(nVarVals).semRate = [];
vals(nVarVals).meanRateHz = [];
vals(nVarVals).semRateHz = [];

% compute mean and sem firing rate for each unique var
for ival = 1:nVarVals
    vals(ival).idx = find(binnedVar==uniqueVarVals(ival));
    vals(ival).val = nanmean(var_interped(vals(ival).idx));
    vals(ival).nTimes = numel(vals(ival).idx);
    vals(ival).binnedSpikes = binnedSpikes(vals(ival).idx)';
    vals(ival).meanRate = mean(binnedSpikes(vals(ival).idx));
    vals(ival).semRate = std(binnedSpikes(vals(ival).idx))./sqrt(vals(ival).nTimes);
    vals(ival).meanRateHz = vals(ival).meanRate.*(1/binSize);
    vals(ival).semRateHz = vals(ival).semRate.*(1/binSize);
end

binnedSpikeArray = {vals.binnedSpikes};

%% tuning curve + error bars
tuningCurveMean = [vals.meanRateHz];
tuningCurveSEM = [vals.semRateHz];
tuningCurveVar = nanvar(tuningCurveMean);
tuningCurveVarNotZero = nanvar(tuningCurveMean(2:end));

[~, peakvalidx] = max(tuningCurveMean);
varAxis = [vals.val];
peakVal = varAxis(peakvalidx);


%% check significance of tuning by comparing to shuffled spikes
pVal = [];
if options.getSig
    if strcmp(options.sigType,'KW')
        
        allSpikes = horzcat(vals.binnedSpikes);
        groupingVals = categorical(repelem([1:numel(vals)], [vals.nTimes]));
        
        pVal = kruskalwallis(allSpikes,groupingVals,'off');
        
    elseif strcmp(options.sigType,'perm')
        
        permTuneVars = nan*ones(options.nPerms,1);
        permTuneVarsNotZero = nan*ones(options.nPerms,1);
        for iperm = 1:options.nPerms
            
            binnedSpikesPerm = binnedSpikes(randperm(length(binnedSpikes)));
            
            % compute mean and sem firing rate for each unique var
            for ival = 1:nVarVals
                vals_shuf(ival).meanRate = mean(binnedSpikesPerm(vals(ival).idx)); % already computed the idx
                vals_shuf(ival).meanRateHz = vals_shuf(ival).meanRate.*(1/binSize);
            end
            
            permTuneVars(iperm) = nanvar([vals_shuf.meanRateHz]);
            tss = [vals_shuf.meanRateHz];
            permTuneVarsNotZero(iperm) = nanvar(tss(2:end));
            
        end
        
        pVal = sum(permTuneVars>tuningCurveVar)*(1/options.nPerms);
        pVal2 = sum(permTuneVarsNotZero>tuningCurveVarNotZero)*(1/options.nPerms);
        
    end
end

%% basic errorbar plot if you want it
% if options.plot
% errorbar(varAxis, tuningCurveMean, tuningCurveSEM, 'LineStyle', 'None',...
%     'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 2);    
% title(num2str(pVal));
% end


%% fitting an custom function to the data.

if options.fitCustomFun
    cvmodel = [];
    func = options.customFun;
    %params0 = options.customFunIni;
    startVal = options.customFunStartVal;

    [maxVal, peakidx] = max(tuningCurveMean(startVal:end));
    params0 = [maxVal, peakVal, 100, -100];
    startVal = options.customFunStartVal;
    k = options.kfold;
    testProp = 1/k;
    ub = [inf inf inf inf];
    lb = [-inf -inf -inf -inf];
%     
%     allx = []; ally = [];
%     for ival = startVal:numel(vals)
%         allx = [allx; repelem(vals(ival).val, numel(vals(ival).binnedSpikes), 1)];
%         ally = [ally; vals(ival).binnedSpikes(:)./binSize];
%     end
    
    
    % fit function onto all data to get good initial conditions
    xvals = varAxis(:); xvals = xvals(startVal:end);
    yvals = tuningCurveMean(:); yvals = yvals(startVal:end);
    fun = @(params)sseexp2sigma(params,xvals,yvals);
    
     problem = createOptimProblem('lsqcurvefit','x0',params0,'objective',func,...
     'lb',lb,'ub',ub,'xdata',xvals,'ydata',yvals);

%     problem = createOptimProblem('fmincon','x0',params0,'objective',fun,...
%     'lb',lb,'ub',ub,'xdata',xvals,'ydata',yvals);
    ms = MultiStart();
    [params,errormulti] = run(ms,problem,500);
   

    hold on
    plot(varAxis(2):0.1:varAxis(end), func(params, varAxis(2):0.1:varAxis(end)), 'k-')
    close
    
    model.params = params;
    model.fun = func;
    model.params
    

    
    %% cross-validated prediction quality
    % get 20% of indexes from each bin
    % get other 80% -> construct mean curve, fit function, test
    if options.CVcustomFun
    % random perm idx and truncate so multiple of k
    for ival = 1:numel(vals)
        cv(ival).idx = 1:numel(vals(ival).idx);
        cv(ival).idx = cv(ival).idx(randperm(numel(cv(ival).idx)));
        cv(ival).idx = cv(ival).idx(1:end-mod(numel(cv(ival).idx),k));
        cv(ival).numTestIdx = floor(numel(cv(ival).idx)*testProp);
        for ik = 1:k
            cv(ival).test_idx(:,ik) =...
                cv(ival).idx(((ik*cv(ival).numTestIdx)-(cv(ival).numTestIdx-1)):(ik*cv(ival).numTestIdx));
            cv(ival).train_idx(:,ik) = ...
                cv(ival).idx(~ismember(cv(ival).idx, cv(ival).test_idx(:,ik)));
        end
    end
    
    for ik = 1:k
        for ival = 1:numel(vals)
            cvmodel(ik).meanTrainRate(ival) =...
                mean(vals(ival).binnedSpikes(round(cv(ival).train_idx(:,ik))))/binSize;
            cvmodel(ik).meanTestRate(ival) =...
                mean(vals(ival).binnedSpikes(round(cv(ival).test_idx(:,ik))))/binSize;
        end
        
        cvmodel(ik).params = lsqcurvefit(func,params,xvals,cvmodel(ik).meanTrainRate(startVal:end)');
        cvmodel(ik).pred = func(cvmodel(ik).params,xvals);
        cvmodel(ik).meanTest = mean(cvmodel(ik).meanTestRate(startVal:end));
        
        cvmodel(ik).SSerr = sum((cvmodel(ik).pred-cvmodel(ik).meanTestRate(startVal:end)').^2);
        cvmodel(ik).SStot = sum((cvmodel(1).meanTest-cvmodel(ik).meanTestRate(startVal:end)').^2);
        cvmodel(ik).R2 = 1-(cvmodel(ik).SSerr/cvmodel(ik).SStot);
        
    end
    end
    
end


end
