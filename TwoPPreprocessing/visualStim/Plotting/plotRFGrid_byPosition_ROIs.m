function plotRFGrid_byPosition_ROIs(...
    sessionFileInfo, RFMappingStimName, prestimWindow, poststimWindow, metric)
% plotRFGrid_byPosition_ROIs  Heatmap of per‑ROI RF metrics across az/el

% defaults
if nargin<3, prestimWindow  = 0.5; end
if nargin<4, poststimWindow = 1; end
if nargin<5, metric         = 'peak'; end
if isscalar(prestimWindow),  prestimWindow  = [-abs(prestimWindow) 0]; end
if isscalar(poststimWindow), poststimWindow = [0 abs(poststimWindow)]; end

% define cmap 

% load PSTH data
iS = find(strcmp(RFMappingStimName,{sessionFileInfo.stimFiles.name}),1);
assert(~isempty(iS),'Stimulus not found.');
response = load(sessionFileInfo.stimFiles(iS).Response,'response');
psth = response.response.psthData;

% compute az/el grid
stimVals   = vertcat(psth.stimValue);
azimuths   = unique(stimVals(:,1));
elevations = unique(stimVals(:,2));
numAz = numel(azimuths);
numEl = numel(elevations);

% sizes & time masks
nPos = numel(psth);
nROI = size(psth(1).alignedResponses,1);
tVec = psth(1).timeVector(:)';
preMask  = tVec>=prestimWindow(1)  & tVec<=prestimWindow(2);
postMask = tVec>=poststimWindow(1) & tVec<=poststimWindow(2);

% compute metric for each ROI×position
allVals = nan(nROI, nPos);
for k = 1:nPos
    AR3d = psth(k).alignedResponses;
    for r = 1:nROI
        tmp = squeeze(AR3d(r,:,:));
        if size(tmp,1)==1
            tmp = tmp';
        end
        mTr = median(tmp,2,'omitnan');
        F0  = median(mTr(preMask),'omitnan'); % change to median? 
        zn  = mTr - F0;
        seg = zn(postMask);
        switch lower(metric)
            case 'mean'
                allVals(r,k) = mean(seg,'omitnan');
            case 'median'
                allVals(r,k) = median(seg,'omitnan');
            case 'peak'
                allVals(r,k) = max(seg,[],'omitnan');
            otherwise
                error('Unknown metric "%s".', metric);
        end
    end
end

% robust color limits (5th–95th percentile)
flat  = allVals(:); flat = flat(~isnan(flat));
clow  = prctile(flat,5);
chigh = prctile(flat,95);
clow  = min(clow,0);
chigh = max(chigh,0);
% absolute alternative:
% clow = min(flat); chigh = max(flat);

fig = figure('Units','normalized','Position',[0 0 1 1],'Color','w');
figDir  = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(figDir,'dir'), mkdir(figDir); end
pdfPath = fullfile(figDir, ...
    [sessionFileInfo.animal_name '_' RFMappingStimName '_RFpositionGrid.pdf']);

% tile layout with gutters
gX = 0.02; gY = 0.04;
w  = (1 - (numAz+1)*gX)/numAz;
h  = (1 - (numEl+1)*gY)/numEl;
gridCols = ceil(sqrt(nROI));
gridRows = ceil(nROI/gridCols);

% draw each az/el tile
for k = 1:nPos
    pos = psth(k).stimValue;
    ai = find(azimuths   == pos(1),1);
    ei = find(elevations == pos(2),1);
    if isempty(ai)||isempty(ei), continue; end

    left   = gX + (ai-1)*(w + gX);
    bottom = gY + (numEl-ei)*(h + gY);
    ax = axes('Position',[left,bottom,w,h]);

    vals = allVals(:,k);
    vals = vals ./ nanmedian(allVals,2);  % normalize by median across all positions
    [sortedVals, sortIdx] = sort(vals, 'descend', 'MissingPlacement', 'last');  % sort ROIs by response
    sortedVals(isnan(sortedVals)) = NaN;  % ensure NaNs stay NaN

    gridData = nan(gridRows, gridCols);
    gridData(1:numel(sortedVals)) = sortedVals;

%     vals     = allVals(:,k);
%     vals = vals./nanmedian(allVals,2); % divide by the mean across all positions. (changed to median)
%     gridData = nan(gridRows,gridCols);
%     gridData(1:numel(vals)) = vals;



%     imagesc(ax, gridData, [clow chigh]);
    imagesc(ax, gridData, [0 3]); % manually set the chigh and clow 
    axis(ax,'off','image');
    colormap(ax, whiteBlueRed);

    title(ax, sprintf('Az %d°, El %d°',pos(1),pos(2)), ...
        'FontSize',14,'Units','normalized','Position',[0.5,1,0],...
        'HorizontalAlignment','center');
end

% add single colorbar on the right, inside figure
% attach to invisible full‐figure axes:
axAll = axes('Position',[0 0 1 1],'Visible','off');
colormap(axAll, whiteBlueRed);
 
% chatgpt.. 
cbWidth  = 0.02;           % 3% of figure width
cbX      = 1 - cbWidth - 0.01;  % leave a 1% gap from the right edge
cbY      = gY;        % same bottom margin as your tiles
cbHeight = 1 - 2*gY;   % top & bottom margins
set(axAll, 'CLim', [0 3]); 
% draw the colorbar in that rectangle
c = colorbar('Position',[cbX, cbY, cbWidth, cbHeight]);
c.LineWidth     = 1.5;     % make the bar edge thicker
c.Label.String  = metric;
c.Label.FontSize= 12;

% get the barss position in normalised figure units
pos = c.Position;  % [x y w h]

% pick a point just to the left, half‐way up
lblX = pos(1) - 0.015;                % 1.5% to the left of the bar
lblY = pos(2) + pos(4)/2;             % vertically centered

% add your label there
text(lblX, lblY, 'Peak DF/F (after stim-onset)', ...
     'Units','normalized', ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','middle', ...
     'FontSize', 15,'Rotation', 90);


% 
text(0.5, 0.98, RFMappingStimName, ...
     'Units','normalized', ...
     'HorizontalAlignment','center', ...
     'FontSize',16, ...
     'FontWeight','bold', ...
     'Interpreter', 'none');


% save PDF
exportgraphics(fig, pdfPath, 'ContentType','vector');
close(fig);
end
