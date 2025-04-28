function plotRFandTemporalPerROI(sessionFileInfo, RFMappingStimName, prestimWindow, ROIsPerPage)
% plotRFandTemporalPerROIsPaged: Plots spatial RF and normalized temporal traces for ROIs.
% Layout: 5 ROIs per row -> 10 tiles (RF + temporal) each row.
% Inputs:
%   sessionFileInfo - struct with .stimFiles and .Directories.save_folder
%   RFMappingStimName - string matching sessionFileInfo.stimFiles.name
%   prestimWindow - [t0 t1] for baseline normalization (default [-0.5 0])
%   ROIsPerPage - number of ROIs per page (multiple of 5; default 30)

if nargin<3, prestimWindow = [-0.5 0]; end
if nargin<4, ROIsPerPage = 30; end

% Load PSTH data
iStim = find(strcmp(RFMappingStimName,{sessionFileInfo.stimFiles.name}),1);
assert(~isempty(iStim),'Stimulus not found');
resp = load(sessionFileInfo.stimFiles(iStim).Response,'response');
psth = resp.response.psthData;

% Azimuth/Elevation grid and important positions
stimVals = vertcat(psth.stimValue);
az = unique(stimVals(:,1));
el = unique(stimVals(:,2));
important = [-70 -22; -70 0; -70 22; -30 -22; -30 0; -30 22];

tVec = psth(1).timeVector;
nROI = size(psth(1).alignedResponses,1);
pages = ceil(nROI/ROIsPerPage);

% Prepare output PDF
figDir = fullfile(sessionFileInfo.Directories.save_folder,'Figures');
if ~exist(figDir,'dir'), mkdir(figDir); end
pdfPath = fullfile(figDir,sprintf('%s_%s_RFpages.pdf',sessionFileInfo.animal_name,RFMappingStimName));

for pg = 1:pages
    idxStart = (pg-1)*ROIsPerPage + 1;
    idxEnd   = min(pg*ROIsPerPage, nROI);
    currentROIs = idxStart:idxEnd;
    nCurr = numel(currentROIs);
    nRows = ceil(nCurr/5);

    % Create figure and layout
    fig = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
    t = tiledlayout(nRows,10,'TileSpacing','compact','Padding','compact');
    sgtitle(t,sprintf('%s | ROIs %d-%d',sessionFileInfo.animal_name,idxStart,idxEnd), 'FontWeight','bold');

    % Loop through selected ROIs
    for ii = 1:nCurr
        roiIdx = currentROIs(ii);
        % Determine tile index for RF map and temporal
        rfTile = (ii-1)*2 + 1;
        tmpTile = rfTile + 1;

        % Spatial RF map
        ax1 = nexttile(rfTile);
        RFmap = nan(numel(el),numel(az));
        postMask = tVec>0;
        for k = 1:numel(psth)
            pos = psth(k).stimValue;
            a = find(az==pos(1)); e = find(el==pos(2));
            if isempty(a)||isempty(e), continue; end
            tr = squeeze(psth(k).alignedResponses(roiIdx,:,:));
            if size(tr,1)>size(tr,2), tr=tr'; end
            mt = mean(tr,1,'omitnan');
            RFmap(e,a) = mean(mt(postMask),'omitnan');
        end
        imagesc(ax1,az,el,RFmap); axis(ax1,'xy');
        xline(ax1,0,'k:'); yline(ax1,0,'k:');
        colormap(ax1,redWhiteBlue);
        set(ax1,'XTick',[],'YTick',[]);
        title(ax1,sprintf('ROI %d RF',roiIdx),'FontSize',8);

        % Temporal trace
        ax2 = nexttile(tmpTile);
        hold(ax2,'on');
        % Shade stimulus period 0 to 1s
        yl = [-1 2];
        patch(ax2,[0 1 1 0],[yl(1) yl(1) yl(2) yl(2)],[0.8 0.8 0.8],...
            'FaceAlpha',0.3,'EdgeColor','none');
        for k = 1:numel(psth)
            pos = psth(k).stimValue;
            tr = squeeze(psth(k).alignedResponses(roiIdx,:,:));
            if size(tr,1)>size(tr,2), tr=tr'; end
            if isempty(tr), continue; end
            mt = mean(tr,1,'omitnan');
            baseMask = tVec>=prestimWindow(1)&tVec<=prestimWindow(2);
            F0 = mean(mt(baseMask),'omitnan');
            nmt = (mt - F0)/F0;
            clr = 'k';
            if any(ismember(important,pos,'rows')), clr='g'; end
            plot(ax2,tVec,nmt,'Color',clr,'LineWidth',0.8);
        end
        xline(ax2,0,':','Color',[0.5 0.5 0.5]);
        ylim(ax2,yl);
        set(ax2,'XTick',[],'YTick',[]);
        title(ax2,'Temporal','FontSize',8);
    end

    % Export page
    exportgraphics(fig,pdfPath,'Append',true,'ContentType','vector');
    close(fig);
end
end
