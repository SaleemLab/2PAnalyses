function summariseRFsession(sessionFileInfo, RFMappingStimName, azRange, elRange, tWindow)
% summariseRFsession: Compute and visualise RF metrics for a session
%   Loads psthData from the Response file and computes for each ROI:
%     Peak amplitude within given az/el window and time window
%     Time-to-peak latency
%   Then plots:
%     1) Scatter of latency vs. amplitude
%     2) Histogram of amplitudes
%     3) Prints count above threshold
%
% Usage:
%   summariseRFsession(sessionFileInfo, RFMappingStimName)
%   summariseRFsession(sessionFileInfo, RFMappingStimName, azRange, elRange)
%   summariseRFsession(sessionFileInfo, RFMappingStimName, azRange, elRange, tWindow)
%
% Inputs:
%   sessionFileInfo   - Struct containing fields .stimFiles and .Directories
%   RFMappingStimName - Name of the RF mapping stimulus (must match a stimFiles.name)
%   azRange           - [min max] Azimuth window in degrees (default: [-80 -25])
%   elRange           - [min max] Elevation window in degrees (default: [-30 30])
%   tWindow           - [t0 t1] Time window (s) relative to stimulus onset (default: [0 0.5])

if nargin < 3, azRange = [-80 -25]; end
if nargin < 4, elRange = [-30 30]; end
if nargin < 5, tWindow = [0 .5]; end

% --- Load Response.psthData ---
iStim = find(strcmp(RFMappingStimName, {sessionFileInfo.stimFiles.name}), 1);
if isempty(iStim)
    error('Stimulus "%s" not found in sessionFileInfo.', RFMappingStimName);
end
resp = load(sessionFileInfo.stimFiles(iStim).Response, 'response');
psthData = resp.response.psthData;  % struct array

% --- Compute RF metrics ---
[amp, ttp] = computeRFMetrics(psthData, azRange, elRange, tWindow);

% --- Plot peak amplitude vs latency scatter ---
figure('Name','RF Strength vs Latency');
scatter(ttp, amp, 30, amp, 'filled');
colormap(redWhiteBlue);
colorbar;
xlabel('Time to Peak (s)');
ylabel('Peak ΔF/F');
title(sprintf('%s: RF in Az[%d,%d] El[%d,%d]', sessionFileInfo.animal_name, azRange(1), azRange(2), elRange(1), elRange(2)));
set(gca,'FontSize',12);

% --- Histogram of amplitudes ---
figure('Name','Distribution of RF Strength');
histogram(amp,20);
xlabel('Peak ΔF/F');
ylabel('Number of ROIs');
title('RF Strength Distribution');
set(gca,'FontSize',12);

% --- Print count above threshold ---
threshold = 0.1;
nResponsive = sum(amp > threshold);
fprintf('ROIs with peak ΔF/F > %.2f: %d out of %d.\n', threshold, nResponsive, numel(amp));

end

%% Subfunction to compute metrics
function [amp, ttp] = computeRFMetrics(psthData, azRange, elRange, tWindow)
% computeRFMetrics: peak amplitude and latency for RF responses
% Inputs:
%   psthData  - struct array with fields .stimValue, .alignedResponses, .timeVector
%   azRange   - [min max] azimuth window
%   elRange   - [min max] elevation window
%   tWindow   - [t0 t1] time window for peak search

tVals = psthData(1).timeVector;
tMask = tVals >= tWindow(1) & tVals <= tWindow(2);

% Determine positions in az/el window
stimVals = vertcat(psthData.stimValue);
azMask = stimVals(:,1) >= azRange(1) & stimVals(:,1) <= azRange(2);
elMask = stimVals(:,2) >= elRange(1) & stimVals(:,2) <= elRange(2);
posMask = azMask & elMask;

nROI = size(psthData(1).alignedResponses, 1);
amp = nan(nROI,1);
ttp = nan(nROI,1);

for r = 1:nROI
    % Gather all trials x time for target positions
    data = [];
    for k = find(posMask)'
        A = squeeze(psthData(k).alignedResponses(r,:,:));
        if size(A,1) > size(A,2)
            A = A';
        end
        data = [data; A];
    end
    if isempty(data), continue; end
    meanTrace = mean(data, 1, 'omitnan');
    sub = meanTrace(tMask);
    if isempty(sub) || all(isnan(sub)), continue; end
    [amp(r), idx] = max(sub);
    tm = tVals(tMask);
    ttp(r) = tm(idx);
end
end