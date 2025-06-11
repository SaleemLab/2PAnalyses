function [response, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, VRStimName)
%  TODO: remove the relative time and
%   For each ROI, lap, and spatial bin (1 cm), this function finds
%   the corresponding two-photon (2P) frame indices and relative
%   times (from lap start). Only includes frames where wheel speed > 1 cm/s.
%   Both cells and non-cells are included here
%
% Inputs:
%   sessionFileInfo : struct
%       Structure containing paths to processed files including
%       Suite2P data, Bonsai data, peripheral signals, and the response struct.
%
% Output:
%   response : struct (updated)
%       Adds:
%         - lapPosition2PFrameIdx{ROI, lap, bin} : 2P frame indices per ROI per lap per bin
%         - lapPositionRelativeTime{lap, bin}    : time relative to lap start (only once per bin)
%
% Example usage:
%   response = get2PFrameLapPositionBins(sessionFileInfo);
%
% Aman and Sonali - April 2025


%% Load processed data
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp(VRStimName, sessionFileInfo.stimFiles(iStim).name);
end
iStim = find(bonsaiData.isVRstim == 1);

if exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).Response, 'file')

    load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'processedTwoPData')
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
    load(sessionFileInfo.stimFiles(iStim).Response, 'response')
else
    error('Missing processed files: Response, BonsaiData, PeripheralData, or TwoPData.');
end

%% Define bins and setup
posBinEdges = 0:140;
posBinCentres = 0.5:1:139.5;
numPosBins = length(posBinCentres);
numROIs = size(processedTwoPData.F, 1);
nLaps = length(response.completedStartTimes);
timeVec = processedTwoPData.(processedTwoPData.resample2PTimeUsed);
speedFilter = response.wheelSpeed > 1 & response.wheelSpeed < 100;

lapPosition2PFrameIdx = cell(numROIs, nLaps, numPosBins);
lapPositionRelativeTime = cell(nLaps, numPosBins);  % Only stored once per lap × bin

%% Loop over laps
for thisLap = 1:nLaps
    disp(['Running Lap: ' num2str(thisLap)]);

    lapStart = response.completedStartTimes(thisLap);
    lapEnd = response.completedEndTimes(thisLap);
    lapFrameIdx = find(timeVec >= lapStart & timeVec <= lapEnd);

    if isempty(lapFrameIdx)
        disp(['No 2P frames found for lap ' num2str(thisLap)]);
        continue;
    end

    lapPosition = response.mouseVirtualPosition(lapFrameIdx);
    positionIdx = discretize(lapPosition, posBinEdges);

    for thisBin = 1:numPosBins
        binMask = (positionIdx == thisBin) & speedFilter(lapFrameIdx);
        frameIdxInBin = lapFrameIdx(binMask);

        % Store relative time only once
        if ~isempty(frameIdxInBin)
            frameTimeInBin = timeVec(frameIdxInBin);
            lapPositionRelativeTime{thisLap, thisBin} = frameTimeInBin - lapStart;
        else
            lapPositionRelativeTime{thisLap, thisBin} = [];
        end

        % Store frame indices for all ROIs (same for all, reused)
        for thisROI = 1:numROIs
            lapPosition2PFrameIdx{thisROI, thisLap, thisBin} = frameIdxInBin; %%%% do you need a cell, or just 3D matrix?
        end
    end
end

%% Save to response
response.lapPosition2PFrameIdx = lapPosition2PFrameIdx;
response.lapPositionRelativeTime = lapPositionRelativeTime;

%% Example Plot (remove)
% roiIdx = 1;
% lapIdx = 5;
% meanSignalPerBin = nan(1, numPosBins);
% for binIdx = 1:numPosBins
%     frameIdx = response.lapPosition2PFrameIdx{roiIdx, lapIdx, binIdx};
%     if ~isempty(frameIdx)
%         signal = processedTwoPData.F(roiIdx, frameIdx);
%         meanSignalPerBin(binIdx) = mean(signal);
%     end
% end
% 
% figure;
% plot(posBinCentres, meanSignalPerBin, 'k', 'LineWidth', 1.5);
% xlabel('Position (cm)');
% ylabel('Mean F');
% title(sprintf('PSTH ROI %d, Lap %d', roiIdx, lapIdx));

%% Save
disp('Saving response...');
save(sessionFileInfo.stimFiles(iStim).Response, 'response', '-v7.3');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

end
