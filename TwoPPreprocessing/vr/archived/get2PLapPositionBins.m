function [response] = get2PLapPositionBins(sessionFileInfo,plotFlag)
% Goal: Index of 2P frames lap position bins for all rois 
%       relative time? 
%% Load processed Bonsai, Peripheral and 2P data 
for iStim = 1:length(sessionFileInfo.stimFiles)
    bonsaiData.isVRstim(iStim) = strcmp('VRCorr',sessionFileInfo.stimFiles(iStim).name);
end

iStim = find(bonsaiData.isVRstim==1);

if exist(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).BonsaiData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'file') && ...
        exist(sessionFileInfo.stimFiles(iStim).Response, 'file')
    load(sessionFileInfo.stimFiles(iStim).processedMergedBonsaiSuite2pData, 'processedTwoPData')
    load(sessionFileInfo.stimFiles(iStim).BonsaiData, 'bonsaiData');
    load(sessionFileInfo.stimFiles(iStim).processedPeripheralData, 'peripheralData');
    load(sessionFileInfo.stimFiles(iStim).Response, 'response')
else
    error('Response, bonsaiData, peripheralData or processedTwoPData not found for this session');
end

%% Find lap position bins 
speedFilter = response.wheelSpeed > 1; 
binEdges = 0:140;  % Position bin edges (0 cm to 140 cm, 1 cm per bin)
binCentres = 0.5:1:139.5;  % Position bin centers 
numROIs = size(processedTwoPData.F, 1);  % Number of neurons; remove 
numBins = length(binCentres);  % Number of position bins
lapPosition2PIdx = zeros(numROIs, length(response.completedStartTimes), numBins);  % Output matrix for neuron x laps x position bins

% Loop over each lap
for thisLap = 1:length(response.completedStartTimes)
    display(['Running Lap: ' num2str(thisLap)]);
    % @Aman - using twoP time here?
    lapTvecIdx = find(processedTwoPData.ArduinoTime >= response.completedStartTimes(thisLap) & ...
        processedTwoPData.ArduinoTime <= response.completedEndTimes(thisLap));
    if isempty(lapTvecIdx)
        disp(['Frame index missing from lap ' thisLap])
        continue;  
    end
    for ineuron = 1:numROIs
        lapPosition = response.mouseVirtualPosition(lapTvecIdx);
%       Assign each frame to a position bin using the defined bin edges
        position_index = discretize(lapPosition, binEdges);
        temp_spatial_responses = zeros(1, length(binCentres));
        
        for iBin = 1:length(binCentres)
            % Get indices for frames in this bin
            bin_index = position_index == (iBin) & speedFilter(lapTvecIdx);
            no_bin_index = sum(bin_index);
            if no_bin_index > 0
                temp_spatial_responses(iBin) = sum(smoothed_responses(ineuron, bin_index)) / no_bin_index;
            end
        end
        lapPosition2PIdx(ineuron, thisLap, :) = temp_spatial_responses;
    end
end

%% Save:
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
save(sessionFileInfo.stimFiles(iStim).Response, 'response'); 
end