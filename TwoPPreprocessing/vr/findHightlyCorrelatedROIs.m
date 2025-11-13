function findHightlyCorrelatedROIs(sessionFileInfo, signalToUse, applySmoothing, plotFlag)

if nargin < 7; signalToUse = 'dFFNeuropilCorrected'; end
if nargin < 8; applySmoothing = true; end
if nargin < 9; plotFlag = false; end
% if nargin < 9; plotFlag = true; end

%% Find the right VRStimIdx and load response structure 
thisVRStim = selectVRStimIndex(sessionFileInfo);
if ~isempty(thisVRStim)
    responseFilePath = sessionFileInfo.stimFiles(thisVRStim).Response;
    
end