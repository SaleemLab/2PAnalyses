Pipeline for z-motion corrected boutons and somas for VR sessions (multi-plane; multi-channel-TODO) 

---- Step 0: Trim tifs and run Suite2P ----- 
Tim tiff files and extract site2p data 

0.1) processRawTifFile.m (general)

ProcessRawTifFile. This loads the last tif file for each stim
and trims the last tif to have equal number of frames for each channel and across planes. 
Runs as an isolated function. Finds root directories and trim tif file
Usage: processRawTifFile('M25012', '20250210')
Outout: animalname_session_stim_0000x_trimmed.tif 
This will also create an archived folder with two copies of the tif file that was 'opened' and the other copy is untouched. 

0.2) Run Suite2p z-registration 
In Python, 'Data' repo to curate a single plane with the z-motion 
correction 

------ Step 1: SessionInfo and 2PFrame Times ------
Get all the session information, pripheral information and metadata 

1.1) get2PsessionFilePaths.m   (general)

Get all the file paths in the session
Usage: [sessionFileInfo] = get2PsessionFilePaths(animal_name, session_name, stim_list, rerun_process, fileNameAddition)
Example: sessionFileInfo = get2PsessionFilePaths('M24048', '20240816')
Output: M24048_20240816_sessionFileInfo.mat 

1.2) get2PMetadata.m (general)

Get 2P metadata from tif file. This will include the imaging parameters. =
Usage: [sessionFileInfo] = get2PMetadata(sessionFileInfo)
Example: sessionFileInfo = get2PMetadata(sessionFileInfo)
Output: M24048_20240816_2PMetaData_StimName.mat (one file for each stim)
 
1.4) get2PFrameTimes_TwoChannels.m (general)

Counts the length of tif files across all stim in one recording session
This is needed for the next steps, but it takes a long time to run.
## TODO: Change strategy here to the one used in processRawTifFiles.m to make
it faster.
Usage: [sessionFileInfo] = get2PFrameTimes_TwoChannels(sessionFileInfo, nPlanesPreZCorrection, nChannels);
Example: sessionFileInfo = get2PFrameTimes(sessionFileInfo, 8, 1);
Output: Updates the stim frame run in the sessionFileInfo  


 ----- Step 2: Load, save and restructure data streams -----
Load Suite2P, 2P-FrameTimes (saved in Bonsai), Bonsai and Peripheral Data 

2.1) processPeripheralFiles.m (general)

Load the peripheral files that are 
common to all stimuli: Wheel, PD, Quad. 
Usage: [sessionFileInfo] = processPeripheralFiles(sessionFileInfo)
Example: sessionFileInfo = processPeripheralFiles(sessionFileInfo)
Output: animalName_session_PeripheralData_stimName.mat 

2.1) getVRBonsaiFiles.m (vr) / temp include VRStim name

Get the bonsai VR mouse position and trial info tables and saves in a 
Usage: [bonsaiData, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, VRStimName)
Example:  [bonsaiData, sessionFileInfo] = getVRBonsaiFiles(sessionFileInfo, 'M25041_VRCorr_20250413_00001')
Output: bonsaiData.mousePos; bonsaiData.trialInfo; bonsaiData.isVRStim 


2.2) mergeBonsaiSuite2pFiles.m (general function)

Suite2P raw and trimmed to match length across planes (if multi)
Usage: [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo) 
Example: sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo) 
Output: animalName_session_2pData_stimName.mat

  
 ----- Step 3: Align to common time base (and correct bonsai lag) -----


(3.1.1) findBonsaiPeripheralLag.m (general)

[bonsaiData, sessionFileInfo] = findBonsaiPeripheralLag(sessionFileInfo, method, sampling_rate)
Usage: [bonsaiData, sessionFileInfo] = findBonsaiPeripheralLag(sessionFileInfo, method, samplingRate, StimName)
Example: [bonsaiData, sessionFileInfo] = findBonsaiPeripheralLag(sessionFileInfo, 1, 60)
Outputs: bonsaiData.LagInfo (depending on method used): 
        .method         - Method used to find lag. 
        .xcorrBestLag   - Frames to shift.  
        .samplingRate   - Use as input in this function.
     

3.1.2) alignVRBonsaiToPeripheralData.m (Currently lag shifts only ArduinoTime)  (vr) / temp include VRStim name

Usage: [bonsaiData, sessionFileInfo]  = alignVRBonsaiToPeripheralData(sessionFileInfo,VRStimName,plotFlag)
Example: [bonsaiData, sessionFileInfo]  = alignVRBonsaiToPeripheralData(sessionFileInfo, 'M25041_VRCorr_20250413_00001')
Outputs: 
   - bonsaiData (struct):
       - .MousePos / .Quadstate / .TrialInfo:
           - .correctedArduinoTime: lag-corrected raw time vector.
       - .LagInfo
           - .lagShift: Lag shift used (in seconds). 


------- Step 4: Interpolate to common times ---------
****** Aman - notes:
A) Define sampling rate and times, and range for analyses 
Start time: First 2p frame 
End time: Last 2p frame 
(calculate this across all planes)
sampleTimes: startTime:(1/samplingRate):endTime

Resample EVERYTHING
Call the interpolation function in sequence

mousepos: linear
trailID: nearest
photodiode:
quad:
Wheel:
for each cell: F, Spks, Fneu
interpolation function: Input: (ArduinoTime (rawSampleTime), Value, sampleTimes, Method): Output: NewValues (at sampleTimes)
******

4.0) resamplAndAlignVR_BonsaiPeripheralSuite2P.m (Resample all data streams to common time base) (vr)

resamplAndAlignVR_BonsaiPeripheralSuite2P.m / temp stim name

Usage: [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo, samplingRate, mainTimeToUse, VRStimName, plotFlag)
Example: [processedTwoPData, bonsaiData, peripheralData, sessionFileInfo] = resamplAndAlignVR_BonsaiPeripheralSuite2P(sessionFileInfo,60,'TwoPFrameTime', 'M25041_VRCorr_20250413_00001', true) 
Outputs:
  processedTwoPData : struct
      Contains resampled F, Fneu, spks, frame times, and ROI metadata.
  bonsaiData : struct
      Bonsai-tracked signals (e.g., mouse position, trial info, quadstate), all corrected for lag and resampled.
  peripheralData : struct
      Peripheral signals (e.g., photodiode, wheel), resampled to the same timebase.
         

---- Step 5: Response: Extract lap position activity ----  


5.0) extractVRBonsaiPeripheralInfo.m (Split 5.0 into two functions? Bonsai / Peripheral and move peripheral to general?) (currently vr)

Extracts wheel speed, virtual position, and lap-related info
from aligned Bonsai and peripheral data during 2P-VR Aman's Classical Corridor.
Handles lap classification (completed/aborted) and optionally plots lap timing.

Usage: [response, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo,  VRStimName, plotFlag)
Example:[response, sessionFileInfo] = extractVRAndPeripheralData(sessionFileInfo,  'M25041_VRCorr_20250413_00001', true)
response. : struct
%       Contains key behavioral and timing variables for downstream lap-by-lap analysis:
%       - wheelSpeed             : real-time wheel speed (cm/s)
%       - mouseVirtualPosition   : virtual track position (1–140 cm)
%       - trackIDFromMousePosition : track IDs inferred from mouse position (usually 1)
%       - mouseRecordedPosition  : raw position signal (-1141 to -1000)
%       - trackIDs               : track ID per lap (all 1s in classical VR corridor)
%       - lapCount               : unified lap index
%       - blockIDs               : cumulative index of track switches (e.g., block transitions; NA)
%       - trialType              : task/trial type ID per lap 
%                                  (0-NoTask; 1-Passive; 2-Hybrid)
%       - completedLaps          : indices of completed laps
%       - abortedLaps            : indices of aborted laps
%       - startTimeAll           : Bonsai start times for all laps
%       - endTimeAll             : parsed lap end time (based on position trace)
%       - completedStartTimes    : start times for completed laps
%       - completedEndTimes      : end times for completed laps
%


5.1.1) get2PFrameLapPositionBins.m  (vr)

   For each ROI, lap, and spatial bin (1 cm), this function finds
   the corresponding two-photon (2P) frame indices and relative
   times (from lap start). Only includes frames where wheel speed > 1 cm/s.
   Both cells and non-cells are included here 

Usage: [response, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, VRStimName)
Example: [response, sessionFileInfo] = get2PFrameLapPositionBins(sessionFileInfo, 'M25041_VRCorr_20250413_00001')
Output:
Output:
  response : struct (updated)
      Adds:
        - lapPosition2PFrameIdx{ROI, lap, bin} : 2P frame indices per ROI per lap per bin
        - lapPositionRelativeTime{lap, bin}    : time relative to lap start (only once per bin) (could exclude)


5.1.2) getLapPositionActivity.m 

  Extracts mean binned fluorescence/activity values per lap from 2P data.
  Only includes ROIs labeled as "cells" and only for completed laps.
  Optionally applies Gaussian temporal smoothing across position bins.

Usage:  [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, signalField, applySmoothing, VRStimName, onlyIncludeROIs)
Example: [response, sessionFileInfo] = getLapPositionActivity(sessionFileInfo, 'F', false, 'M25041_VRCorr_20250413_00001', true)
Outputs:
  response : struct (updated)
      Adds the following fields:
        - lapPositionActivity : [nCells x nLaps x nBins] mean activity per bin
        - cellROIs             : indices of cell ROIs
        - signalUsed           : which signal type was used
        - smoothingApplied     : whether smoothing was applied

5.1.2.1) plotSortedPopulationResponse.m 

plotSortedPopulationResponse(sessionFileInfo, response, applySmoothing)