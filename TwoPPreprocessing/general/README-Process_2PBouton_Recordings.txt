Pipeline for z-motion corrected boutons (multi-plane; multi-channel) 

---- Step 0 ----- 
Tim tiff files and extract site2p data 

0.1) ProcessRawTifFile. This loads the last tif file for each stim
and trims the last tif to have equal number of frames for each channel and across planes. 
Runs as an isolated function. Finds root directories and trim tif file
processRawTifFile.m 
Usage: processRawTifFile('M25012', '20250210')
Outout: animalname_session_stim_0000x_trimmed.tif 
This will also create an archived folder with two copies of the tif file that was 'opened' and the other copy is untouched. 

0.2) Run Suite2p z-registration (in Python, 'Data' repo) to curate a single plane with the z-motion 
correction 

------ Step 1 ------
Get all the session information, pripheral information and metadata 

1.1) Get all the file paths in the session
get2PsessionFilePaths.m  
Usage: [sessionFileInfo] = get2PsessionFilePaths(animal_name, session_name)
Example: sessionFileInfo = get2PsessionFilePaths('M24048', '20240816')
Output: M24048_20240816_sessionFileInfo.mat 

1.2) Get 2P metadata from tif file. This will include the imaging parameters. 
get2PMetadata.m
Usage: [sessionFileInfo] = get2PMetadata(sessionFileInfo)
Example: sessionFileInfo = get2PMetadata(sessionFileInfo)
Output: M24048_20240816_2PMetaData_StimName.mat (one file for each stim)
 
1.4) Counts the length of tif files across all stim in one recording session
This is needed for the next steps, but it takes a long time to run.
## TODO: Change strategy here to the one used in processRawTifFiles.m to make
it faster.
get2PFrameTimes.m 
Usage: [sessionFileInfo] = get2PFrameTimes(sessionFileInfo, nPlanesPreZCorrection);
Example: sessionFileInfo = get2PFrameTimes(sessionFileInfo, 8);
Output: Updates the stim frame run in the sessionFileInfo

Version 2: @Aman
% [sessionFileInfo] = get2PFrameTimes_TwoChannels(sessionFileInfo, nPlanesPreZCorrection, nChannels)
%
% This function calculates how many frames belong to each imaging plane
% for each stimulus, taking into account the number of imaging channels.
%
% Arguments:
%   - sessionFileInfo: struct containing session and stimulus information
%   - nPlanesPreZCorrection: number of planes (optional)
%   - nChannels: number of imaging channels (optional, default = 1)  


1.5) Load the peripheral files that are 
common to all stimuli: Wheel, PD, Quad. Saves as a new .m file for each 
stimulus  
processPeripheralFiles.m 
Usage: [sessionFileInfo] = processPeripheralFiles(sessionFileInfo)
Example: sessionFileInfo = processPeripheralFiles(sessionFileInfo)
Output: animalName_session_PeripheralData_stimName.m 

Step 2:
Get intermediate files with events for each stimulus and save

2.1) Gets all the times for all planes, also getting the F/fneu/ops/spks/iscell from suit2p for each stimulus; 
Trims bonsai time to match the suite2p plane lengths. Saves one .m file for each stimulus.

mergeBonsaiSuite2pFiles.m 
Usage: [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo) 
Example: sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo) 
Output: animalName_session_2pData_stimName.mat

2.1.1) Uses this function: [plane_data, plane_data_new] = get_bonsai_twopframetimes_by_planes(filepath, nplanes)
to get the planetimes 

 ----- Step 3 -----
Stim specfific Bonsai Data 

3.2) VR
getVRBonsaiFiles.m
Get the bonsai VR mouse position and trial info tables and saves in a 
Usage: [bonsaiData] = getVRBonsaiFiles(sessionFileInfo)
Output: bonsaiData.mousePos; bonsaiData.trialInfo; bonsaiData.isVRStim 

3.3) getTuningStimEventsBonsaiFile.m 

DirTuning, DotMotion_SpeedTuning and other tuning files
Get the Bonsai stim-events/Trial Parameters table (if saved) and add them to a .mat file.
CAUTION! New stimuli will need to be defined as a new case. 
Add all stim identities to bonsaiData.stimType 
Usage: [bonsaiData, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, stimName, stimTypeTableName)  
Example: [bonsaiData, sessionFileInfo] = getTuningStimEventsBonsaiFile(sessionFileInfo, 'M25041_RFMapping_20250417_00001', 'StimulusParams');
Output: bonsaiData, sessionFileInfo  

3.4) Sparse Noise 
getSparseNoiseBonsaiFiles.m
Usage: [bonsaiData] = getSparseNoiseBonsaiFiles(sessionFileInfo, gridSize)
% Import and save stimulus matrix (position on screen) from Sparse Noise Bin file 

3.5) getStimTimes.m 
% Extracts stimulus ON and OFF times using PD (and QuadState as a reference, if required)
% and removes stimulus events after two-photon imaging has stopped.
[bonsaiData, sessionFileInfo] = getStimTimes(sessionFileInfo, stimName, thresholdOn, thresholdOff, useQuadState, plotFlag)
(Look up .m file for more information on inputs and usage) 
(@Aman - move onARDTimes and offARDTimes to a new .mat file called response or save in bonsaiData?)

---- Step 4 ---- CURRENTLY HERE 23.02.25
response

4.1) @Aman ? 
[responseFrameIdx, responseFrameRelTimesIdx] = get2PFramesByTrial(onARDTimes, twoPData, preStimTime, postStimTime, offARDTimes)
OR
[responseFrameIdx, responseFrameRelTimes] = get2PFramesByTrialV2(onARDTimes, twoPData, preStimTime, postStimTime, offARDTimes, nExpected)
OR 
[response] = get2PFramesByTrialV3(sessionFileInfo, stimName, useoffARDTimes, preStimTime, postStimTime)


4.2) getTrialGroups.m % Changed from response to bonsaiData (as response is currently structures per plane) 
Groups trials based on stimulus type from Bonsai data. 
[bonsaiData] = getTrialGroups(sessionFileInfo, stimName)
@Aman - currently saving this in response? 
Outputs:
trialGroups.value - Identifier for each stimulus type (e.g., 0, 45, 90, etc.).
trialGroups..stimTypeName - Name of the stimulus category (uses `stimName` if not in Bonsai data)
trialGroups..trials - Indices of trials corresponding to each stimulus type.
 


