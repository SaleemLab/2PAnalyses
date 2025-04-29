# 2PAnalyses
Code for analysing 2 photon imaging data
#### By: Sonali
#### Last update: ???

## Step 0
#### 0.1 Run Suite2p (concatenated across all stimuli recorded in one session)


#### 0.2 Get all the session information.  
<br> <br>

## Step 1: Get all the session information.

#### 1.1 Get all the file paths in the session
get2PsessionFilePaths.m  
Usage: [sessionFileInfo] = get2PsessionFilePaths(animal_name, session_name)
Example: sessionFileInfo = get2PsessionFilePaths('M24048', '20240816')
Output: M24048_20240816_sessionFileInfo.mat 
<br> <br>

#### 1.2 Get 2P metadata from tif file. This will include the imaging parameters. 
get2PMetadata.m
Usage: [sessionFileInfo] = get2PMetadata(sessionFileInfo)
Example: sessionFileInfo = get2PMetadata(sessionFileInfo)
Output: M24048_20240816_2PMetaData_StimName.mat (one file for each stim)
<br> <br>

#### 1.3 Counts the length of tif files across all stim in one recording session
This is needed for the next steps, but it takes a long time to run.
### "TODO: Change strategy here to the one used in processRawTifFiles.m to make it faster."
get2PFrameTimes.m 
Usage: [sessionFileInfo] = get2PFrameTimes(sessionFileInfo);
Example: sessionFileInfo = get2PFrameTimes(sessionFileInfo);
Output: Updates the stim frame run in the sessionFileInfo  
<br> <br>

#### 1.4 Load the peripheral files that are common to all stimuli: Wheel, PD, Quad, 2P. Saves as a new .m file for each stimulus  
processPeripheralFiles.m 
Usage: [sessionFileInfo] = processPeripheralFiles(sessionFileInfo)
Example: sessionFileInfo = processPeripheralFiles(sessionFileInfo)
Output: animalName_session_PeripheralData_stimName.m 
<br> <br>


## Step 2: Get intermediate files with events for each stimulus and save

#### 2.1 Gets all the times for all planes, also getting the F/fneu/ops/spks/iscell from suit2p for each stimulus; 
Trims bonsai time to match the suite2p plane lengths. Saves one .m file for each stimulus.

mergeBonsaiSuite2pFiles.m 
Usage: [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo) 
Example: sessionFileInfo = mergeBonsaiSuite2pFiles(sessionFileInfo) 
Output: animalName_session_2pData_stimName.mat

#### 2.1.1 Uses this function: [plane_data, plane_data_new] = get_bonsai_twopframetimes_by_planes(filepath, nplanes) to get the planetimes  


#### 2.2.1 Gets additional bonsai data from csv files. Currently still raw and not interpolated. 
getVRBonsaiFiles.m 
Usage: VRInfo = getVRBonsaiFiles(sessionFileInfo)
Output: VRInfo.mousePos; VRInfo.trialInfo; VRInfo.isVRStim

------- CURRENTLY DONE UP TO HERE 18.2.25 ------------------------------


#### 2.2.2 For DirTuning 

#### 2.2.3 For Dot fields

#### 2.2.4 For Sparse Noise 

## Step 3: Getting into stim specific analyses

For VR - interpolation is the first step
