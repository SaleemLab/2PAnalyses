function [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo) 
% Example usage: [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo) 
% Outputs: Updates sessionFileInfo with filepaths and mouse_session_2pData.mat file 
% Sonali and Aman - Jan 2025

for iStim = 1:length(sessionFileInfo.stimFiles)
    
    % Add the merged .mat file to the sessionFileInfo 
    stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_2pData' '_' sessionFileInfo.stimFiles(iStim).name '.mat'];
    sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData = fullfile(sessionFileInfo.Directories.save_folder,stimFileName);
    
    % find the last plane's size to trim bonsai and suite2p data
    if strcmp(sessionFileInfo.suite2pFiles(end).planeName,'plane')
        lastPlane_FrameRun = sessionFileInfo.stim_framerun{[sessionFileInfo.suite2pFiles(end).planeName '0'],sessionFileInfo.stimFiles(iStim).name}{1};
        numPlanes = sessionFileInfo.origNPlanes;
    else
        lastPlane_FrameRun = sessionFileInfo.stim_framerun{sessionFileInfo.suite2pFiles(end).planeName,sessionFileInfo.stimFiles(iStim).name}{1};
        numPlanes = sessionFileInfo.numPlanes;
    end
    trimLength = lastPlane_FrameRun(2)-lastPlane_FrameRun(1) + 1;
    
%     % Pick out VRCorr if present and process differently.
%     bonsaiData.isVRstim(iStim) = strcmp('VRCorr',sessionFileInfo.stimFiles(iStim).name);
%     if iStim == bonsaiData.isVRstim
% 
%     end 
    % load the Bonsai 2p plane times 
    twop_filepath     = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, '2P');
    planeTimes_table  = get_bonsai_twopframetimes_by_planes(twop_filepath, numPlanes);
   
    %% Suite2p and Bonsai
    for iPlane = 1:sessionFileInfo.numPlanes
        planeName = sessionFileInfo.suite2pFiles(iPlane).planeName;
        twoPData(iPlane).planeName = planeName;
        
        % ---------- Suite2p -----------
        % This filepath contains all the suite2p data stored in a .mat file
        fAll_filepath = findFile(sessionFileInfo.suite2pFiles(iPlane).planes, 'fall');
          
        % Check if the is empty (i.e., file not found)
        if isempty(fAll_filepath)
           warning([twoPData(iPlane).planeName, ' fall.mat file has not been generated from suite2p or this plane has no data']);
           continue;  % Skip this plane and move on to the next one
        end 
        
        % Load the Suite2p data for the current plane
        fAll                    = load(fAll_filepath);
        twoPData(iPlane).ops    = fAll.ops;
        twoPData(iPlane).iscell = fAll.iscell;
        twoPData(iPlane).stat   = fAll.stat;
        twoPData(iPlane).redcell= fAll.redcell;

        % Pick out the stimulus frameRun for the corresponding stimulus 
        if strcmp(planeName,'plane')
            twoPData(iPlane).frameRun = sessionFileInfo.stim_framerun{[planeName '0'],sessionFileInfo.stimFiles(iStim).name}{1};
            currentPlaneInfo = twoPData.ops.current_plane(twoPData(iPlane).frameRun(1):(twoPData(iPlane).frameRun(1)+trimLength-1));
        else
            twoPData(iPlane).frameRun = sessionFileInfo.stim_framerun{planeName,sessionFileInfo.stimFiles(iStim).name}{1};
        end
    
        
        % Trim across all planes using the last plane's max length and save
        trimIndices             = twoPData(iPlane).frameRun(1):(twoPData(iPlane).frameRun(1)+trimLength-1);
        twoPData(iPlane).F      =    fAll.F(:, trimIndices);
        twoPData(iPlane).spks   = fAll.spks(:, trimIndices);
        twoPData(iPlane).Fneu   = fAll.Fneu(:, trimIndices);  

        % % ---------- Bonsai -----------
        if strcmp(planeName,'plane')
            planeTimes_trim = planeTimes_table.plane0(1:trimLength,:);
            for iTime = 1:trimLength
                planeTimes_trim(iTime,:) = planeTimes_table.([planeName num2str(currentPlaneInfo(iTime))])(iTime, :);
            end
        else
            planeTimes_trim = planeTimes_table.(planeName)(1:trimLength, :);
        end
        twoPData(iPlane).TwoPFrameTime = planeTimes_trim.TwoPFrameTime;
        twoPData(iPlane).BonsaiTime = planeTimes_trim.BonsaiTime;
        twoPData(iPlane).RenderFrameCount = planeTimes_trim.RenderFrameCount;
        twoPData(iPlane).LastSyncPulseTime = planeTimes_trim.LastSyncPulseTime;
        twoPData(iPlane).ArduinoTime = planeTimes_trim.ArduinoTime;
    end

save(sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData, 'twoPData');
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');

end