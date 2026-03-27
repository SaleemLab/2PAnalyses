function [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo)
% Merges Bonsai and Suite2p data for each stimulus in a session.
% 
% NOTE TO SELF: Data is split here because the rest of the pipeline 
% requires this specific structure for plane-wise analysis.
% CHECK this function 
for iStim = 1:length(sessionFileInfo.stimFiles)
    stimName = sessionFileInfo.stimFiles(iStim).name;
    
    try
        fprintf('Processing mergeBonsaiSuite2PFiles for %s\n', stimName);
        
        stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_2pData' '_' stimName '.mat'];
        merged_data_filepath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
        
        % Determine Plane/Trim logic % previously added NAN padding and
        % this caused many issues further up the pipeline during smoothning
        % step for example. 
        if strcmp(sessionFileInfo.suite2pFiles(end).planeName, 'plane_z') % changed to planez instead of plane 
            lastPlane_FrameRun = sessionFileInfo.stim_framerun{['plane' '0'], stimName}{1};
            numPlanes = sessionFileInfo.origNPlanes;
            if numPlanes == 8
                sessionFileInfo.flybackPlanes = [0, 1, 2];
            end
        else
            lastPlane_FrameRun = sessionFileInfo.stim_framerun{sessionFileInfo.suite2pFiles(end).planeName, stimName}{1};
            numPlanes = sessionFileInfo.numPlanes;
        end
        trimLength = lastPlane_FrameRun(2) - lastPlane_FrameRun(1) + 1;
        
        % Search for 2P Times file; legacy version 
        twop_filepath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, '2P');
        if isempty(twop_filepath) % ucl open version 
            twop_filepath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'MatrixArduino');
        end
        
        if isempty(twop_filepath)
            error('No 2P or MatrixArduino file found for %s', stimName);
        end
        
        planeTimes_table = get_bonsai_twopframetimes_by_planes(twop_filepath, numPlanes);
        twoPData = struct(); 
        
        
        for iPlane = 1:sessionFileInfo.numPlanes
            planeName = sessionFileInfo.suite2pFiles(iPlane).planeName;
            fAll_filepath = findFile(sessionFileInfo.suite2pFiles(iPlane).planes, 'fall');
            fAll = load(fAll_filepath);
            
            % Get frame run indices
            if strcmp(planeName, 'plane_z')
                frameRun = sessionFileInfo.stim_framerun{['plane' '0'], stimName}{1}; %take the first one; does not matter they will all be the same 
                referencePlaneInfo = fAll.ops.planes_across_time(frameRun(1):(frameRun(1) + trimLength - 1));% was previoursly called currentplane
                
            else 
                frameRun = sessionFileInfo.stim_framerun{planeName, stimName}{1};
            end
            trimIndices = frameRun(1):(frameRun(1) + trimLength - 1);

            fullBadFrames = fAll.ops.badframes;
            
            % Populate struct with Suite2p data
            twoPData(iPlane).badframes = fullBadFrames(trimIndices)'; % this is new 
            twoPData(iPlane).planeName = planeName;
            twoPData(iPlane).ops = fAll.ops;
            twoPData(iPlane).iscell = fAll.iscell;
            twoPData(iPlane).stat = fAll.stat;
            twoPData(iPlane).redcell = fAll.redcell;
            twoPData(iPlane).frameRun = frameRun;
            twoPData(iPlane).F = fAll.F(:, trimIndices);
            twoPData(iPlane).spks = fAll.spks(:, trimIndices);
            twoPData(iPlane).Fneu = fAll.Fneu(:, trimIndices);
            
            % Extract timing data per plane
            if strcmp(planeName,'plane_z')
                planeTimes_trim = planeTimes_table.plane0(1:trimLength,:);
                for iTime = 1:trimLength
                    planeTimes_trim(iTime,:) = planeTimes_table.(['plane' num2str(referencePlaneInfo(iTime))])(iTime, :);
                end
            else
                planeTimes_trim = planeTimes_table.(planeName)(1:trimLength, :);
            end
            
            % Map timing fields 
            twoPData(iPlane).TwoPFrameTime = planeTimes_trim.TwoPFrameTime;
            twoPData(iPlane).ArduinoTime = planeTimes_trim.ArduinoTime;
            twoPData(iPlane).LastSyncPulseTime = planeTimes_trim.LastSyncPulseTime;
            
            % Safe check for missing fields that were previously saved
            if ismember('BonsaiTime', planeTimes_trim.Properties.VariableNames)
                twoPData(iPlane).BonsaiTime = planeTimes_trim.BonsaiTime;
            else
                twoPData(iPlane).BonsaiTime = NaN(height(planeTimes_trim), 1);
            end
            
            if ismember('RenderFrameCount', planeTimes_trim.Properties.VariableNames)
                twoPData(iPlane).RenderFrameCount = planeTimes_trim.RenderFrameCount;
            else
                twoPData(iPlane).RenderFrameCount = NaN(height(planeTimes_trim), 1);
            end
        end
        
        save(merged_data_filepath, 'twoPData');
        sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData = merged_data_filepath;
        
    catch ME
        warning('Failed to process stimulus "%s". Error: %s', stimName, ME.message);
        sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData = []; 
    end
end
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end

% function [sessionFileInfo] = mergeBonsaiSuite2pFiles(sessionFileInfo)
% % Merges Bonsai and Suite2p data for each stimulus in a session.
% % This version is robust to errors, continuing to the next stimulus
% % even if one fails to process. It saves all successful filepaths.
% %
% % Sonali and Aman - Jan 2025
% % Modified for robustness - Oct 2025
% 
% for iStim = 1:length(sessionFileInfo.stimFiles)
% 
%     stimName = sessionFileInfo.stimFiles(iStim).name;
% 
%     % --- Wrap the processing for each stimulus in a try...catch block ---
%     try
%         fprintf('Processing mergeBonsaiSuite2PFiles for %s\n', stimName);
% 
%         % Define the output filepath for the merged .mat file
%         stimFileName = [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_2pData' '_' stimName '.mat'];
%         merged_data_filepath = fullfile(sessionFileInfo.Directories.save_folder, stimFileName);
% 
%         %  Determine the trim length based on the last plane's frame run
%         if strcmp(sessionFileInfo.suite2pFiles(end).planeName, 'plane')
%             lastPlane_FrameRun = sessionFileInfo.stim_framerun{['plane' '0'], stimName}{1};
%             numPlanes = sessionFileInfo.origNPlanes;
%         else
%             lastPlane_FrameRun = sessionFileInfo.stim_framerun{sessionFileInfo.suite2pFiles(end).planeName, stimName}{1};
%             numPlanes = sessionFileInfo.numPlanes;
%         end
%         trimLength = lastPlane_FrameRun(2) - lastPlane_FrameRun(1) + 1;
% 
%         % Load the Bonsai 2p plane times
%         twop_filepath = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, '2P');
% 
%         planeTimes_table = get_bonsai_twopframetimes_by_planes(twop_filepath, numPlanes);
% 
%         % Initialize a temporary struct for this stimulus's data
%         twoPData = struct(); 
% 
%         %% Loop through each plane for the current stimulus
%         for iPlane = 1:sessionFileInfo.numPlanes
%             planeName = sessionFileInfo.suite2pFiles(iPlane).planeName;
% 
%             fAll_filepath = findFile(sessionFileInfo.suite2pFiles(iPlane).planes, 'fall');
%             if isempty(fAll_filepath)
%                warning('fall.mat not found for %s. Skipping this plane.', planeName);
%                continue;
%             end
%             fAll = load(fAll_filepath);
% 
%             % Get stimulus-specific frame run and trim indices
%             if strcmp(planeName, 'plane')
%                 frameRun = sessionFileInfo.stim_framerun{['plane' '0'], stimName}{1};
%                 currentPlaneInfo = fAll.ops.current_plane(frameRun(1):(frameRun(1) + trimLength - 1));
%             else
%                 frameRun = sessionFileInfo.stim_framerun{planeName, stimName}{1};
%             end
%             trimIndices = frameRun(1):(frameRun(1) + trimLength - 1);
% 
%             % Populate the twoPData struct for this plane
%             twoPData(iPlane).planeName = planeName;
%             twoPData(iPlane).ops = fAll.ops;
%             twoPData(iPlane).iscell = fAll.iscell;
%             twoPData(iPlane).stat = fAll.stat;
%             twoPData(iPlane).redcell = fAll.redcell;
%             twoPData(iPlane).frameRun = frameRun;
%             twoPData(iPlane).F = fAll.F(:, trimIndices);
%             twoPData(iPlane).spks = fAll.spks(:, trimIndices);
%             twoPData(iPlane).Fneu = fAll.Fneu(:, trimIndices);
% 
%             % Trim and add Bonsai data
%             if strcmp(planeName,'plane')
%                 planeTimes_trim = planeTimes_table.plane0(1:trimLength,:);
%                 for iTime = 1:trimLength
%                     planeTimes_trim(iTime,:) = planeTimes_table.([planeName num2str(currentPlaneInfo(iTime))])(iTime, :);
%                 end
%             else
%                 planeTimes_trim = planeTimes_table.(planeName)(1:trimLength, :);
%             end
%             twoPData(iPlane).TwoPFrameTime = planeTimes_trim.TwoPFrameTime;
%             twoPData(iPlane).BonsaiTime = planeTimes_trim.BonsaiTime;
%             twoPData(iPlane).ArduinoTime = planeTimes_trim.ArduinoTime; % why was this missing.. 
%         end
% 
%         % If all planes processed successfully, save the data
%         save(merged_data_filepath, 'twoPData');
% 
%         % Crucially, save the successful filepath to sessionFileInfo
%         sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData = merged_data_filepath;
% 
%     catch ME
%         %  This block runs ONLY if an error occurred in the 'try' block ---
%         warning('Failed to process stimulus "%s". Error: %s', stimName, ME.message);
%         fprintf(2, '    -> Error occurred in function %s at line %d.\n', ME.stack(1).name, ME.stack(1).line);
% 
%         % Set the filepath to empty to indicate failure
%         sessionFileInfo.stimFiles(iStim).mergedBonsai2PSuite2pData = []; 
%     end
% end
% 
% % --- Save the final sessionFileInfo ONCE after the loop is finished ---
% disp('All stimuli processed. Saving final sessionFileInfo...');
% save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
% disp('Done.');
% end