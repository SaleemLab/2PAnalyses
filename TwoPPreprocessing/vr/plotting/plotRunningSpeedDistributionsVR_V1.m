function plotRunningSpeedDistributionsVR_V1(MouseID)
filteredTable = filterMasterTable('Exclude', 0, 'Suite2PPreprocessing', 1, 'MouseID', MouseID);
mouseInfo = sessionsToProcess(filteredTable);
sessionNames = mouseInfo{1, 2};

% Pre-calculate total number of required plots (tiles) ---
allVrStimNames = {};
for thisSession = 1:length(sessionNames)
    sessionName = sessionNames{thisSession};
    try
        stimList = getStimList(MouseID, sessionName);
        
        vrStimNames = {};
        for thisStim = 1:length(stimList)
            if any(contains(stimList{thisStim}, 'VRCorr', 'IgnoreCase', true))
                vrStimNames{end+1} = stimList{thisStim};
            end
        end
        
        if contains(vrStimNames, 'CombinedRun', 'IgnoreCase',true)
            combinedRunIdx = find(contains(vrStimNames, 'CombinedRun', 'IgnoreCase', true), 1);
            vrStimNames = vrStimNames(combinedRunIdx);
        end
        
        allVrStimNames = [allVrStimNames, vrStimNames];
    catch
        continue;
    end
end

% Set up the tiled layout based on the total number of stimuli
numTotalPlots = length(allVrStimNames);
if numTotalPlots == 0
    fprintf('No VR stimuli found for mouse %s.\n', MouseID);
    return;
end

% Use sqrt to determine a roughly square layout (e.g., 9 plots = 3x3)
numCols = ceil(sqrt(numTotalPlots));
numRows = ceil(numTotalPlots / numCols);

figure;
tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');
fprintf('Processing Mouse: %s\n', MouseID);

% --- Plotting Loop ---
for thisSession = 1:length(sessionNames)
    sessionName = sessionNames{thisSession};
    
    fprintf('\n-- Processing Session: %s --\n', sessionName);
    try
        stimList = getStimList(MouseID, sessionName);
        
        vrStimNames = {};
        for thisStim = 1:length(stimList)
            if any(contains(stimList{thisStim}, 'VRCorr', 'IgnoreCase', true))
                vrStimNames{end+1} = stimList{thisStim};
            end
        end
        
        if contains(vrStimNames, 'CombinedRun', 'IgnoreCase',true)
            combinedRunIdx = find(contains(vrStimNames, 'CombinedRun', 'IgnoreCase', true), 1);
            vrStimNames = vrStimNames(combinedRunIdx);
        end 
        
        sessionFileInfo = get2PsessionFilePaths(MouseID, sessionName, stimList);
        
        if ~isempty(vrStimNames)
            fprintf('Found %d VR stimulus file(s).\n', length(vrStimNames));
            
            for thisVRStim = 1:length(vrStimNames)
                vrStimName = vrStimNames{thisVRStim};
                fprintf('Processing VR Stim: %s\n', vrStimName);
                
                stimIdx = find(strcmp(vrStimName, {sessionFileInfo.stimFiles.name}), 1); 
                if isempty(stimIdx), error('Specified VRStimName not found in sessionFileInfo.'); end
                
                load(sessionFileInfo.stimFiles(stimIdx).Response, 'response');
                sessionWheelSpeed = response.wheelSpeed; 
                % Call nexttile for EACH unique stimulus/plot
                nexttile; 
                percentageRun = round(100*sum(sessionWheelSpeed > 1) / length(sessionWheelSpeed));
                histogram(sessionWheelSpeed(sessionWheelSpeed > 1), 'Normalization', 'count', 'DisplayName', vrStimName, 'BinEdges', 1:5:60); 
                
                % Combine session and stim name in the title
                correctedVRStimName = replace(vrStimName, '_', '\_');
                title([correctedVRStimName ' ' num2str(percentageRun) '% Run']);
                
                ylabel('Count');
            end
        end
    catch ME
        fprintf('Error processing session %s: %s\n', sessionName, ME.message);
    end
end
sgtitle(sprintf('Running Speed Distributions for Mouse: %s', MouseID));
end