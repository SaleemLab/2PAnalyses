% Load RSP sessions 
RSPSessions = filterMasterTable('MouseID', {'M25132', 'M25133', 'M26003'}, 'Suite2PPreprocessing', 1, 'HasStimulus', {'LandManipCorridor', 'BaselineCorridor'});
uniqueMice = unique(RSPSessions.MouseID);
numMice = length(uniqueMice);

% Find the maximum number of sessions any one mouse has to set the grid width
sessionCounts = groupsummary(RSPSessions, 'MouseID');
maxSessions = max(sessionCounts.GroupCount);

fileName = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\FOVs\FOV_ComparisonRSP_BaselineAndManip.png'; 

% Create figure
figure('Name', 'FOV Sequential Sessions', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout(numMice, maxSessions, 'TileSpacing', 'compact', 'Padding', 'tight');

for m = 1:numMice
    currMouse = uniqueMice(m);
    
    % Get all rows for this specific mouse, sorted by Session
    mouseSessions = RSPSessions(RSPSessions.MouseID == currMouse, :);
    mouseSessions = sortrows(mouseSessions, 'Session'); 
    
    numMouseSessions = height(mouseSessions);
    
    for s = 1:maxSessions
        nexttile;
        
        if s <= numMouseSessions
            % Extract data
            mouseID = char(mouseSessions.MouseID(s));
            sessionStr = char(mouseSessions.Session(s)); % e.g., '20260219'
            dayVal  = mouseSessions.DayOfExperience(s);
            
            % Load and find data (same as before)
            sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID, sessionStr);
            load(sessionFileInfoFilePath); 
            
            try
                stimNames = {sessionFileInfo.stimFiles.name};
                stimIdx = find(contains(stimNames, 'RFMapping', 'IgnoreCase', true), 1);
                
                dataPath = sessionFileInfo.stimFiles(stimIdx).mergedBonsai2PSuite2pData;
                data = load(dataPath);
                
                if isfield(data.twoPData.ops, 'refImg')
                    imagesc(data.twoPData.ops.refImg);
                    colormap gray;
                    axis image off;
                    
                    % UPDATED TITLE: Includes Session ID
                    title({sprintf('%s | %s', mouseID, sessionStr); ...
                           sprintf('S%d (Day %d)', s, dayVal)}, 'FontSize', 8);
                else
                    text(0.5, 0.5, 'No refImg', 'HAlign', 'center');
                    axis off;
                end
                
            catch ME
                axis off;
                title(sprintf('Error %s', sessionStr));
            end
        else
            axis off;
        end
    end
end

% Save high-res
exportgraphics(gcf, fileName, 'Resolution', 600);