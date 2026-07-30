%% Plot learning FOVs 
% Load RSP sessions 
RSPSessions = filterMasterTable('MouseID', {'M25132', 'M25133', 'M26003'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', [1 2 3 4 5]);


uniqueMice = unique(RSPSessions.MouseID);
numMice = length(uniqueMice);
maxDays = max(RSPSessions.DayOfExperience);

fileName = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\FOVs\FOV_ComparisonRSP.png'; 

% Create a figure with a grid: Rows = Mice, Cols = Days
figure('Name', 'FOV refImage', 'Color', 'w');
tiledlayout(numMice, maxDays, 'TileSpacing', 'compact', 'Padding', 'tight');

for m = 1:numMice
    currMouse = uniqueMice(m);
    
    for d = 1:maxDays
        nexttile;
        
        % Find the row in RSPSessions for this mouse and day
        sessionRow = RSPSessions(RSPSessions.MouseID == currMouse & ...
                                 RSPSessions.DayOfExperience == d, :);
        
        if isempty(sessionRow)
            axis off; title(sprintf('%s: Day %d (N/A)', currMouse, d));
            continue;
        end
        
        % Extract Session info 
        mouseID = char(sessionRow.MouseID);
        session = char(sessionRow.Session);
        
        % Load sessionFileInfo
        sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID,session);
        load(sessionFileInfoFilePath); 
        
        try
            % Find the index for 'Darkness'
            stimNames = {sessionFileInfo.stimFiles.name};
            darknessIdx = find(contains(stimNames, 'RFMapping', 'IgnoreCase', true), 1);
            
            if isempty(darknessIdx)
                error('No Darkness stim found');
            end
            
            % Load the 2pData file
            dataPath = sessionFileInfo.stimFiles(darknessIdx).mergedBonsai2PSuite2pData;
            data = load(dataPath);
            
            % Plot the max projection
            if isfield(data.twoPData.ops, 'refImg')
                imagesc(data.twoPData.ops.refImg);
                colormap gray;
                axis image off;
                title(sprintf('%s Day %d', mouseID, d));
            else
                text(0.5, 0.5, 'No max\_proj', 'HAlign', 'center');
                axis off;
            end
            
        catch ME
            fprintf('Could not load Day %d for %s: %s\n', d, mouseID, ME.message);
            axis off;
            title(sprintf('%s Day %d (Error)', mouseID, d));
        end
    end
end

exportgraphics(gcf, fileName, 'Resolution', 600);
saveFigureFormats(fig_grid, fullfile(suite2p_dir, 'all_plane_mean_images'));
