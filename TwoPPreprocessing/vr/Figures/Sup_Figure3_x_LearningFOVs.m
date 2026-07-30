%% Plot learning FOVs 
% Load RSP sessions 
RSPSessions = filterMasterTable('MouseID', {'M25132', 'M25133', 'M26003'}, 'Suite2PPreprocessing', 1, 'DayOfExperience', [1 2 3 4 5]);


uniqueMice = unique(RSPSessions.MouseID);
numMice = length(uniqueMice);
maxDays = max(RSPSessions.DayOfExperience);

fileName = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.4.1_Learning_Supp\allFOVs'; 

% Create a figure with a grid: Rows = Mice, Cols = Days
fig=figure('Name', 'FOV refImage', 'Color', 'w');
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


saveFigureFormats(fig, fullfile(fileName, 'all_days_fovs_mean_images'));
