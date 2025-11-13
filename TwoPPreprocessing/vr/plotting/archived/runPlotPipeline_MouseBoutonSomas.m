
days_to_plot = [1, 2, 3, 4, 5]; 
mice_to_plot = {'M25041', 'M25040', 'M25037'};

%% Create ONE filtered table
% Create ONE table that contains all data (boutons and somas) for these mice
disp('Creating the master table for plotting...');
full_table_for_plotting = filterMasterTable(...
    'Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'ImagedType', {'Boutons', 'Somas'}, ...
    'MouseID', mice_to_plot);

disp('Table created successfully.');
disp(head(full_table_for_plotting)); 

%% Save directory
saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures';

%% Loop through mice and plot
disp('Looping through mice to generate figures...');

for i = 1:length(mice_to_plot)
    mouseID = mice_to_plot{i};
    
    % Create a unique save path for this mouse
    fileName = sprintf('%s_Bouton_Soma_Tuning.png', mouseID);
    savePath = fullfile(saveDir, fileName);
    
    fprintf('Plotting figure for %s...\n', mouseID);
    
    % Call function for this mouse
    plotMouseBoutonAndSomaTuning(full_table_for_plotting, mouseID, ...
        'ImagedTypeColumn', 'TypeImaged', ... % Tell it the column name
        'SavePath', savePath);
        
end

disp('All mouse-specific plots saved.');
