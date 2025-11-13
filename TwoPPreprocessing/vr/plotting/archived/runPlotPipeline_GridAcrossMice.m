
%% Define days of learning 
% Define all the days of experience you want to see on the grid
days_to_plot = [1, 2, 3, 4, 5]; 

% Define any mice you want to skip
mice_to_skip = {'M24046', 'M24043', 'M24048', 'M24049'};



disp('Creating the master table for plotting...');
table_for_plotting = filterMasterTable(...
    'Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'DayOfExperience', days_to_plot, ...  
    'ImagedType', 'Boutons');

disp('Table created successfully.');
disp(table_for_plotting); % Check it


saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures';
fileName = 'All_Mice_Tuning_Grid_DaysOfLearning_OddSortedByEven-Medians.png';
fullSavePath = fullfile(saveDir, fileName);

% Call the function. You give it the big table, the exclude list,
% and now the full path where you want it saved.
disp('Calling the grid plotting function...');
plotPositionTuning_GridAcrossMice(table_for_plotting, ...
    'ExcludeMice', mice_to_skip, ...
    'SavePath', fullSavePath); 
disp(['Plotting complete. Figure saved to: ' fullSavePath]);