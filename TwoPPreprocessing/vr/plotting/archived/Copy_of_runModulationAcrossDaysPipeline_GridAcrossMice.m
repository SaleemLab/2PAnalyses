%% Define days of learning 
% Define all the days of experience you want to see on the grid
daysToPlot = [1, 2, 3, 4, 5]; 

% Define any mice you want to skip
miceToSkip = {'M24046', 'M24043', 'M24048', 'M24049'};



disp('Creating the master table for plotting...');
tableForPlotting = filterMasterTable(...
    'Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'DayOfExperience', daysToPlot, ...  
    'ImagedType', 'Boutons');

disp('Table created successfully.');
disp(tableForPlotting); % Check it


saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures';
fileName = 'ModulatioAcrossSpace_AcorssMiceAcrossDaysOfLearning_CDF_means.png';
fullSavePath = fullfile(saveDir, fileName);

% Call the function. You give it the big table, the exclude list,
% and now the full path where you want it saved.
disp('Calling the grid plotting function...');
plotModulationCDF_SuperimposedDays(tableForPlotting, ...
    'ExcludeMice', miceToSkip, ...
    'SavePath', fullSavePath); 
disp(['Plotting complete. Figure saved to: ' fullSavePath]);