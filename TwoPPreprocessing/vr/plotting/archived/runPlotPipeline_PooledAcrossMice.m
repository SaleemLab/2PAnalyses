
% days_to_plot = [1, 2, 3, 4, 5]; 
days_to_plot = 'PostLearning';

%% Create TWO filtered table
disp('Creating the master table for plotting...');
tableForPlotting_RSP = filterMasterTable(...
    'Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'DayOfExperience', days_to_plot, ...
    'TargetArea', 'RSP');

disp('Table created successfully.');

tableForPlotting_Entm = filterMasterTable(...
    'Exclude', 0, ...
    'Suite2PPreprocessing', 1, ...
    'DayOfExperience', days_to_plot, ...
    'TargetArea', 'ENTm');

% tableForPlotting_V1 = filterMasterTable(...
%     'Exclude', 0, ...
%     'Suite2PPreprocessing', 1, ...
%     'DayOfExperience', 'PostLearning', ...
%     'TargetArea', 'V1');

%% Define save directory
saveDir = 'Z:\ibn-vision\USERS\Sonali\Figures';

%% RSP
disp('Plotting RSP Population...');
rspSavePath = fullfile(saveDir, 'RSP_Population_Tuning_OddLapsSortedByEven.png');

plotPopulationTuning_PooledAcrossMice(tableForPlotting_RSP, 'RSP', ...
    'SavePath', rspSavePath);

%% Entm
disp('Plotting ENTm Population...');
entmSavePath = fullfile(saveDir, 'ENTm_Population_Tuning_OddLapsSortedByEven.png');

plotPopulationTuning_PooledAcrossMice(tableForPlotting_Entm, 'ENTm', ...
    'SavePath', entmSavePath);

disp('All population plots saved.');

%% V1 
disp('Plotting ENTm Population...');
entmSavePath = fullfile(saveDir, 'V1_Population_Tuning_OddLapsSortedByEven.png');

plotPopulationTuning_PooledAcrossMice(tableForPlotting_Entm, 'V1', ...
    'SavePath', entmSavePath);

disp('All population plots saved.');
