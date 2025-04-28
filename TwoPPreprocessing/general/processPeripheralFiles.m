function [sessionFileInfo] = processPeripheralFiles(sessionFileInfo)
% Process the peripheral files and save the processed info into a
% directory; Only include the periphereals that are common 

% Aman and Sonali - Dec 2024

for iStim = 1:length(sessionFileInfo.stimFiles)
    %% Photodiode file
    peripheralData = struct() ;
    sessionFileInfo.stimFiles(iStim).processedPeripheralData = fullfile(sessionFileInfo.Directories.save_folder,...
        [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_PeripheralData' '_' sessionFileInfo.stimFiles(iStim).name '.mat']);

    photodiode_path     = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Photodiode');
    photodiode_table    = readtable(photodiode_path);
    if ismember('ArduinoTime', photodiode_table.Properties.VariableNames) && ismember('PDOutput', photodiode_table.Properties.VariableNames)
        rows_to_remove = (photodiode_table.ArduinoTime == 0) | (photodiode_table.PDOutput == 0);
        photodiode_table(rows_to_remove, :) = [];
    end
    
    % to remove repeated time measures
    [~,keep_idx,~] = unique(photodiode_table.ArduinoTime);

    peripheralData.Photodiode.rawArduinoTime  = photodiode_table.ArduinoTime(keep_idx)./1000;
    peripheralData.Photodiode.rawBonsaiTime   = photodiode_table.BonsaiTime(keep_idx);
    peripheralData.Photodiode.rawValue        = photodiode_table.PDOutput(keep_idx);
    peripheralData.Photodiode.rawRenderFrameCount = photodiode_table.RenderFrameCount(keep_idx); 
    peripheralData.Photodiode.rawLastSyncPulseTime = photodiode_table.LastSyncPulseTime(keep_idx);
    
    %% Quad File
    quad_path = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Quad');
    if exist(quad_path, 'file')
        fprintf('Found Quad file: %s\n', quad_path);  % Debugging
        
        quadstate_table = readtable(quad_path);
        if ismember('ArduinoTime', quadstate_table.Properties.VariableNames)
            rows_to_remove = (quadstate_table.ArduinoTime == 0);
            quadstate_table(rows_to_remove, :) = [];
        end
    
        % Remove repeated time measures
        [~, keep_idx, ~] = unique(quadstate_table.ArduinoTime);
    
        peripheralData.Quadstate.rawArduinoTime  = quadstate_table.ArduinoTime(keep_idx)./1000;
        peripheralData.Quadstate.rawBonsaiTime   = quadstate_table.BonsaiTime(keep_idx);
        peripheralData.Quadstate.rawRenderFrameCount = quadstate_table.RenderFrameCount(keep_idx); 
        peripheralData.Quadstate.rawLastSyncPulseTime = quadstate_table.LastSyncPulseTime(keep_idx);
        peripheralData.Quadstate.rawValue = quadstate_table.QuadState(keep_idx);
        
    else
        fprintf('Quad file missing for stimulus: %s\n', sessionFileInfo.stimFiles(iStim).name);
        peripheralData.quadstate = [];
    end

    %% Wheel File
    wheel_path      = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Wheel');
    wheel_table     = readtable(wheel_path);% Read the Bonsai wheel table

    % Remove rows where ArduinoTime or Wheel are 0
    if ismember('ArduinoTime', wheel_table.Properties.VariableNames) && ismember('Wheel', wheel_table.Properties.VariableNames)
        rows_to_remove = (wheel_table.ArduinoTime == 0) | (wheel_table.Wheel == 0);
        wheel_table(rows_to_remove, :) = [];
    end
    % wheel_table = wheel_table(wheel_table.Wheel ~= 0, :);
    % Extract the wheel iutput (raw input) and corresponding wheel time
    % to remove repeated time measures
    [~,keep_idx,~] = unique(wheel_table.ArduinoTime);

    peripheralData.Wheel.rawArduinoTime  = wheel_table.ArduinoTime(keep_idx)./1000;
    peripheralData.Wheel.rawBonsaiTime   = wheel_table.BonsaiTime(keep_idx);
    peripheralData.Wheel.rawValue        = wheel_table.Wheel(keep_idx);
    peripheralData.Wheel.rawRenderFrameCount = wheel_table.RenderFrameCount(keep_idx); 
    peripheralData.Wheel.rawLastSyncPulseTime = wheel_table.LastSyncPulseTime(keep_idx);

    save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData")
    save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
end