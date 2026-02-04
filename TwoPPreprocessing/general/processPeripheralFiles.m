function [sessionFileInfo] = processPeripheralFiles(sessionFileInfo)
% Process the peripheral files and save the processed info into a directory.
% Aman and Sonali - Dec 2024
% Modified by Sonali Jan 2026 to integrate UCL-Open logging structure.
%
% NOTE TO SELF: 
% We are explicitly splitting and saving peripheral data into individual 
% .mat files here because the rest of the downstream pipeline expects 
% this specific directory/struct architecture for synchronization and 
% plotting. Do not consolidate (yet) unless updating the entire pipeline!

for iStim = 1:length(sessionFileInfo.stimFiles)
    %% Initialization
    peripheralData = struct();
    sessionFileInfo.stimFiles(iStim).processedPeripheralData = fullfile(sessionFileInfo.Directories.save_folder,...
        [sessionFileInfo.animal_name '_' sessionFileInfo.session_name '_PeripheralData' '_' sessionFileInfo.stimFiles(iStim).name '.mat']);
    
    % Find file paths
    matrixArduino_path = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'MatrixArduino');
    photodiode_path    = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Photodiode');
    quad_path          = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Quad');
    wheel_path         = findFile(sessionFileInfo.stimFiles(iStim).bonsai_filepaths, 'Wheel');

    %% 1. Photodiode processing
    if exist(photodiode_path, 'file')
        % Preserve naming to prevent warnings about dots in headers
        photodiode_table = readtable(photodiode_path, 'VariableNamingRule', 'preserve');
        
        if ismember('ArduinoTime', photodiode_table.Properties.VariableNames)
            rows_to_remove = (photodiode_table.ArduinoTime == 0) | (photodiode_table.PDOutput == 0);
            photodiode_table(rows_to_remove, :) = [];
            
            [~,keep_idx,~] = unique(photodiode_table.ArduinoTime);
            peripheralData.Photodiode.rawArduinoTime       = photodiode_table.ArduinoTime(keep_idx)./1000;
            peripheralData.Photodiode.rawBonsaiTime        = photodiode_table.BonsaiTime(keep_idx);
            peripheralData.Photodiode.rawValue             = photodiode_table.PDOutput(keep_idx);
            peripheralData.Photodiode.rawRenderFrameCount  = photodiode_table.RenderFrameCount(keep_idx); 
            peripheralData.Photodiode.rawLastSyncPulseTime = photodiode_table.LastSyncPulseTime(keep_idx)./1000; % Has always been saved in ms 
        end
    
    elseif exist(matrixArduino_path, 'file')
        matrixArduino_table = readtable(matrixArduino_path, 'VariableNamingRule', 'preserve');
        
        if ismember('Seconds', matrixArduino_table.Properties.VariableNames)
            rows_to_remove = (matrixArduino_table.Seconds == 0) | (matrixArduino_table.("Value.PhotodiodeValue") == 0);
            matrixArduino_table(rows_to_remove, :) = [];
            
            [~,keep_idx,~] = unique(matrixArduino_table.Seconds);
            peripheralData.Photodiode.rawArduinoTime       = matrixArduino_table.Seconds(keep_idx); 
            peripheralData.Photodiode.rawValue             = matrixArduino_table.("Value.PhotodiodeValue")(keep_idx);
            peripheralData.Photodiode.rawLastSyncPulseTime = matrixArduino_table.("Value.LastSyncPulseTime")(keep_idx)./1000; % saved in ms 
        end
    else
        peripheralData.Photodiode = [];
    end   

    %% 2. QuadState Processing
    if exist(quad_path, 'file')
        quadstate_table = readtable(quad_path, 'VariableNamingRule', 'preserve');
        if ismember('ArduinoTime', quadstate_table.Properties.VariableNames)
            rows_to_remove = (quadstate_table.ArduinoTime == 0);
            quadstate_table(rows_to_remove, :) = [];
            [~, keep_idx, ~] = unique(quadstate_table.ArduinoTime);
    
            peripheralData.Quadstate.rawArduinoTime       = quadstate_table.ArduinoTime(keep_idx)./1000;
            peripheralData.Quadstate.rawBonsaiTime        = quadstate_table.BonsaiTime(keep_idx);
            peripheralData.Quadstate.rawRenderFrameCount  = quadstate_table.RenderFrameCount(keep_idx); 
            peripheralData.Quadstate.rawLastSyncPulseTime = quadstate_table.LastSyncPulseTime(keep_idx);
            peripheralData.Quadstate.rawValue             = quadstate_table.QuadState(keep_idx);
        
        elseif ismember('Seconds', quadstate_table.Properties.VariableNames)
            rows_to_remove = (quadstate_table.Seconds == 0);
            quadstate_table(rows_to_remove, :) = [];
            [~, keep_idx, ~] = unique(quadstate_table.Seconds);
            
            peripheralData.Quadstate.rawArduinoTime = quadstate_table.Seconds(keep_idx); 
            peripheralData.Quadstate.rawValue       = quadstate_table.Value(keep_idx);
        end
    else
        peripheralData.Quadstate = [];
    end

    %% 3. Wheel Processing
    if exist(wheel_path, 'file')
        wheel_table = readtable(wheel_path, 'VariableNamingRule', 'preserve');
        if ismember('ArduinoTime', wheel_table.Properties.VariableNames)
            rows_to_remove = (wheel_table.ArduinoTime == 0);
            wheel_table(rows_to_remove, :) = [];
            
            [~,keep_idx,~] = unique(wheel_table.ArduinoTime);
            peripheralData.Wheel.rawArduinoTime       = wheel_table.ArduinoTime(keep_idx)./1000;
            peripheralData.Wheel.rawBonsaiTime        = wheel_table.BonsaiTime(keep_idx);
            peripheralData.Wheel.rawValue             = wheel_table.Wheel(keep_idx);
            peripheralData.Wheel.rawRenderFrameCount  = wheel_table.RenderFrameCount(keep_idx); 
            peripheralData.Wheel.rawLastSyncPulseTime = wheel_table.LastSyncPulseTime(keep_idx);
        end

    elseif exist(matrixArduino_path, 'file')
        matrixArduino_table = readtable(matrixArduino_path, 'VariableNamingRule', 'preserve');
        
        if ismember('Seconds', matrixArduino_table.Properties.VariableNames) && ...
           ismember('Value.EncoderPos', matrixArduino_table.Properties.VariableNames)
            
            rows_to_remove = (matrixArduino_table.Seconds == 0);
            matrixArduino_table(rows_to_remove, :) = [];
            
            [~,keep_idx,~] = unique(matrixArduino_table.Seconds);
            peripheralData.Wheel.rawArduinoTime       = matrixArduino_table.Seconds(keep_idx); 
            peripheralData.Wheel.rawValue             = matrixArduino_table.("Value.EncoderPos")(keep_idx);
            
            if ismember('Value.LastSyncPulseTime', matrixArduino_table.Properties.VariableNames)
                peripheralData.Wheel.rawLastSyncPulseTime = matrixArduino_table.("Value.LastSyncPulseTime")(keep_idx);
            end
        end
    else
        peripheralData.Wheel = [];
    end
    
    %% Save individual stim peripheral data
    save(sessionFileInfo.stimFiles(iStim).processedPeripheralData, "peripheralData")
    fprintf('Saved Processed Data: %s\n', sessionFileInfo.stimFiles(iStim).name);
end

% Final update to sessionFileInfo
save(sessionFileInfo.sessionFileInfo_filepath, 'sessionFileInfo');
disp('Full peripheral processing complete.');

end