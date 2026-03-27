function taskParams = parseTaskLogic(jsonPath)
    % Reads and decodes the specific UclOpenVrCorridor2p JSON file
    
    str = fileread(jsonPath);
    data = jsondecode(str);
    
    trials = data.taskParameters.blocks.availableTrials;
    

    taskParams.movementVisualGain = trials(1).movementVisualGain;
    
    % Mapping 'inner' and 'outer' bounds:
    % In these VR tasks, boundaryThreshold is often where the trial starts/triggers
    % and endTrialThreshold is where the corridor ends.
    taskParams.innerBound = trials(1).boundaryThreshold;
    taskParams.outerBound = trials(1).endTrialThreshold;
    
    % Additional useful metadata
    taskParams.corridorWidth = data.taskParameters.corridorWidth;
    taskParams.maxTrials = data.taskParameters.blocks.maxTrials;
    
    fprintf('Parsed Task: Gain=%.4f', taskParams.movementVisualGain);
end