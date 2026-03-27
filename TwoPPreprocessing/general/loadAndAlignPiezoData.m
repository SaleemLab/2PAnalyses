function alignedData = loadAndAlignPiezoData(pzFilePath, twoPFrameTimesFilePath)
    % Load raw data
    pzData = readtable(pzFilePath); 
    twoPData = readtable(twoPFrameTimesFilePath);

    % Find the first sync pulse recorded on both devices
    idxPz = find(pzData.LastSyncPulsetime > 0, 1, 'first');
    idx2P = find(twoPData.LastSyncPulseTime > 0, 1, 'first');

    if isempty(idxPz) || isempty(idx2P)
        error('Sync pulse not found. Check Arduino logs for trigger signals.');
    end

    % Define the Anchor points (Internal Clock values at the moment of pulse)
    tPulsePz = pzData.LastSyncPulsetime(idxPz);
    tPulse2P = twoPData.LastSyncPulseTime(idx2P);

    % Shift the Piezo timeline to the 2P clock
    % pzTimeIn2PClock = PiezoTime - (Time of Pulse on Pz - Time of Pulse on 2P)
    pzTimeIn2PClock = pzData.ArduinoTime - (tPulsePz - tPulse2P);

    % Create the output table (preserving all 1000Hz Piezo samples)
    alignedData = table();
    alignedData.Time_s = pzTimeIn2PClock / 1000; % Global time in seconds
    alignedData.PiezoValue = pzData.Piezo;

    % Map frame numbers to the Piezo samples using 'nearest' 
    % This tells you which frame was being captured at this specific millisecond
    alignedData.FrameNum = interp1(twoPData.TwoPFrameTime, (1:height(twoPData))', ...
                                   pzTimeIn2PClock, 'nearest', 'extrap');

    % Quick check: Set FrameNum to NaN for Piezo data recorded before/after the 2P session
    outOfBounds = pzTimeIn2PClock < min(twoPData.TwoPFrameTime) | ...
                  pzTimeIn2PClock > max(twoPData.TwoPFrameTime);
    alignedData.FrameNum(outOfBounds) = NaN;

    % Visualization
    figure('Color', 'w');
    plot(alignedData.Time_s, alignedData.PiezoValue, 'Color', [0.2 0.2 0.2]);
    hold on;
    xline(tPulse2P/1000, '--r', 'Sync Pulse', 'LineWidth', 1.5);
    xlabel('Time relative to 2P Arduino Boot (s)');
    ylabel('Piezo Signal (Raw)');
    title('1000Hz Piezo Signal Aligned to 2P Timeline');
    grid on;
end


