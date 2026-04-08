function plotThresholdSensitivity(allData)
    % Define the target directory
    saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\DistibutionsAllCritera';
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end

    % 1. Filter for Day 4
    if isnumeric(allData(1).Day)
        dayIdx = [allData.Day] == 4;
    else
        dayIdx = str2double({allData.Day}) == 4;
    end
    data = allData(dayIdx);
    
    if isempty(data), error('No data found for Day 4'); end

    % Pre-calculate SigShuffle Mask (Must pass Peak OR Range)
    % This is our "Base Population"
    validIdx = cell(length(data), 1);
    totalSigCells = 0;
    for s = 1:length(data)
        validIdx{s} = (data(s).isSignificantByPeakShuffling == 1) | ...
                      (data(s).isSignificantByRange == 1);
        totalSigCells = totalSigCells + sum(validIdx{s});
    end

    % Define Threshold Ranges
    range_Corr = 0:0.02:0.9;             % For Halves and Odd-Even
    range_VarVar = 1:1:60;               % Var/TuningVar 
    range_VarRange = 0.1:0.1:5;          % Var/TuningRange 


    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.3 0.9 0.4]);
    t = tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'loose');
    title(t, sprintf('Day 4 Threshold Sensitivity (Starting N = %d SigShuffle ROIs; pooled across mice)', totalSigCells), ...
        'FontSize', 14, 'FontWeight', 'bold');

     % @gemini
    getCounts = @(field, thresholds, isGreater) ...
        arrayfun(@(thr) calculateRemaining(data, validIdx, field, thr, isGreater), thresholds);

    % Halves Correlation 
    nexttile;
    plot(range_Corr, getCounts('lapCorr_HalvesRho', range_Corr, true), 'LineWidth', 2, 'Color', [0.2 0.6 0.2]);
    hold on; xline(0.4, '--r', '0.4', 'LabelVerticalAlignment', 'bottom');
    xlabel('Rho'); ylabel('Remaining ROIs'); title('Halves Correlation'); grid on;

    % Odd-Even Correlation
    nexttile;
    plot(range_Corr, getCounts('lapCorr_OddEvenRho', range_Corr, true), 'LineWidth', 2, 'Color', [0.2 0.5 0.5]);
    hold on; xline(0.4, '--r', '0.4', 'LabelVerticalAlignment', 'bottom');
    xlabel('Rho'); title('Odd-Even Correlation'); grid on;

    %Ratio Var to Tuning Var
    nexttile;
    plot(range_VarVar, getCounts('ratioVarToTuningVar', range_VarVar, false), 'LineWidth', 2, 'Color', [0.6 0.2 0.2]);
    hold on; xline(20, '--r', '20', 'LabelVerticalAlignment', 'bottom');
    xlabel('Ratio'); title('Var / TuningVar'); grid on;

    % Ratio Var to Tuning Range
    nexttile;
    plot(range_VarRange, getCounts('ratioVarToTuningRange', range_VarRange, false), 'LineWidth', 2, 'Color', [0.2 0.2 0.6]);
    hold on; xline(1.0, '--r', '1.0', 'LabelVerticalAlignment', 'bottom');
    xlabel('Ratio'); title('Var / TuningRange'); grid on;

    % Save result
    saveas(fig, fullfile(saveDir, 'Threshold_Sensitivity_Day4.png'));
end

% Sub-function to handle the counting loop safely 
function count = calculateRemaining(data, validIdx, field, thr, isGreater)
    count = 0;
    for s = 1:length(data)
        vals = data(s).(field);
        mask = validIdx{s};
        if isGreater
            count = count + sum(vals(mask) >= thr);
        else
            count = count + sum(vals(mask) <= thr);
        end
    end
end