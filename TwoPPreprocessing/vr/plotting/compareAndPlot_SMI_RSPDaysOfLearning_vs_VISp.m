function compareAndPlot_SMI_RSPDaysOfLearning_vs_VISp(RSPData, VISpData)

    daysOfInterest = [1,2, 3,4, 5];  % day 200 merged into day 5

    % RSP: black to light grey across days
    dayColors_RSP = [...
        0.84, 0.19, 0.15;   % day 1 - red
        0.99, 0.55, 0.24;   % day 2 - orange
        0.13, 0.47, 0.71;   % day 3 - blue
        0.17, 0.63, 0.17;   % day 4 - green
        0.58, 0.15, 0.68];  % day 5/200 - purple

    % VISp: single comparison color (dashed)
    vispColor = [0.5, 0.5, 0.5];  % muted red to distinguish from RSP greys

    %% --- Pool SMI per day for RSP ---
    SMI_RSP_byDay       = cell(length(daysOfInterest), 1);
    daysPresent_RSP     = cell(length(daysOfInterest), 1);

    for s = 1:length(RSPData)
        sess = RSPData(s);
        if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMI   = sess.SMI.SMI;
            cleanSMI = rawSMI(sess.FilteredROIs);
            cleanSMI = cleanSMI(~isnan(cleanSMI));

            % RSP may not have dayofexperience — treat missing or 200 as day 5



            if isfield(sess, 'Day')
                day = sess.Day;
                dayIdx = find(daysOfInterest == day, 1);
                if ~isempty(dayIdx)
                    SMI_RSP_byDay{dayIdx}   = [SMI_RSP_byDay{dayIdx}; cleanSMI(:)];
                    daysPresent_RSP{dayIdx} = unique([daysPresent_RSP{dayIdx}, sess.Day]);
                end
            end
        end
    end

    %% --- Pool ALL SMI for VISp (single comparison line) ---
    pooledSMI_VISp = [];
    for s = 1:length(VISpData)
        sess = VISpData(s);
        if isfield(sess, 'SMI') && isfield(sess, 'FilteredROIs') && ~isempty(sess.FilteredROIs)
            rawSMI   = sess.SMI.SMI;
            cleanSMI = rawSMI(sess.FilteredROIs);
            cleanSMI = cleanSMI(~isnan(cleanSMI));
            pooledSMI_VISp = [pooledSMI_VISp; cleanSMI(:)];
        end
    end

    %% --- Print stats ---
    fprintf('\n--- RSP SMI by Day ---\n');
    for d = 1:length(daysOfInterest)
        rsp = SMI_RSP_byDay{d};
        if isempty(rsp), continue; end

        % build label from actual days that contributed
        if isfield(RSPData(1), 'dayofexperience')
            dStr = strjoin(arrayfun(@num2str, daysPresent_RSP{d}, 'UniformOutput', false), '+');
        else
            dStr = '200';
        end
        fprintf('Day %s:  n=%d,  median SMI = %.3f\n', dStr, length(rsp), median(rsp));
    end

    if ~isempty(pooledSMI_VISp)
        fprintf('\n--- VISp SMI (pooled, comparison) ---\n');
        fprintf('All days:  n=%d,  median SMI = %.3f\n', length(pooledSMI_VISp), median(pooledSMI_VISp));
    end

    % Ranksum: RSP day 5/200 vs VISp pooled
    rsp_last = SMI_RSP_byDay{daysOfInterest == 5};
    if ~isempty(rsp_last) && ~isempty(pooledSMI_VISp)
        [pVal, ~] = ranksum(rsp_last, pooledSMI_VISp);
        fprintf('\nRanksum RSP day 5/200 vs VISp pooled: p = %.4e\n', pVal);
    end

    %% --- Figure ---
    figHandle = figure('Name', 'SMI by Day: RSP vs VISp', ...
                       'Color', [1 1 1], 'Position', [100 100 700 700]);
    hold on;

    % Plot VISp as single dashed comparison line
    if ~isempty(pooledSMI_VISp)
        [fV, xV] = ecdf(pooledSMI_VISp);
        plot(xV, fV, '--', 'LineWidth', 1.5, 'Color', vispColor, ...
             'DisplayName', sprintf('VISp somas pooled (n=%d)', length(pooledSMI_VISp)));
    end

    % Plot RSP per day
    for d = 1:length(daysOfInterest)
        rsp = SMI_RSP_byDay{d};
        if isempty(rsp), continue; end

        [f, x] = ecdf(rsp);

        % Build legend label from actual days present
        if ~isempty(daysPresent_RSP{d}) && isfield(RSPData(1), 'dayofexperience')
            dStr = strjoin(arrayfun(@num2str, daysPresent_RSP{d}, 'UniformOutput', false), '+');
        else
            dStr = '200';
        end
        if daysOfInterest(d) == 5
            dayLabel = sprintf('RSP Day %s (n=%d)', dStr, length(rsp));
        else
            dayLabel = sprintf('RSP Day %d (n=%d)', daysOfInterest(d), length(rsp));
        end

        plot(x, f, 'LineWidth', 2, 'Color', dayColors_RSP(d,:), 'DisplayName', dayLabel);
    end

    % Reference lines
    xline(0,   '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');

    xlabel('Spatial modulation index');
    ylabel('Cumulative probability');
    title('RSP Boutons — SMI across days (VISp pooled for comparison)');
    xlim([-1.1, 1.1]);
    ylim([0, 1.02]);
    yticks([0, 0.5, 1]);
    legend('Location', 'best', 'Interpreter', 'none', 'Box', 'off');
    set(gca, 'Box', 'off');
    defaultAxesProperties(gca);
    offsetAxes(gca);
    axis square;

    %% --- Save ---
    outputDir = 'Z:\ibn-vision\USERS\Sonali\Figures\ThesisFigs\ResultsChapter2-RSP-PostExp\Section2_Fig3.3\smi_RSPLearning_200_VISp';
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    baseFileName = 'rsp_vs_visp_smi_byDay';
    fullSavePath = fullfile(outputDir, baseFileName);
    saveFigureFormats(figHandle, fullSavePath);
end