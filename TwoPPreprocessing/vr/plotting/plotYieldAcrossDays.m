function plotYieldAcrossDays(allData)

    saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\DistibutionsAllCritera';
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    
    % Thresholds
    T_CORR = 0.4;
    days = 1:5;
    

    miceIDs = unique({allData.MouseID});
    numMice = length(miceIDs);
    
    % NaN so missing data doesn't plot as zero
    mousePct = nan(numMice, 5);
    pooledTotal = zeros(5, 1);
    pooledPassing = zeros(5, 1);

    for d = 1:5

        if isnumeric(allData(1).Day)
            dayIdx = [allData.Day] == d;
        else
            dayIdx = str2double({allData.DayOfExperience}) == d;
        end
        dayData = allData(dayIdx);
        
        if isempty(dayData), continue; end
        

        for i = 1:numMice
            thisMouse = miceIDs{i};
            mSessIdx = strcmp({dayData.MouseID}, thisMouse);
            mSessions = dayData(mSessIdx);
            
            if isempty(mSessions), continue; end % Skip if mouse not in this day
            
            mTotal = 0; mPass = 0;
            for s = 1:length(mSessions)
                sess = mSessions(s);
                try
                    % Logic: (Significant Peak OR Range) AND (Halves Rho >= 0.4)
                    isSig = (sess.isSignificantByPeakShuffling == 1) | (sess.isSignificantByRange == 1);
                    isStable = isSig & (sess.lapCorr_HalvesRho >= T_CORR);
                    
                    mTotal = mTotal + sess.NumCells;
                    mPass  = mPass  + sum(isStable);
                catch
       
                end
            end
            if mTotal > 0
                mousePct(i, d) = (mPass / mTotal) * 100;
            end
        end
        
        %pooled totals for the day 
        for s = 1:length(dayData)
            try
                sess = dayData(s);
                isSig = (sess.isSignificantByPeakShuffling == 1) | (sess.isSignificantByRange == 1);
                isStable = isSig & (sess.lapCorr_HalvesRho >= T_CORR);
                
                pooledTotal(d)   = pooledTotal(d) + sess.NumCells;
                pooledPassing(d) = pooledPassing(d) + sum(isStable);
            catch
            end
        end
    end
    

    pooledPct = nan(1, 5);
    validDays = pooledTotal > 0; %sanity check
    pooledPct(validDays) = (pooledPassing(validDays) ./ pooledTotal(validDays)) * 100;


    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.6]);
    hold on;
    
    %
    hMouse = plot(days, mousePct', '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
    
    %  pooled population as a thick black line
    hPooled = plot(days, pooledPct, '-ko', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    

    xlabel('Day of Experience', 'FontSize', 12);
    ylabel('% of ROIs Passing Criteria', 'FontSize', 12);
    title({'Population yield across experience', 'Criteria: SigShuffle (Peak|Range) & Halves Corr rho \geq 0.4'}, 'FontSize', 14);
    
    xticks(days);
    xlim([0.8 5.2]);
    ylim([0 100]); % Set Y-axis 0-100 for percentage
    grid on;
    

    legend([hPooled, hMouse(1)], {'Pooled Population', 'Individual Mice'}, 'Location', 'best');
    

    saveName = fullfile(saveDir, 'Yield_Trend_Sig_Halves.png');
    exportgraphics(fig, saveName, 'Resolution', 300);
    fprintf('Figure saved to: %s\n', saveName);
end