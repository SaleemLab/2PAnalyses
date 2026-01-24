function plotDay4SurvivorCurves(allData)

    saveDir = '\\rdp.arc.ucl.ac.uk\ritd-ag-project-rd01ie-asale69\ibn-vision\USERS\Sonali\Figures\DistibutionsAllCritera';
    if ~exist(saveDir, 'dir'), mkdir(saveDir); end


    T_CORR = 0.4;
    T_VAR_RANGE = 1;
    T_VAR_VAR = 20;


    if isnumeric(allData(1).Day)
        dayIdx = [allData.Day] == 4;
    else
        dayIdx = str2double({allData.Day}) == 4;
    end
    day4Data = allData(dayIdx);

    if isempty(day4Data)
        error('No data found for DayOfExperience == 4. Please verify the data structure.');
    end

    mice = unique({day4Data.MouseID});
    numMice = length(mice);


    filterLabels = {'Total ROIs', 'SigShuffle Peak|Range', '+ Halves Corr', '+ Tuning Range Ratio', '+ Tuning Var Ratio'};
    

    mouseCounts = zeros(numMice, 5);

    for i = 1:numMice
        mouseID = mice{i};
        mouseIdx = strcmp({day4Data.MouseID}, mouseID);
        mouseSessions = day4Data(mouseIdx);
        
        m_total = 0; m_sig = 0; m_stable = 0; m_range = 0; m_tuning = 0;

        for s = 1:length(mouseSessions)
            sess = mouseSessions(s);
            try
                % total rois
                totalCount = sess.NumCells;
                
               
                isSig = (sess.isSignificantByPeakShuffling == 1) | (sess.isSignificantByRange == 1);
                
           
                isStable = isSig & (sess.lapCorr_HalvesRho >= T_CORR);
                
          
                isRangeRatio = isStable & (sess.ratioVarToTuningRange <= T_VAR_RANGE);
                
              
                isTuningVar = isRangeRatio & (sess.ratioVarToTuningVar <= T_VAR_VAR);

            
                m_total  = m_total  + totalCount;
                m_sig    = m_sig    + sum(isSig);
                m_stable = m_stable + sum(isStable);
                m_range  = m_range  + sum(isRangeRatio);
                m_tuning = m_tuning + sum(isTuningVar);
            catch ME
                warning('Error processing Mouse %s: %s', mouseID, ME.message);
            end
        end
        mouseCounts(i, :) = [m_total, m_sig, m_stable, m_range, m_tuning];
    end

  
    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.6]);
    hold on;
    
    colors = lines(numMice);
    
    
    
    for i = 1:numMice
     
        plot(1:5, mouseCounts(i, :), '-o', 'Color', colors(i,:), ...
            'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
            'DisplayName', sprintf('Mouse %s', mice{i}));
 
        text(5.05, mouseCounts(i, 5), sprintf('n=%d', mouseCounts(i, 5)), ...
            'Color', colors(i,:), 'FontWeight', 'bold', 'FontSize', 9);
    end


    xticks(1:5);
    xticklabels(filterLabels);
    xtickangle(25); % angle
    ylabel('Number of ROIs');
    title('ROI Sequential Filtering: Day 4 (Per Animal)', 'FontSize', 14);
    
    grid on;
    legend('Location', 'northeastoutside');
    

    ylim([0 max(mouseCounts(:)) * 1.1]);
    xlim([0.7 5.6]);


    saveName = fullfile(saveDir, 'Day4_ROI_Yield_Curves.png');
    exportgraphics(fig, saveName, 'Resolution', 300);
    fprintf('Analysis complete. Day 4 yield plot saved to: %s\n', saveName);
end