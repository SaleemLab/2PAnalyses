function pooledData = poolSMIAcrossSessions(allData)
    % poolSMIResults: Aggregates SMI values from all sessions.
    % Groups results by TargetArea and DayOfExperience.
    
    pooledData = struct();
    
    for s = 1:length(allData)
        sess = allData(s);
        
        % SMI exists for the Baseline condition
        if ~isfield(sess.ConditionData, 'Baseline') || ~isfield(sess.ConditionData.Baseline, 'SMI')
            continue;
        end
        
        area = sess.TargetArea;
        day = sess.Day;
        smiValues = sess.ConditionData.Baseline.SMI;
        
        % Remove NaNs (neurons that didn't meet activity/peak criteria)
        validSMI = smiValues(~isnan(smiValues));
        
        % 
        if ~isfield(pooledData, area)
            pooledData.(area).AllSMI = [];
            pooledData.(area).Days = struct();
        end
        
        % Pool into global area list
        pooledData.(area).AllSMI = [pooledData.(area).AllSMI; validSMI];
        
        % Pool into specific day list for that area
        dayLabel = sprintf('Day%d', day);
        if ~isfield(pooledData.(area).Days, dayLabel)
            pooledData.(area).Days.(dayLabel) = [];
        end
        pooledData.(area).Days.(dayLabel) = [pooledData.(area).Days.(dayLabel); validSMI];
    end
end