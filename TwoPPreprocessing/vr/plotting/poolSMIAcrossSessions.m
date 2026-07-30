function pooledData = poolSMIAcrossSessions(allData)
    % poolSMIResults: Aggregates SMI values from all sessions.
    % Groups results by TargetArea and DayOfExperience.
    pooledData = struct();
    for s = 1:length(allData)
        sess = allData(s);
        % SMI exists for the Baseline condition
        if ~isfield(sess, 'SMI')
            continue;
        end
        area = sess.TargetArea;
        day = sess.Day;

        % Day 200 is an alias for Day 5 — merge it in rather than letting
        % it form its own separate 'Day200' bucket.
        if day == 200
            day = 5;
        end

        smiValues = sess.SMI.SMI;
        % Remove NaNs (neurons that didn't meet activity/peak criteria)
        validSMI = smiValues(~isnan(smiValues));

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