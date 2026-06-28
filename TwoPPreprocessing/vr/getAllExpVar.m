function ExpVar = getAllExpVar(activity, nfold)
% activity : trial x position [changed from position x trials]
% Aman Saleem - May 2026
% @Aman - should the training set accumulate trials across folds? 

if nargin<2
    nfold = 5;
end

% Since activity is trial x position, number of trials is the 1st dimension
numTrials = size(activity, 1); 
ExpVar = ones(1,nfold)*nan;
CV = newCrossVal(numTrials,nfold);

for m = 1:CV.nfold
    trainMeanCurve  = calMeanCurve(activity, CV.trainTrials{m});
    testMeanCurve   = calMeanCurve(activity, CV.testTrials{m});
    ExpVar(m) = newCalExpVar(trainMeanCurve, testMeanCurve);
end

    %% helper functions
    function ExpVar = newCalExpVar(trainMeanCurve, testMeanCurve)
        %Calculate explained variance from trainmean and test mean curves -
        %FYI: this is Explained variance by position rather than time -
        %values should be higher
        % totSS — how much the test curve varies around its own mean (baseline variance)
        % resSS — how much the test curve deviates from the train mean
        % curve 
        
        trainMean = mean(trainMeanCurve, 'omitnan'); % included omitnans
        
        % Calculate sum of squares while omitting NaNs frame-by-frame
        % does the test curve match the train curve shape + magnitude
        resSS = sum((testMeanCurve - trainMeanCurve).^2, 'omitnan');
        % does the test curve deviate from a flat baseline of overall train activity
        totSS = sum((testMeanCurve - trainMean).^2, 'omitnan');
        
        ExpVar = 1 - resSS/totSS;
        % if ExpVar<0
        %     ExpVar=0;
        % end 
    end

    function meanCurve = calMeanCurve(activity, ntrialList)
        % Simple function to calculate the mean curve from a list of trials
        % activity : trial x position
        % Select target rows (trials) and average down column dimension 1
        meanCurve = mean(activity(ntrialList, :), 1, 'omitnan'); 
    end

    %This version accumulates trials across folds; the first fold also
    %looks to be 0 it
    % function CV = newCrossVal(numTrials,nfold)
    %     % make a cross validation set based on the number of trials and with nfold
    % 
    %     CV.numTrials = numTrials;
    %     CV.nfold = nfold;
    % 
    % 
    %     perms = randperm(numTrials);
    % 
    %     kidx = floor(length(perms)./nfold).*(0:nfold);
    %     kidx(end) = length(perms);
    % 
    %     for n = 1:nfold
    %         % accumulate trials across folds? why? @Aman 
    %         CV.trainTrials{n} = perms([1:kidx(n) (kidx(n+1)+1):kidx(end)]);
    %         CV.testTrials{n} = perms((kidx(n)+1):kidx(n+1));
    %     end
    % end


    function CV = newCrossVal(numTrials,nfold)
        % make a cross validation set based on the number of trials and with nfold
        CV.numTrials = numTrials;
        CV.nfold = nfold;
        perms = randperm(numTrials);
        % kidx = floor(length(perms)./nfold).*(0:nfold);
        % kidx(end) = length(perms); % this seemed to generate unequal number
        % of trials for the last fold
        kidx = round(linspace(0, numTrials, nfold + 1));
        for n = 1:nfold
            currentTest = perms(kidx(n)+1:kidx(n+1));
            CV.testTrials{n} = currentTest;

            currentTrain = perms;
            % ismember to find and drop the test trial
            currentTrain(ismember(currentTrain, currentTest)) = [];
            CV.trainTrials{n} = currentTrain(:)';

            % archived
            % For fold 1, kidx(1) equals 0, making both bracket ranges
            % evaluate as empty (1:0 and kidx(end)+1:kidx(end)), resulting
            % in an empty trainTrials array and subsequent NaNs.
            % CV.trainTrials{n} = perms([kidx(1)+1:kidx(n) kidx(nfold+1)+1:kidx(end)]);
            % CV.testTrials{n} = perms(kidx(n)+1:kidx(n+1));
            % CV.trainTrials{n} = perms([1:kidx(n) (kidx(n+1)+1):kidx(end)]);
            % CV.testTrials{n} = perms((kidx(n)+1):kidx(n+1));
        end
    end

end
