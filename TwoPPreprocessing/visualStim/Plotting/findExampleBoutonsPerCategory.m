function exampleTable = findExampleBoutonsPerCategory(allDotUnits, nPerCategory, r2Thresh, enforceSessionDiversity)
% exampleTable = findExampleBoutonsPerCategory(allDotUnits, nPerCategory, r2Thresh, enforceSessionDiversity)
%
% Scans across ALL boutons (all sessions/mice) in allDotUnits and, for
% each of the 4 Gaussian tuning categories (1=low-pass, 2=high-pass,
% 3=band-pass, 4=trough), ranks candidates by R^2 and returns the top N
% -- optionally enforcing that examples come from different
% sessions/mice, so your example figure isn't accidentally showing 4
% boutons from the same animal.
%
% ASSUMES allDotUnits has these fields (adjust names below if yours
% differ):
%   .gaussChar    - category code (1-4) from fitGaussianTemplates_tuning
%   .bestR2       - R^2 for that fit (descriptive or cross-validated,
%                   whichever you use for classification/reporting)
%   .mouseID      - string
%   .sessionName  - string
%   .roiIdx       - numeric ROI/bouton index within session
%
% INPUTS
%   allDotUnits             : your main struct array, one row per bouton
%   nPerCategory            : how many example candidates to return per
%                             category (default 5, so you have a
%                             shortlist to eyeball rather than just one)
%   r2Thresh                : minimum R^2 to be considered at all
%                             (default 0.3 -- a good example should fit
%                             well above your inclusion threshold, not
%                             just barely pass it)
%   enforceSessionDiversity : if true, skips candidates from a
%                             session/mouse already represented in that
%                             category's shortlist (default true)
%
% OUTPUT
%   exampleTable : table with columns
%     Category, CategoryName, BoutonIdx, MouseID, SessionName, RoiIdx, R2
%   sorted by Category then descending R2 -- use BoutonIdx to index
%   directly into allDotUnits, or MouseID/SessionName/RoiIdx to call
%   plotDotFieldsExampleBouton(...) for a full trial-level look.

if nargin < 2 || isempty(nPerCategory), nPerCategory = 5; end
if nargin < 3 || isempty(r2Thresh), r2Thresh = 0.3; end
if nargin < 4 || isempty(enforceSessionDiversity), enforceSessionDiversity = true; end

categoryNames = {'Low-Pass', 'High-Pass', 'Band-Pass', 'Trough'};

allCategory   = cat(1, allDotUnits.gaussChar);
allR2         = cat(1, allDotUnits.bestR2);
allMouseID    = {allDotUnits.mouseID}';
allSessionName = {allDotUnits.sessionName}';
allRoiIdx     = cat(1, allDotUnits.roiIdx);
allBoutonIdx  = (1:numel(allDotUnits))';

Category = []; CategoryName = {}; BoutonIdx = []; MouseID = {}; SessionName = {}; RoiIdx = []; R2 = [];

for cat_i = 1:4
    candIdx = find(allCategory == cat_i & allR2 >= r2Thresh);
    if isempty(candIdx)
        warning('No boutons found for category %d (%s) at R^2 >= %.2f -- lowering threshold or check gaussChar coding.', ...
            cat_i, categoryNames{cat_i}, r2Thresh);
        continue
    end

    % sort candidates for this category by R^2 descending
    [~, sortOrder] = sort(allR2(candIdx), 'descend');
    candIdx = candIdx(sortOrder);

    selectedForCat = [];
    usedSessions = {};

    for k = 1:numel(candIdx)
        thisIdx = candIdx(k);
        thisSessionKey = [allMouseID{thisIdx}, '_', allSessionName{thisIdx}];

        if enforceSessionDiversity && ismember(thisSessionKey, usedSessions)
            continue % skip -- already have an example from this session for this category
        end

        selectedForCat(end+1) = thisIdx; %#ok<AGROW>
        usedSessions{end+1} = thisSessionKey; %#ok<AGROW>

        if numel(selectedForCat) >= nPerCategory
            break
        end
    end

    for k = 1:numel(selectedForCat)
        idx = selectedForCat(k);
        Category(end+1,1)     = cat_i; %#ok<AGROW>
        CategoryName{end+1,1} = categoryNames{cat_i}; %#ok<AGROW>
        BoutonIdx(end+1,1)    = allBoutonIdx(idx); %#ok<AGROW>
        MouseID{end+1,1}      = allMouseID{idx}; %#ok<AGROW>
        SessionName{end+1,1}  = allSessionName{idx}; %#ok<AGROW>
        RoiIdx(end+1,1)       = allRoiIdx(idx); %#ok<AGROW>
        R2(end+1,1)           = allR2(idx); %#ok<AGROW>
    end
end

exampleTable = table(Category, CategoryName, BoutonIdx, MouseID, SessionName, RoiIdx, R2);

fprintf('\n--- Example bouton candidates (R^2 >= %.2f, session diversity = %d) ---\n', r2Thresh, enforceSessionDiversity);
for cat_i = 1:4
    nFound = sum(exampleTable.Category == cat_i);
    fprintf('%-12s: %d candidates found\n', categoryNames{cat_i}, nFound);
end

end
