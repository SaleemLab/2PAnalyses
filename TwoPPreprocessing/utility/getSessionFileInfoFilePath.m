function sessionFileInfoFilePath = getSessionFileInfoFilePath(mouseID, session)

rootDir = ['Z:' filesep fullfile('ibn-vision','DATA','SUBJECTS',mouseID)];
analysisFolder = fullfile(rootDir, 'Analysis/', session);

if ~exist(analysisFolder, 'dir')
    mkdir(analysisFolder);
end

sessionFileInfoFileName = [mouseID '_' session '_sessionFileInfo.mat'];
sessionFileInfoFilePath = fullfile(analysisFolder, sessionFileInfoFileName);
end 