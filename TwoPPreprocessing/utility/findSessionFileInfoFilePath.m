function sessionFileInfoFilePath = findSessionFileInfoFilePath(mouseID,session)

rootDir = ['Z:' filesep fullfile('ibn-vision','DATA','SUBJECTS',mouseID)];
sessionFileInfoFileName = [mouseID '_' session '_sessionFileInfo.mat'];
sessionFileInfoFilePath = fullfile(rootDir, 'Analysis', session, sessionFileInfoFileName);

end