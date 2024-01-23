function folders = get_folders_ephys(rootDir)
% get experiment folder list from directory

    filelist = dir(fullfile(rootDir, '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist([filelist.isdir]);  % only folders from list
    filelistNum = length(filelist);
    correct = [];
    for f = 1:filelistNum
      baseFolderName = filelist(f).name;
      if strcmp(baseFolderName,'..') || ~strcmp(baseFolderName,'.')
          present = false;
      else
          parentFolderName = filelist(f).folder;
          fullFolderName = fullfile(parentFolderName, baseFolderName);
          filelist(f).folder = fullFolderName;
          contents = fullfile(fullFolderName,'trialData.mat'); 
          present = exist(contents,'file')>0;
          listingTiff = dir(fullfile(fullFolderName,'pro_trialData.mat'));
          present = numel(listingTiff) > 0 && present; 
      end
      correct(end+1) = present;
    end 
    folders = filelist(logical(correct));
end