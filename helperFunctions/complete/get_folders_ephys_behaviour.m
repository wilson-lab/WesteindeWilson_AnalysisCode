function folders = get_folders_ephys_behaviour(rootDir, overwrite)

%% get folder list

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
          contents_ephys = fullfile(fullFolderName,'trialData.mat'); % contains ephys data
          contents_CL = fullfile(fullFolderName,'behaviorData.mat'); % contains behaviour data
          present = exist(contents_ephys,'file')>0 && exist(contents_CL,'file')>0;
          processed_behaviourdata = fullfile(fullFolderName,'pro_behaviourData.mat');
          processed_trialdata = fullfile(fullFolderName,'pro_trialData.mat');
          processed = exist(processed_behaviourdata, 'file') > 0 && exist(processed_trialdata, 'file') > 0 ;  % already processed?
          if ~overwrite 
              present = ~processed && present;
          end
      end
      correct(end+1) = present;
    end 
    folders = filelist(logical(correct));
end