function [allFtData, allFtData_DAQ, allFtData_dat] = load_ft_data(expList, parentDir, DAQ, dat)
% ==================================================================================================   
%  Loads FicTrac data files for a set of expIDs and compiles them into a single table. 
%
%  INPUTS: 
%       expList                 = cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files

%       DAQ                     = load in DAQ fictrac data
%       dat                     = load in dat fictrac data
%
%  OUTPUTS:
%       allFtData               = table with these columns of FicTrac data for each trial:
%                                   expID
%                                   trialNum
%                                   intX
%                                   intY
%                                   intHD
%                                   moveSpeed
%                                   intFwMove
%                                   intSideMove
%                                   yawSpeed
%                                   fwSpeed
%                                   sideSpeed
%                                   frameTimes   
%                                   badVidFrames    
%                                   meanFlow                                       
%
% ==================================================================================================

% Parse/process input arguments
if nargin < 2
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
elseif nargin == 1
    DAQ = 0 ;
    dat = 0 ;
end

if strcmp(parentDir(end),'.')
    parentDir = parentDir(1:end-2); 
end 
    
    
if isa(expList, 'table')
    expList = expList.expID;
end

% Load FicTrac data file for each expID if it exists
allFtData = [];
allFtData_dat = [];
allFtData_DAQ = [];

disp('------------------------------------------');
disp('Loading FicTrac data files...')
for iExp = 1:numel(expList)
   currExpID = expList{iExp};
   
   if dat && DAQ
       
       ftDataFile_dat = fullfile(parentDir, [currExpID, '_ficTracData_dat.mat']);
       if exist(ftDataFile_dat, 'file')
           disp(['Loading ', currExpID, '...'])
           load(ftDataFile_dat, 'ftData_dat');
           allFtData_dat = [allFtData_dat; ftData_dat];
       else
           disp(['Skipping ', currExpID, '...file not found']);
       end
       
       ftDataFile_DAQ = fullfile(parentDir, [currExpID, '_ficTracData_DAQ.mat']);
       if exist(ftDataFile_DAQ, 'file')
           disp(['Loading ', currExpID, '...'])
           load(ftDataFile_DAQ, 'ftData_DAQ');
           allFtData_DAQ = [allFtData_DAQ; ftData_DAQ];
       else
           disp(['Skipping ', currExpID, '...file not found']);
       end
   elseif dat
       ftDataFile_dat = fullfile(parentDir, [currExpID, '_ficTracData_dat.mat']);
       if exist(ftDataFile_dat, 'file')
           disp(['Loading ', currExpID, '...'])
           load(ftDataFile_dat, 'ftData_dat');
           allFtData_dat = [allFtData_dat; ftData_dat];
       else
           disp(['Skipping ', currExpID, '...file not found']);
       end
       
   elseif DAQ
       ftDataFile_DAQ = fullfile(parentDir, [currExpID, '_ficTracData_DAQ.mat']);
       if exist(ftDataFile_DAQ, 'file')
           disp(['Loading ', currExpID, '...'])
           load(ftDataFile_DAQ, 'ftData_DAQ');
           allFtData_DAQ = [allFtData_DAQ; ftData_DAQ];
       else
           disp(['Skipping ', currExpID, '...file not found']);
       end
       
   else
      ftDataFile = fullfile(parentDir, [currExpID, '_ficTracData.mat']);
      if exist(ftDataFile, 'file')
           disp(['Loading ', currExpID, '...'])
           load(ftDataFile, 'ftData');
           allFtData = [allFtData; ftData];
       else
           disp(['Skipping ', currExpID, '...file not found']);
       end
   end
   
end
disp('FicTac data loaded')
end