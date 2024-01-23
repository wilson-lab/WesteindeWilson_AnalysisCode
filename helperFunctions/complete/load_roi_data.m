function allRoiData = load_roi_data(expList, parentDir, channel)
% ==================================================================================================   
%  Loads ROI data files for a set of expIDs and compiles them into a single table. 
%  
%  INPUTS: 
%       expList                 = cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files
%
%
%  OUTPUTS:
%       allRoiData              = table with these columns of ROI data for each trial:
%                                   expID
%                                   trialNum
%                                   roiName
%                                   subROIs
%                                   rawFl
%                                   dffData
%                                   expDffData
%
% ==================================================================================================

% Parse/process input arguments

%% Arguments
arguments
    expList
    parentDir char
    channel char = 'ch1'
end

if isa(expList, 'table')
    expList = expList.expID;
end

% Load ROI data file for each expID if it exists
allRoiData = [];
disp('------------------------------------------');
disp('Loading ROI data files...')
for iExp = 1:numel(expList)
   currExpID = expList{iExp};
   if strcmp('ch1', channel)
       roiDataFile = fullfile(parentDir, [currExpID, '_roiData_' channel '.mat']);
       if ~exist(roiDataFile, 'file')
           roiDataFile = fullfile(parentDir, [currExpID, '_roiData.mat']);
       end
   else
       roiDataFile = fullfile(parentDir, [currExpID, '_roiData_' channel '.mat']);
   end
   if exist(roiDataFile, 'file') 
       disp(['Loading ', currExpID, '...'])
       load(roiDataFile, 'roiData');
       allRoiData = [allRoiData; roiData];
   else
       disp(['Skipping ', currExpID, '...file not found']);
   end
end
disp('ROI data loaded');
end