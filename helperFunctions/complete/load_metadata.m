function [expMetadata, trialMetadata, patternMetadata, fictracMetadata] =  load_metadata(expList, parentDir)
% ==================================================================================================   
%  Loads the experiment and trial metadata files for a set of expIDs and compiles them into two 
%  tables.
%
%  INPUTS: 
%       expList                 = cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files
%
%
%  OUTPUTS:
%       expMetadata             = table with these columns of experiment-specific info:
%                                       expID
%                                       expName
%                                       daqSampRate
%                                       panelsDisplayRate
%                                       volumeRate
%                                       nPlanes
%                                       nTrials
%
%       trialMetadata           = table with these columns of trial-specific info:
%                                       expID
%                                       trialNum
%                                       trialDuration
%                                       nVolumes
%                                       nDaqSamples
%                                       nPanelsFrames
%                                       usingOptoStim
%                                       optoStimTiming
%                                       usingPanels
%                                       using2P
%                                       originalTrialCount
%                                       pmtShutoffVols
%
% ==================================================================================================


% Parse/process input arguments
if nargin < 2
   parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments'; 
end
if isa(expList, 'table')
   expList = expList.expID; 
end

% Loop through each expID
expMetadata = [];
trialMd = [];
disp('------------------------------------------');
disp('Loading metadata files...')
for iExp = 1:numel(expList)
   currExpID = expList{iExp};
   expMdFile = fullfile(parentDir,'csv', 'expMd.csv');       
   
   % Load pattern metaData
   try
       mdFiles = dir(fullfile(parentDir,['*' currExpID,'_metadata_*.mat']));
       mdFile = mdFiles.name; 
       load(fullfile(parentDir,mdFile),'mD');
       if isfield(mD.trialSettings,'pattern')
           
%            if isempty(mD.trialSettings.pattern)
%                 disp("---Warning Pattern metaData empty, this is a temp fix, line 65 load metadata---")
%                 filepath = 'Z:\2photon_data\data2process\needs_python\20220413\20220413-1_PFL2_Fly1_LAL'; 
%                 mdFile2 = dir(fullfile(filepath, ['*metadata*trial_', pad(num2str(1), 3, 'left', '0'), '.mat']));
%                 old_mD = load(fullfile(filepath, mdFile2.name), 'mD');
%                 mD.trialSettings.pattern = old_mD.mD.trialSettings.pattern; 
%            end
            
           patternMetadata = mD.trialSettings.pattern;
           patternMetadata.initialPos = mD.experiment.panels.initialPosition;   
       end
       if isfield(mD,'experiment')
            if isfield(mD,'panels')
                for fn = fieldnames(mD.experiment.panels)'
                    fieldData = mD.experiment.panels.(fn{1});
                    if isempty(fieldData)
                        fieldData = NaN;
                    end
                    patternMetadata.(fn{1}) = fieldData;
                 end 
            end
       end
   catch
       disp("couldn't load pattern metaData"); 
   end
   
   % load fictrac metaData
   try
    fictracMetadata = mD.fictrac; 
   catch
       disp("fictracmetdaData no found")
   end
   if exist(expMdFile, 'file')
       disp(['Adding ', currExpID, '...'])
       
       % Load expMetadata file if it exists
       currExpMd = readtable(expMdFile, 'delimiter', ',');
       expMetadata = [expMetadata; currExpMd];
       
        % Also load trial metadata file 
        trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.mat']);
        if exist(trialMdFile,'file')
            load(trialMdFile, 'trialMetadata');
            currTrialMd = trialMetadata;
            if ~ismember(fieldnames(currTrialMd), 'optoStimTiming')
            currTrialMd.optoStimTiming = repmat({[]}, size(currTrialMd, 1), 1);
            end
            trialMd = [trialMd; currTrialMd];
        else
            disp('Skipping trialMetadata...file not found');
        end
   else
       disp(['Skipping ', currExpID, '...file not found']);
   end
   
end
trialMetadata = trialMd;
disp('All metadata loaded')
end