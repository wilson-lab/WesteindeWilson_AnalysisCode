function [MenoData, nMenoData] = Group_Meno_notMeno_ephys(rootDir,saveData)
% split electrophysiology data into menotaxing vs not menotaxing chunks
% based on path straightness & movement. Collect data across all trials 

varNames = ["Folder","numTrial","Goal","rho","timeMov","Indices"];
varTypes = ["string","double","double","double","double","cell"];
MenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
nMenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);


folders = get_folders_ephys(rootDir);

mCount = 1; 
nCount = 1; 

window = 30; % window length (s)
highThres = 0.70; % threshold when crossed with a + slope indicates start of meno bought
lowThres = 0.70; % threshold when crossed with a - slope indicates end of meno bought

minVel = 1.5; 

% iterate through each trial folder
for f = 1:size(folders,1)
    
    folder = folders(f).folder; 
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end

    % load in data
    processedDir = fullfile(folder,'processedData');
    try
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        load(fullfile(processedDir,'pro_trialData.mat'))
    catch
        load(fullfile(folder,'pro_behaviourData.mat'))
        pro_behaviourData = processed_behaviourData;
        load(fullfile(folder,'pro_trialData.mat'))
        pro_trialData = processed_trialData;
    end
    
    numTrials = length(pro_behaviourData);
        
    for nTrial = 1:numTrials
    mCountFly = 1;
    nCountFly = 1;
    MenoDataFly = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
    nMenoDataFly = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
  
    Meno_chunks = [];
    not_Meno_chunks = []; 
    bData = pro_behaviourData{nTrial};
    tData = pro_trialData{nTrial};
    
    % perform index by index categorizartion
    [~, ~, Meno_chunks, not_Meno_chunks,bData_interp] = SegmentMenovsNotMeno_ephys(bData,window, minVel,highThres,lowThres,folder);
    
    % assign datapoints to meno or not menotaxing chunks based on segmentation & collect in cell
    if ~isempty(Meno_chunks) || ~isempty(not_Meno_chunks)
        Meno_chunks = Meno_chunks';
        not_Meno_chunks = not_Meno_chunks';

        count = 1; 
        for chunk = 1:length(Meno_chunks)
            interpIdx = Meno_chunks{chunk,1};
            interpTime = bData_interp.time(interpIdx);
            idx = find(bData.time >= interpTime(1) & bData.time <= interpTime(end));
            [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, idx);
            if ~isnan(theta)
                Meno_sum{count,1} = idx'; 
                Meno_sum{count,2} = theta; 
                Meno_sum{count,3} = rho;
                count = count + 1; 
            end
        end
        count = 1; 
        for chunk = 1:length(not_Meno_chunks)
            interpIdx = not_Meno_chunks{chunk,1};
            interpTime = bData_interp.time(interpIdx);
            idx = find(bData.time > interpTime(1) & bData.time < interpTime(end));
            [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, idx);
            if ~isnan(theta)
                not_Meno_sum{count,1} = idx'; 
                not_Meno_sum{count,2} = theta; 
                not_Meno_sum{count,3} = rho;
                count = count + 1; 
            end
        end
    

        % collect chunks with sim goals (no longer doing, set accepted goal 
        % diff to 0 but left code here as an option)
        if ~isempty(Meno_chunks)
            while ~isempty(Meno_sum)
                g = Meno_sum{1,2}; 
                other_g = [Meno_sum{2:end,2}];
                diffg = abs(wrapToPi(g - other_g));
                simg = find(diffg < 0) + 1; % allow a +/- x degree diff around (old option)
                MenoData.Indices{mCount} = [Meno_sum{[1,simg],1}];
                [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, MenoData.Indices{mCount});
                speed = sqrt(bData.vel_side(MenoData.Indices{mCount}).^2 + bData.vel_for(MenoData.Indices{mCount}).^2);
                timeMov = sum(speed > minVel)/1000;
                MenoData.Goal(mCount) = theta; 
                MenoData.rho(mCount) = rho;
                MenoData.timeMov(mCount) = timeMov; 
                MenoData.Folder(mCount) = folder; 
                MenoData.numTrial(mCount) = nTrial;

                MenoDataFly.Indices{mCountFly} = [Meno_sum{[1,simg],1}];
                MenoDataFly.Goal(mCountFly) = theta; 
                MenoDataFly.rho(mCountFly) = rho;
                MenoDataFly.timeMov(mCountFly) = timeMov; 
                MenoDataFly.Folder(mCountFly) = folder; 
                MenoDataFly.numTrial(mCountFly) = nTrial;
                Meno_sum([1, simg],:) = [];
                mCount = mCount + 1;
                mCountFly = mCountFly + 1;
            end
        end
        
        if ~isempty(not_Meno_chunks)
            while ~isempty(not_Meno_sum)
                g = not_Meno_sum{1,2}; 
                other_g = [not_Meno_sum{2:end,2}];
                diffg = abs(wrapToPi(g - other_g));
                simg = find(diffg < 0) + 1; % allow a +/- x degree diff around  (old option)
                nMenoData.Indices{nCount} = [not_Meno_sum{[1,simg],1}];
                [rho, theta] = CalculateAverageHeading_ephys(bData,minVel, nMenoData.Indices{nCount});
                speed = sqrt(bData.vel_side(nMenoData.Indices{nCount}).^2 + bData.vel_for(nMenoData.Indices{nCount}).^2);
                timeMov = sum(speed > minVel)/1000;
                nMenoData.Goal(nCount) = theta; 
                nMenoData.rho(nCount) = rho;
                nMenoData.timeMov(nCount) = timeMov;
                nMenoData.Folder(nCount) = folder;
                nMenoData.numTrial(nCount) = nTrial;

                nMenoDataFly.Indices{nCountFly} = [not_Meno_sum{[1,simg],1}];
                nMenoDataFly.Goal(nCountFly) = theta; 
                nMenoDataFly.rho(nCountFly) = rho;
                nMenoDataFly.timeMov(nCountFly) = timeMov;
                nMenoDataFly.Folder(nCountFly) = folder;
                nMenoDataFly.numTrial(nCountFly) = nTrial;
                not_Meno_sum([1, simg],:) = [];
                nCount = nCount + 1;
                nCountFly = nCountFly + 1;
            end
        end
        
        if saveData
            save(fullfile(folder,['MenoDataFly_trial_',num2str(nTrial),'.mat']),'MenoDataFly');
            save(fullfile(folder,['nMenoDataFly_trial_',num2str(nTrial),'.mat']),'nMenoDataFly');
        end
        
    end
    end
    
end

if saveData
    save(fullfile(rootDir,'MenoData.mat'),'MenoData')
    save(fullfile(rootDir,'nMenoData.mat'),'nMenoData')
end
    
end
