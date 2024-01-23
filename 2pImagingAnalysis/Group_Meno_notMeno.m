function [MenoData, nMenoData] = Group_Meno_notMeno(rootDir, PFL3, PFL2, window, minVel, highThres, lowThres)
% Iterates through trials and segments data into likely menotaxis vs not
% menotaxing groups. Recalculates estimated goal direction and rho for each
% completed segment. 

% INPUTS

% PFL3: include PFL3 experiments, boolean
% PFL2: include PFL2 experiments, boolean
% rootDir: directory containing experiment folders, string
% window: time window to calculate rho over for each datapoint (s), float
% minVel: velocity threshold total velocity must be over to count as moving, float
% highThres: rho value to be above to start goal directed segment, float
% low Thres: rho value when below that ends goal directed segment, float

% OUTPUTS

% MenoData: table summarizing goal directed segments
% nMenoData: table summarizing non-goal directed segments
    
    
    % setup output tables
    varNames = ["Folder","numTrial","Goal","rho","timeMov","Indices"];
    varTypes = ["string","double","double","double","double","cell"];
    MenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
    nMenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
    
    % gather experiment folders & iterate through them
    folders = get_folders(rootDir, PFL2, PFL3);

    mCount = 1; 
    nCount = 1; 

    for f = 1:size(folders,1)
        folder = folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        expID = get_expID(folder);
        expList = {expID};

        % Load metadata 
        [~, trialMd] = load_metadata(expList, folder);

        % Load imaging data
        roiData = load_roi_data(expList, folder);

        processedData_dir = fullfile(folder,'processed_data');
        try
            numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum));
        catch
            numTrials = 1; 
        end

        for nTrial = 1:numTrials

            data_filelist = dir(processedData_dir);
            % load correct ftT file
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
            end

            Meno_sum = {};
            not_Meno_sum = {};

            % categorize data into 'menotaxing' vs 'non-menotaxing' chunks
            [~, ~, Meno_chunks, not_Meno_chunks] = SegmentMenovsNotMeno(ftT, window, minVel,highThres,lowThres,folder,nTrial,trialMd.volumeRate,1);

            if ~isempty(Meno_chunks) || ~isempty(not_Meno_chunks)
                Meno_chunks = Meno_chunks';
                not_Meno_chunks = not_Meno_chunks';

                % catch any chunks with problems (unable to calculate an average HD)
                count = 1; 
                for chunk = 1:length(Meno_chunks)
                    [~, theta] = CalculateAverageHeading(ftT,minVel, Meno_chunks{chunk,1});
                    if ~isnan(theta)
                        Meno_sum{count,1} = Meno_chunks{chunk,1}; 
                        count = count + 1; 
                    end
                end
                count = 1; 
                for chunk = 1:length(not_Meno_chunks)
                    [~, theta] = CalculateAverageHeading(ftT,minVel, not_Meno_chunks{chunk,1});
                    if ~isnan(theta)
                        not_Meno_sum{count,1} = not_Meno_chunks{chunk,1}; 
                        count = count + 1; 
                    end
                end

                % recalculate final segment goals & collect into summary tables
                while ~isempty(Meno_sum)
                    MenoData.Indices{mCount} = Meno_sum{1,1};
                    [rho, theta] = CalculateAverageHeading(ftT,minVel, MenoData.Indices{mCount});
                    speed = sqrt(ftT.velSide{1}(MenoData.Indices{mCount}).^2 + ftT.velFor{1}(MenoData.Indices{mCount}).^2);
                    timeMov = sum(speed > 1.5)/60;
                    MenoData.Goal(mCount) = theta; 
                    MenoData.rho(mCount) = rho;
                    MenoData.timeMov(mCount) = timeMov; 
                    MenoData.Folder(mCount) = folder; 
                    MenoData.numTrial(mCount) = nTrial;
                    Meno_sum(1,:) = [];
                    mCount = mCount + 1;
                end


                while ~isempty(not_Meno_sum)
                    nMenoData.Indices{nCount} = not_Meno_sum{1,1};
                    [rho, theta] = CalculateAverageHeading(ftT,minVel, nMenoData.Indices{nCount});
                    speed = sqrt(ftT.velSide{1}(nMenoData.Indices{nCount}).^2 + ftT.velFor{1}(nMenoData.Indices{nCount}).^2);
                    timeMov = sum(speed > 1.5)/60;
                    nMenoData.Goal(nCount) = theta; 
                    nMenoData.rho(nCount) = rho;
                    nMenoData.timeMov(nCount) = timeMov;
                    nMenoData.Folder(nCount) = folder;
                    nMenoData.numTrial(nCount) = nTrial;
                    not_Meno_sum(1,:) = []; 
                    nCount = nCount + 1;
                end
            end
        end
    end
end
