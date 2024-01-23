function summarizePFL2_bumpAmp(summaryArray_all) 
    % Averages data across flies & trials to make bump amplitude vs
    % behaviour plots
    
    % INPUTS
    
    % summaryArray_all: table describing segmented PFL2 trials & each segment's estimated goal HD

    edges_vy = [-200:10:200];
    edges_vf = [-4:1:10];
    edges_angle = [-180:10:180];     
    threshold=3;
    fcount = 1;
    
    % only look at FB & PB trials (where the bump amplitude can be estimated)
    summaryArraya = summaryArray_all(contains(summaryArray_all.Folder,'FB'),:) ;
    summaryArrayb = summaryArray_all(contains(summaryArray_all.Folder,'PB'),:);  
    summaryArray = [summaryArraya; summaryArrayb];
    
    all_folders = summaryArray.Folder;
    
    % collect fly ID associated with each trial
    for f = 1:size(all_folders,1)
        [start,finish] = regexp(all_folders(f),'_fly', 'ignorecase');
        dateIdx = regexp(all_folders(f),'\');
        dateIdx = [dateIdx(3)+1:dateIdx(4)-1]; 
        Date = char(all_folders(f));
        Date = Date(dateIdx);
        if isempty(finish)
            [start,finish] = regexp(all_folders(f),'_Fly');
        end
        fly_temp = char(all_folders(f));
        fly = fly_temp(start:finish + 1);
        flyID = strcat(Date,fly);
        flies(f) = string(flyID);
    end
        flies = flies';
        uniqueFlies = unique(flies); 
        flyCount = zeros(size(all_folders));
        for fly = 1:length(uniqueFlies)
            flyCount(flies == uniqueFlies(fly)) = fly;
        end
        
    trial_vy = cell(size(uniqueFlies));
    trial_angle = cell(size(uniqueFlies));
    trial_vf = cell(size(uniqueFlies)); 


    % gets rid of trials w/ no heading change --> indicates problem
    summaryArray = summaryArray(~ismembertol(summaryArray.rho,1,10^-10),:); 
    % only look at trials where fly's vel was above threshold for at least 2 seconds
    summaryArray = summaryArray(summaryArray.timeMov > 2,:); 
 
    % iterate through each trial
    for trial = 1:size(summaryArray,1)
        try
            folder = table2array(summaryArray(trial,1)); 
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

            fly = flyCount(trial);

            % Load bump parameters & behaviour data
            clear bump_params ftT
            processedData_dir = fullfile(folder,'processed_data');
            nTrial = summaryArray.numTrial(trial);
            bump_params = [];
            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name)); 
                end
            end
            
            load(fullfile(processedData_dir,['bump_parameters_Trial00',num2str(nTrial),'.mat']))

            % Gather relevant data & remove idx where the fly isn't moving
            total_mov_mm = abs(ftT.velFor{1}(summaryArray.Indices{trial}) + abs(ftT.velSide{1}(summaryArray.Indices{trial})) + abs(ftT.velYaw{1}(summaryArray.Indices{trial}))*4.5);
            no0vel_idx = find(total_mov_mm > threshold);
            vf = ftT.velFor{1}(summaryArray.Indices{trial});
            vf = vf(no0vel_idx); 
            vy = ftT.velYaw{1}(summaryArray.Indices{trial});
            vy = vy(no0vel_idx); 
            angle = ftT.cueAngle{1}(summaryArray.Indices{trial}); 
            angle = angle(no0vel_idx); 
            
            % shift HD to be relative to the estimated goal HD
            angle = wrapTo180(wrapTo180(angle)-wrapTo180(rad2deg(summaryArray.Goal(trial)))); 
            % convert yaw velocity to deg/sec
            vy = (vy/ (2*pi) ) * 360; 

            sum_mean = cell(3,1); 
            sum_mean{1} = zeros(length(edges_vy)-1,1);
            sum_mean{2} = zeros(length(edges_angle)-1,1);
            sum_mean{3} = zeros(length(edges_vf)-1,1);

            activity = bump_params.amp;
            activity = activity(summaryArray.Indices{trial}); 
            activity = activity(no0vel_idx);

            % bin neural activity by associated behavioural datapoints
            
            %vf 
            behaviour = vf; 
            [binnedActivity, centers_vf, ~] = binData(activity, behaviour, edges_vf);
            sum_mean{3} = binnedActivity;

           % vy 
             behaviour = vy; 
            [binnedActivity, centers_vy, ~] = binData(activity, behaviour, edges_vy);
            sum_mean{1} = binnedActivity; 

            %angle
            [binnedActivity, centers_angle, ~] = binData(activity, angle, edges_angle);
            sum_mean{2} = binnedActivity; 

            % keep count of how many trials came from the same fly
            trialCount = size(trial_vf{fly},1);
            trialCount = trialCount + 1; 

            trial_vf{fly}(trialCount,:) = sum_mean{3};
            trial_vy{fly}(trialCount,:) = sum_mean{1};
            trial_angle{fly}(trialCount,:) = sum_mean{2};
        
        catch
            disp(['folder ',folder,' failed'])
            failedFolders{fcount} = folder; 
            fcount = fcount + 1; 
            trialCount = size(trial_vf{fly},1);
            trialCount = trialCount + 1; 
            
            sum_mean{1} = nan(length(edges_vy)-1,1);
            sum_mean{2} = nan(length(edges_angle)-1,1);
            sum_mean{3} = nan(length(edges_vf)-1,1);
        
            trial_vf{fly}(trialCount,:) = sum_mean{3};
            trial_vy{fly}(trialCount,:) = sum_mean{1};
            trial_angle{fly}(trialCount,:) = sum_mean{2};
            
        end
    end  

    vfCount = 1; 
    vyCount = 1; 
    angleCount = 1; 
    % if num trials associated with each fly is > 0 calculate average
    % binned bump amplitude values across trials within each fly
    for t = 1:size(trial_vf,1)
        if ~isempty(trial_vf{t})
            flyAve_vf(vfCount,:) = squeeze(mean(trial_vf{t},1,'omitnan'));
            vfCount = vfCount + 1;
        end
        if ~isempty(trial_vy{t})
            flyAve_vy(vfCount,:) = squeeze(mean(trial_vy{t},1,'omitnan'));
            vyCount = vyCount + 1;
        end
        if ~isempty(trial_angle{t})
            flyAve_angle(angleCount,:) = squeeze(mean(trial_angle{t},1,'omitnan'));
            angleCount = angleCount + 1;
        end
    end
            
    % calculate the SEM for each relationship to plot
    SEM_vf = std(flyAve_vf,[],1,'omitnan') / sqrt(size(flyAve_vf,1));
    SEM_vy = std(flyAve_vy,[],1,'omitnan') / sqrt(size(flyAve_vy,1));
    SEM_angle = std(flyAve_angle,[],1,'omitnan') / sqrt(size(flyAve_angle,1));

    % average binned bump amplitude values across flies
    aveTrial_vf = mean(flyAve_vf,1,'omitnan');
    aveTrial_vy = mean(flyAve_vy,1,'omitnan');
    aveTrial_angle = mean(flyAve_angle,1,'omitnan');

    % plot data
    figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [aveTrial_vy(keepIndex) + SEM_vy(keepIndex)]; 
    SEMlow = [aveTrial_vy(keepIndex) - SEM_vy(keepIndex)];
    ax1 = subplot(3,1,1);
    hold on
    plot(centers_vy, aveTrial_vy,'b')
    patch([centers_vy(keepIndex) fliplr(centers_vy(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlabel('vy (deg/s)')
    ylabel('bump amplitude')

    keepIndex = ~isnan(SEM_angle);
    SEMhigh = [aveTrial_angle(keepIndex) + SEM_angle(keepIndex)]; 
    SEMlow = [aveTrial_angle(keepIndex) - SEM_angle(keepIndex)];
    ax2 = subplot(3,1,2);
    hold on
    plot(centers_angle, aveTrial_angle,'b')
    patch([centers_angle(keepIndex) fliplr(centers_angle(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ax1.YAxis.Exponent = 0;
    ax2.YAxis.Exponent = 0;
    xlabel('cue pos rel goal')
    ylabel('bump amplitude')

    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [aveTrial_vf(keepIndex) + SEM_vf(keepIndex)]; 
    SEMlow = [aveTrial_vf(keepIndex) - SEM_vf(keepIndex)];
    ax3 = subplot(3,1,3);
    hold on
    plot(centers_vf, aveTrial_vf,'b')
    patch([centers_vf(keepIndex) fliplr(centers_vf(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlabel('vf (mm/s)')
    ylabel('bump amplitude')
    ax3.YAxis.Exponent = 0;

end