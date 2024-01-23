function summarizePFL2_meanActivity(summaryArray)       
    % Averages data across flies & trials to make mean activity vs
    % behaviour plots
    
    % INPUTS
    
    % summaryArray_all: table describing segmented PFL2 trials & each segment's estimated goal HD


    edges_vy = [-200:10:200];
    edges_vf = [-4:1:10];
    edges_angle = [-180:10:180]; 
    
    threshold=3;
    trial_count = 1;
    
    trial_vy = [];
    trial_angle = [];
    trial_vf = []; 
    speed_sum = [];

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

        % Load imaging & behaviour data
        clear ZData ftT
        processedData_dir = fullfile(folder,'processed_data');
        nTrial = summaryArray.numTrial(trial);
        data_filelist = dir(processedData_dir);
        for files = 1:length(data_filelist)
            if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                load(fullfile(processedData_dir,data_filelist(files).name));
            end            
        end
        load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))
        
        all_roi = [];
        for roi = 1:size(ZData,1)
            all_roi(roi,:) = ZData.Z{roi};
        end
        
        % calculate mean population activity for each timepoint
        Z_ave = mean(all_roi,1,'omitnan');

        % Remove idx where the fly isn't moving
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

        activity = Z_ave;
        activity = activity(summaryArray.Indices{trial}); 
        activity = activity(no0vel_idx);
        
        % bin neural activity by associated behavioural datapoints

        %vf 
        behaviour = vf; 
        [zscore, centers_vf, ~] = binData(activity, behaviour, edges_vf);
        sum_mean{3} = zscore;

       % vy 
         behaviour = vy; 
        [zscore, centers_vy, ~] = binData(activity, behaviour, edges_vy);
        sum_mean{1} = zscore; 

        %angle
        [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
        sum_mean{2} = zscore; 
        
        trial_vf(trial_count,:) = sum_mean{3};
        trial_vy(trial_count,:) = sum_mean{1};
        trial_angle(trial_count,:) = sum_mean{2};

        trial_count = trial_count + 1; 
        
        catch
            disp(['folder ',folder,' failed'])
        end
    end
    
    % calculate the SEM for each relationship to plot
    SEM_vf = std(trial_vf,[],1,'omitnan') / sqrt(size(trial_vf,1));
    SEM_vy = std(trial_vy,[],1,'omitnan') / sqrt(size(trial_vy,1));
    SEM_angle = std(trial_angle,[],1,'omitnan') / sqrt(size(trial_angle,1));

    % average binned bump amplitude values across trials
    aveTrial_vf = mean(trial_vf,1,'omitnan');
    aveTrial_vy = mean(trial_vy,1,'omitnan');
    aveTrial_angle = mean(trial_angle,1,'omitnan');


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
    ylabel('mean z-scored df/f')

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
    ylabel('mean z-scored df/f')

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
    ylabel('mean z-scored df/f')
    ax3.YAxis.Exponent = 0;
end