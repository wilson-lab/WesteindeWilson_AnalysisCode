function summarizePFL3(summaryArray_all)
    % Averages data across flies & trials to make R-L activity vs
    % behaviour plots
    
    % INPUTS
    
    % summaryArray_all: table describing segmented PFL3 trials & each segment's estimated goal HD
       
    summaryArray = summaryArray_all(contains(summaryArray_all.Folder,'LAL'),:);     
    
    edges_vy = [-200:10:200];
    edges_vf = [-4:1:10];
    edges_angle = [-180:10:180];  
    
    trial_vyL = zeros(size(summaryArray,1),length(edges_vy)-1);
    trial_vyR = zeros(size(summaryArray,1),length(edges_vy)-1);
    trial_angleL = zeros(size(summaryArray,1),length(edges_angle)-1);
    trial_angleR = zeros(size(summaryArray,1),length(edges_angle)-1);
    trial_vy_LR = zeros(size(summaryArray,1),length(edges_vy)-1);
    trial_angle_LR = zeros(size(summaryArray,1),length(edges_angle)-1);
    speed_sum = [];

    threshold= 3;
    trial_count = 1;
    
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

            expID = get_expID(folder);
            expList = {expID};

            % Load imaging & behaviour data
            roiData = load_roi_data(expList, folder);

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

            try
                trial_roiData = roiData(roiData.trialNum == nTrial,:);
            catch 
                trial_roiData = [1;2];
            end

            sum_mean = cell(3,1);
            sum_mean{1} = zeros(length(edges_vy)-1,1);
            sum_mean{2} = zeros(length(edges_angle)-1,1);
            sum_mean{3} = zeros(length(edges_vf)-1,1);

            activityTable = ZData; 

            % bin neural activity by associated behavioural datapoints
            
            for roi = 1:size(trial_roiData,1)
                activity = activityTable.(3){roi}(summaryArray.Indices{trial});
                activity = activity(no0vel_idx);

                % vf 
                behaviour = vf; 
                [zscore, centers_vf, ~] = binData(activity, behaviour, edges_vf);
                sum_mean{3}(:,roi) = zscore;
                
                % vy 
                 behaviour = vy; 
                [zscore, centers_vy, ~] = binData(activity, behaviour, edges_vy);
                sum_mean{1}(:,roi) = zscore; 

                % angle
                [zscore, centers_angle, ~] = binData(activity, angle, edges_angle);
                sum_mean{2}(:,roi) = zscore; 

            end
            
            % calculate right left difference in PFL3 activity
            vf_LR = (sum_mean{3}(:,2) - sum_mean{3}(:,1));
            vy_LR = (sum_mean{1}(:,2) - sum_mean{1}(:,1));
            angle_LR = (sum_mean{2}(:,2) - sum_mean{2}(:,1)) ;

            % collect binned trial data
            trial_vfL(trial_count,:) = sum_mean{3}(:,1);
            trial_vfR(trial_count,:) = sum_mean{3}(:,2);
            trial_vf_LR(trial_count,:) = vf_LR;
            trial_vyL(trial_count,:) = sum_mean{1}(:,1);
            trial_vyR(trial_count,:) = sum_mean{1}(:,2);
            trial_vy_LR(trial_count,:) = vy_LR;
            trial_angleL(trial_count,:) = sum_mean{2}(:,1);
            trial_angleR(trial_count,:) = sum_mean{2}(:,2);
            trial_angle_LR(trial_count,:) = angle_LR; 


            trial_count = trial_count + 1; 
        catch
            disp(['folder ',folder,' failed'])

        end
    end
    
    % average binned activity values across trials
    aveTrial_vfL = mean(trial_vfL,1,'omitnan');
    aveTrial_vfR = mean(trial_vfR,1,'omitnan');
    aveTrial_vfLR = mean(trial_vf_LR,1,'omitnan');
    aveTrial_vyL = mean(trial_vyL,1,'omitnan');
    aveTrial_vyR = mean(trial_vyR,1,'omitnan');
    aveTrial_vyLR = mean(trial_vy_LR,1,'omitnan');
    aveTrial_angleL = mean(trial_angleL,1,'omitnan');
    aveTrial_angleR = mean(trial_angleR,1,'omitnan');
    aveTrial_angleLR = mean(trial_angle_LR,1,'omitnan');

    % calculate the SEM for each relationship to plot
    SEM_vfL = std(trial_vfL,[],1,'omitnan') / sqrt(size(trial_vfL,1));
    SEM_vfR = std(trial_vfR,[],1,'omitnan') / sqrt(size(trial_vfR,1));
    SEM_vfLR = std(trial_vf_LR,[],1,'omitnan') / sqrt(size(trial_vf_LR,1));
    SEM_vyL = std(trial_vyL,[],1,'omitnan') / sqrt(size(trial_vyL,1));
    SEM_vyR = std(trial_vyR,[],1,'omitnan') / sqrt(size(trial_vyR,1));
    SEM_vyLR = std(trial_vy_LR,[],1,'omitnan') / sqrt(size(trial_vy_LR,1));
    SEM_angleL = std(trial_angleL,[],1,'omitnan') / sqrt(size(trial_angleL,1));
    SEM_angleR = std(trial_angleR,[],1,'omitnan') / sqrt(size(trial_angleR,1));
    SEM_angleLR = std(trial_angle_LR,[],1,'omitnan') / sqrt(size(trial_angle_LR,1));


    % plot data
    data = 'Z';
    figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    keepIndex = ~isnan(SEM_vyL);
    SEMhigh = [aveTrial_vyL(keepIndex) + SEM_vyL(keepIndex)]; 
    SEMlow = [aveTrial_vyL(keepIndex) - SEM_vyL(keepIndex)];
    ax1 = subplot(2,1,1);
    hold on
    plot(centers_vy, aveTrial_vyL,'b')
    patch([centers_vy(keepIndex) fliplr(centers_vy(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')

    keepIndex = ~isnan(SEM_vyR);
    SEMhigh = [aveTrial_vyR(keepIndex) + SEM_vyR(keepIndex)]; 
    SEMlow = [aveTrial_vyR(keepIndex) - SEM_vyR(keepIndex)];
    subplot(2,1,1)
    plot(centers_vy, aveTrial_vyR,'r')
    patch([centers_vy(keepIndex) fliplr(centers_vy(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xlabel('vy (deg/s)')
    ylabel(data)

    keepIndex = ~isnan(SEM_angleL);
    SEMhigh = [aveTrial_angleL(keepIndex) + SEM_angleL(keepIndex)]; 
    SEMlow = [aveTrial_angleL(keepIndex) - SEM_angleL(keepIndex)];
    ax2 = subplot(2,1,2);
    hold on
    plot(centers_angle, aveTrial_angleL,'b')
    patch([centers_angle(keepIndex) fliplr(centers_angle(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')

    keepIndex = ~isnan(SEM_angleR);
    SEMhigh = [aveTrial_angleR(keepIndex) + SEM_angleR(keepIndex)]; 
    SEMlow = [aveTrial_angleR(keepIndex) - SEM_angleR(keepIndex)];
    subplot(2,1,2)
    plot(centers_angle, aveTrial_angleR,'r')
    patch([centers_angle(keepIndex) fliplr(centers_angle(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xlabel('cue pos relative to goal')
    ylabel(data)
    ax1.YAxis.Exponent = 0;
    ax2.YAxis.Exponent = 0;

    figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    keepIndex = ~isnan(SEM_vfL);
    SEMhigh = [aveTrial_vfL(keepIndex) + SEM_vfL(keepIndex)]; 
    SEMlow = [aveTrial_vfL(keepIndex) - SEM_vfL(keepIndex)];
    ax1 = subplot(2,1,1);
    hold on
    plot(centers_vf, aveTrial_vfL,'b')
    patch([centers_vf(keepIndex) fliplr(centers_vf(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')

    keepIndex = ~isnan(SEM_vfR);
    SEMhigh = [aveTrial_vfR(keepIndex) + SEM_vfR(keepIndex)]; 
    SEMlow = [aveTrial_vfR(keepIndex) - SEM_vfR(keepIndex)];
    subplot(2,1,1)
    plot(centers_vf, aveTrial_vfR,'r')
    patch([centers_vf(keepIndex) fliplr(centers_vf(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0,0],'FaceAlpha',0.25,'EdgeColor','none')
    xlabel('vf (mm/s)')
    ylabel(data)

    keepIndex = ~isnan(SEM_vfLR);
    SEMhigh = [aveTrial_vfLR(keepIndex) + SEM_vfLR(keepIndex)]; 
    SEMlow = [aveTrial_vfLR(keepIndex) - SEM_vfLR(keepIndex)];
    ax2 = subplot(2,1,2);
    hold on
    yline(0,'--r')
    plot(centers_vf, aveTrial_vfLR,'k')
    patch([centers_vf(keepIndex) fliplr(centers_vf(keepIndex))],[SEMhigh fliplr(SEMlow)],'k','FaceAlpha',0.25,'EdgeColor','none')
    xlabel('vf mm/s)')
    ylabel(data)
    ax1.YAxis.Exponent = 0;
    ax2.YAxis.Exponent = 0;

    figure();
    set(gcf,'color','w')
    set(gcf,'Renderer','painters')
    keepIndex = ~isnan(SEM_vyLR);
    SEMhigh = [aveTrial_vyLR(keepIndex) + SEM_vyLR(keepIndex)]; 
    SEMlow = [aveTrial_vyLR(keepIndex) - SEM_vyLR(keepIndex)];
    ax1 = subplot(2,1,1);
    hold on
    yline(0,'--r')
    plot(centers_vy, aveTrial_vyLR,'k')
    patch([centers_vy(keepIndex) fliplr(centers_vy(keepIndex))],[SEMhigh fliplr(SEMlow)],'k','FaceAlpha',0.25,'EdgeColor','none')
    xlabel('vy (deg/s)')
    ylabel(data)

    keepIndex = ~isnan(SEM_angleLR);
    SEMhigh = [aveTrial_angleLR(keepIndex) + SEM_angleLR(keepIndex)]; 
    SEMlow = [aveTrial_angleLR(keepIndex) - SEM_angleLR(keepIndex)];
    ax2 = subplot(2,1,2);
    hold on
    yline(0,'--r')
    plot(centers_angle, aveTrial_angleLR,'k')
    patch([centers_angle(keepIndex) fliplr(centers_angle(keepIndex))],[SEMhigh fliplr(SEMlow)],'k','FaceAlpha',0.25,'EdgeColor','none')
    xlabel('cue pos relative to goal')
    ylabel(data)
    ax1.YAxis.Exponent = 0;
    ax2.YAxis.Exponent = 0;
            
end