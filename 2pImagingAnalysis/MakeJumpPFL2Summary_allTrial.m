function MakeJumpPFL2Summary_allTrial(rootDir, startAtMov)       
% Summarize PFL2 responses following cue jumps 

% INPUTs

% rootDir: directory to look for experimental folders
% startAtMov: Boolean, find first + movement transition following the
% jump and collect window starting at that point. 
    
    % Collect experimental folders
     folders = get_folders(rootDir,1,0);
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','Z','bumpAmp','region'};
     jump_summary = cell2table(cell(1,13),'VariableNames',headers); 
     count = 1; 

    % Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    for ff = 1:folderNum
      
      % Get folder 
      folder = folders(ff).folder;
      
       
      % Load in data
        try 
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end
            
            if startAtMov
                [startIdx, ~] = FindMovementTransitions_2p(folder, 0);
            end

             processedData_dir = fullfile(folder,'processed_data');

            % Get data files
            expID = get_expID(folder);
            expList = {expID};

            % Load metadata 
            [~, trialMd] = load_metadata(expList, folder);

            % Load imaging data
            roiData = load_roi_data(expList, folder);
            try
                numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
            catch
                numTrials = 1; 
            end
            
            % iterate through trials
            for nTrial = 1:numTrials
                % get processed data
                clear ZData dffData ftT 
                bump_params = [];
                data_filelist = dir(processedData_dir);
                for files = 1:length(data_filelist)
                    if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                        load(fullfile(processedData_dir,data_filelist(files).name));
                    end
                end            
                load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))

                % Calculate full trial straightness metric
                [rho_fullTrial, ~] = CalculateAverageHeading(ftT,1.5, 'all');

                if ~ismembertol(rho_fullTrial,1,10^-10) % if rho = 1 something is wrong

                    % detect jumps & datapoints within a window around them
                    window = 15; 
                    [~, jump_array, ~] = detect_jumps(ftT, window, window,1,0);
                    % calculating rho over a longer window than I want to
                    % average & plot
                    [~, jump_array_rho, ~] = detect_jumps(ftT, 30, 30,1,0);
                    all_roi = []; 


                    % clean up jumps to remove those where the window of
                    % interest extends beyond time = 0 or time > length of
                    % trial
                    jumpIdx = [jump_array(1,1):jump_array(1,3)]';  
                    jumpIdx_rho = [jump_array_rho(1,1):jump_array_rho(1,3)]';

                    if jumpIdx(1) < 1 || jumpIdx_rho(1) < 1
                        jump_array(1,:) = []; 
                        jump_array_rho(1,:) = []; 
                    end

                    for jump = 1:size(jump_array,1)
                        jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
                        if jumpIdx(end) > size(ftT.velFor{1},1)
                            jump_array(jump,:) = []; 
                            jump_array_rho(jump,:) = []; 
                        end
                    end 

                    % iterate through remaining jumps
                    for jump = 1:size(jump_array,1)                   
                        if startAtMov % collect window after jump starting when fly started moving
                            % find first movement bout following jump (if any)
                            jumpmovStart = startIdx(startIdx > jump_array(jump,2) & startIdx < jump_array(jump,3));

                            totalMov = smoothdata(abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs((ftT.velYaw{1})*4.5),'gaussian',1);

                            if totalMov(jump_array(jump,2) + 1) >= 1.5
                                % was the fly's cumulative velocity above the movement threshold at the timepoint following the jump?
                                % if so collect jump data
                                jumpIdx = [jump_array(jump,1):jump_array(jump,3)]';   
                                jump_vf = ftT.velFor{1}(jumpIdx);
                                jump_vy = (ftT.velYaw{1}(jumpIdx)/(2*pi))* 360;
                                jump_vs = (ftT.velSide{1}(jumpIdx));
                                jump_angle = ftT.cueAngle{1}(jumpIdx); 

                                for roi = 1:size(ZData,1)
                                    all_roi(roi,:) = ZData.(3){roi}(jumpIdx);
                                end

                                jump_Z = mean(all_roi,1,'omitnan');
                                jump_bumpAmp = bump_params.amp(jumpIdx);

                            elseif ~isempty(jumpmovStart)
                                % if the fly wasn't moving right after the jump
                                % collect data starting from when the first movement bout started 
                                jumpmovStart = min(jumpmovStart);
                                jumpEnd = jumpmovStart + window * 60;
                                preJump = [jump_array(jump,1):jump_array(jump,2)];
                                postJump = [jumpmovStart:jumpEnd-1];   
                                jumpIdx = [preJump,postJump];
                                jump_vf = ftT.velFor{1}(jumpIdx);
                                jump_vy = (ftT.velYaw{1}(jumpIdx)/(2*pi))* 360;
                                jump_vs = (ftT.velSide{1}(jumpIdx));
                                jump_angle = ftT.cueAngle{1}(jumpIdx);

                                for roi = 1:size(ZData,1)
                                    all_roi(roi,:) = ZData.(3){roi}(jumpIdx);
                                end

                                jump_Z = mean(all_roi,1,'omitnan');

                                jump_bumpAmp = bump_params.amp(jumpIdx);
                            else
                                % if the fly didn't move after the jump skip
                                % the jump
                                jumpIdx = [];  
                                jump_vf = [];
                                jump_vy = [];
                                jump_vs = [];
                                jump_angle = [];
                                jump_Z = []; 
                                jump_bumpAmp = []; 
                            end
                        else 
                            jumpIdx = [jump_array(jump,1):jump_array(jump,3)]';   
                            jump_vf = ftT.velFor{1}(jumpIdx);
                            jump_vy = (ftT.velYaw{1}(jumpIdx)/(2*pi))* 360;
                            jump_vs = (ftT.velSide{1}(jumpIdx));
                            jump_angle = ftT.cueAngle{1}(jumpIdx); 
                            jumpmovStart = [];

                            for roi = 1:size(ZData,1)
                                    all_roi(roi,:) = ZData.(3){roi}(jumpIdx);
                            end

                            jump_Z = mean(all_roi,1,'omitnan');
                            jump_bumpAmp = bump_params.amp(jumpIdx);
                        end


                        % gather info from the jump data
                        jumpSize = jump_array(jump,4);
                        preJumpIdx = [jump_array(jump,1):jump_array(jump,2)-1];
                        postJumpIdx = [jump_array(jump,2):jump_array(jump,3)-1];
                        preJumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,2)-1];

                        timeMov = (sum(abs(ftT.velFor{1}(preJumpIdx)) + abs(ftT.velSide{1}(preJumpIdx)) + abs((ftT.velYaw{1}(preJumpIdx))*4.5) > 3))/60; % seconds
                        timeMov_post = (sum(abs(ftT.velFor{1}(postJumpIdx)) + abs(ftT.velSide{1}(postJumpIdx)) + abs((ftT.velYaw{1}(postJumpIdx))*4.5) > 3))/60; % seconds
                        [rho, ~] = CalculateAverageHeading(ftT,1.5, preJumpIdx_rho);
                        [~, preJumpAngle] = CalculateAverageHeading(ftT,0, preJumpIdx);

                        postJumpIdx = jump_array(jump,2)+1;
                        sizePostJump = size(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)));
                        preJump = ones(sizePostJump) * preJumpAngle;
                        cueDiff = angdiff(preJump,deg2rad(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)))); 


                        % if fly returned cue position within a minimum range
                        % of its pre jump position within the time window the 
                        % jump is categorized as corrected 
                        if jumpSize == 180
                            correct = sum(abs(cueDiff) < deg2rad(60)); 
                        else
                            correct = sum(abs(cueDiff) < deg2rad(30));
                        end

                        % add jump info to summary table
                        jump_summary.folder(count) = {folder}; 
                        jump_summary.trial(count) = {nTrial}; 
                        jump_summary.jumpSize(count) = {jumpSize}; 
                        jump_summary.jumpIdx(count) = {jumpIdx};
                        jump_summary.velFor(count) = {jump_vf};
                        jump_summary.velYaw(count) = {jump_vy};
                        jump_summary.velSide(count) = {jump_vs}; 
                        jump_summary.cueAngle(count) = {jump_angle};
                        jump_summary.timeMov(count) = {timeMov}; 
                        jump_summary.Z(count) = {jump_Z}; 
                        jump_summary.bumpAmp(count) = {jump_bumpAmp}; 
                        jump_summary.postJumpMoveIdx(count) = {jumpmovStart};
                        isPB = regexp(folder,'PB');
                        isFB = regexp(folder,'FB');

                        if ~isempty(isPB)
                            jump_summary.region(count) = {'PB'};
                        elseif ~isempty(isFB)
                            jump_summary.region(count) = {'FB'};
                        else
                            jump_summary.region(count) = {'LAL'};
                        end


                        % set jump category value based on overall movement,
                        % straightness prior to the jump, and whether the jump
                        % was corrected for. 
                        if  timeMov >= 1 && timeMov_post > 0 && correct >= 1  && rho >= 0.88
                            jump_summary.corrected(count) = {1};
                        elseif correct < 1 && timeMov >= 5 && rho >= 0.88
                            jump_summary.corrected(count) = {2}; 
                            figure(1);clf;plot(jump_angle)
                        elseif  timeMov >= 3 && timeMov_post > 0 && rho < 0.8   && correct < 1
                            jump_summary.corrected(count) = {0}; 

                        else
                            jump_summary.corrected(count) = {3}; 
                        end

                        count = count + 1; 

                    end
                end
            end
        catch
            disp([folder, ' failed'])
        end
    end  
    % organize data into appropriate categories
    Sum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 1,:); 
    Sum180_not = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 0,:); 
    Sum90 = jump_summary((cell2mat(jump_summary.jumpSize)) == 90& cell2mat(jump_summary.corrected) == 1,:);
    Sum90_not = jump_summary((cell2mat(jump_summary.jumpSize)) == 90& cell2mat(jump_summary.corrected) == 0,:);
    Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 1,:);
    Summ90_not = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 0,:);

    %% Create jump summary plots
    % smooth activity & behaviour data for visualization
    smoothWin    = 150;  
    smoothWin_v = 200; 
    jumpTime = [-window+(1/60):1/60:window+(1/60)]; 
    
    % 180 jump plot
    ZSum = zeros(size(Sum180,1),size(jump_vf,1));
    bAmpSum = zeros(size(Sum180,1),size(jump_vf,1));
    vfSum = zeros(size(Sum180,1),size(jump_vf,1));
    vySum = zeros(size(Sum180,1),size(jump_vf,1));
    figure();
    sgtitle('180 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    ax1 = subplot(4,1,1); 
    hold on
    ax2 = subplot(4,1,3);
    hold on
    ax3 = subplot(4,1,4); 
    hold on
    ax4 = subplot(4,1,2); 
    hold on

    for jump = 1:size(Sum180,1)  
        if ~isempty(Sum180.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(Sum180.Z{jump},'loess',smoothWin);
            bAmpSum(jump,:) = smoothdata(Sum180.bumpAmp{jump},'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum180.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum180.velYaw{jump}'),'loess',smoothWin_v);
        end
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(ax1,jumpTime,mean(ZSum,'omitnan'),'b','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(ax1,'Z')
    xlabel(ax1,'Time (s)')
    xlim(ax1,[min(jumpTime),max(jumpTime)])
    
    keepIndex = ~isnan(SEM_bAmp);
    SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
    SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
    plot(ax4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
    patch(ax4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax4,[min(jumpTime),max(jumpTime)])
    ylabel(ax4,'bump amp')
    xlabel(ax4,'Time (s)')
    
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vf(keepIndex)];
    plot(ax2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
    patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax2,[min(jumpTime),max(jumpTime)])
    ylabel(ax2,'forward velocity (mm/s)')
    xlabel(ax2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
    plot(ax3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
    patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax3,[min(jumpTime),max(jumpTime)])
    ylabel(ax3,'rotational speed (deg/s)')
    xlabel(ax3,'Time (s)')
    
    ZSum = zeros(size(Sum180_not,1),size(jump_vf,1));
    bAmpSum = zeros(size(Sum180_not,1),size(jump_vf,1));
    vfSum = zeros(size(Sum180_not,1),size(jump_vf,1));
    vySum = zeros(size(Sum180_not,1),size(jump_vf,1));

    for jump = 1:size(Sum180_not,1) 
        if ~isempty(Sum180_not.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(Sum180_not.Z{jump},'loess',smoothWin);
            bAmpSum(jump,:) = smoothdata(Sum180_not.bumpAmp{jump},'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum180_not.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum180_not.velYaw{jump})','loess',smoothWin_v);
        end
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(ax1,jumpTime,mean(ZSum,'omitnan'),'k','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(ax1,'Z')
    xlabel(ax1,'Time (s)')
    xlim(ax1,[min(jumpTime),max(jumpTime)])
    
    keepIndex = ~isnan(SEM_bAmp);
    SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
    SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
    plot(ax4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(ax4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax4,[min(jumpTime),max(jumpTime)])
    ylabel(ax4,'bump amp')
    xlabel(ax4,'Time (s)')
    
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vf(keepIndex)];
    plot(ax2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax2,[min(jumpTime),max(jumpTime)])
    ylabel(ax2,'forward velocity (mm/s)')
    xlabel(ax2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
    plot(ax3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax3,[min(jumpTime),max(jumpTime)])
    ylabel(ax3,'rotational speed (deg/s)')
    xlabel(ax3,'Time (s)')

    linkaxes([ax1,ax2,ax3,ax4],'x') 
    
    % +/-90 jump plot
    vfSum = zeros(size(Sum90,1),size(jump_vf,1));
    vySum = zeros(size(Sum90,1),size(jump_vf,1));
    ZSum = zeros(size(Sum90,1),size(jump_vf,1));
    bAmpSum = zeros(size(Sum90,1),size(jump_vf,1));
    figure();
    sgtitle('90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    bx1 = subplot(4,1,1);
    hold on
    bx2 = subplot(4,1,3);
    hold on
    bx3 = subplot(4,1,4); 
    hold on
    bx4 = subplot(4,1,2); 
    hold on
    for jump = 1:size(Sum90,1) 
        if ~isempty(Sum90.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(Sum90.Z{jump},'loess',smoothWin);
            bAmpSum(jump,:) = smoothdata(Sum90.bumpAmp{jump},'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum90.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum90.velYaw{jump})','loess',smoothWin_v);
        end
    end
    
        startCount = jump;
    for jump = 1:size(Summ90,1)
        if ~isempty(Summ90.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(Summ90.Z{jump},'loess',smoothWin);
            bAmpSum(jump,:) = smoothdata(Summ90.bumpAmp{jump},'loess',smoothWin);
            vfSum(startCount,:) = smoothdata(Summ90.velFor{jump}','loess',smoothWin_v); 
            vySum(startCount,:) = smoothdata(abs(Summ90.velYaw{jump})','loess',smoothWin_v);
            startCount = startCount + 1; 
        end
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(bx1,jumpTime,mean(ZSum,'omitnan'),'b','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(bx1,'Z')
    xlabel(bx1,'Time (s)')
    xlim(bx1,[min(jumpTime),max(jumpTime)])
    
    keepIndex = ~isnan(SEM_bAmp);
    SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
    SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
    plot(bx4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
    patch(bx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx4,[min(jumpTime),max(jumpTime)])
    ylabel(bx4,'bump amp')
    xlabel(bx4,'Time (s)')
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex)),'gaussian',5) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex)),'gaussian',5) - SEM_vf(keepIndex)];
    plot(bx2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx2,[min(jumpTime),max(jumpTime)])
    ylabel(bx2,'forward velocity (mm/s)')
    xlabel(bx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
    plot(bx3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx3,[min(jumpTime),max(jumpTime)])
    ylabel(bx3,'rotational velocity (deg/s)')
    xlabel(bx3,'Time (s)')
    

    for jump = 1:size(Sum90_not,1) 
        if ~isempty(Sum90_not.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(Sum90_not.Z{jump},'loess',smoothWin);
            bAmpSum(jump,:) = smoothdata(Sum90_not.bumpAmp{jump},'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum90_not.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum90_not.velYaw{jump})','loess',smoothWin_v);
        end
    end
    
        startCount = jump;
    for jump = 1:size(Summ90_not,1)
        if ~isempty(Summ90_not.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(Summ90_not.Z{jump},'loess',smoothWin);
            bAmpSum(jump,:) = smoothdata(Summ90_not.bumpAmp{jump},'loess',smoothWin);
            vfSum(startCount,:) = smoothdata(Summ90_not.velFor{jump}','loess',smoothWin_v); 
            vySum(startCount,:) = smoothdata(abs(Summ90_not.velYaw{jump})','loess',smoothWin_v);
            startCount = startCount + 1; 
        end
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(bx1,jumpTime,mean(ZSum,'omitnan'),'k','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(bx1,'Z')
    xlabel(bx1,'Time (s)')
    xlim(bx1,[min(jumpTime),max(jumpTime)])
    
    keepIndex = ~isnan(SEM_bAmp);
    SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
    SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
    plot(bx4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(bx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx4,[min(jumpTime),max(jumpTime)])
    ylabel(bx4,'bump amp')
    xlabel(bx4,'Time (s)')
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex)),'gaussian',5) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex)),'gaussian',5) - SEM_vf(keepIndex)];
    plot(bx2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx2,[min(jumpTime),max(jumpTime)])
    ylabel(bx2,'forward velocity (mm/s)')
    xlabel(bx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
    plot(bx3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx3,[min(jumpTime),max(jumpTime)])
    ylabel(bx3,'rotational velocity (deg/s)')
    xlabel(bx3,'Time (s)')
    
    linkaxes([bx1,bx2,bx3,bx4],'x')

end