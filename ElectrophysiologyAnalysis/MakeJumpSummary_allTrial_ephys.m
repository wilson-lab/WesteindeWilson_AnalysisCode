function MakeJumpSummary_allTrial_ephys(rootDir)       
% Makes summary plots between neural activity & behaviour around cue jumps
 
    % gather experiment folders
    folders = get_folders_ephys(rootDir);

    all_folders = []; 
    for f = 1:length(folders)
    all_folders{f} = (folders(f).folder);
    end
    all_folders = all_folders';

    [flies, ~, uniqueFlies] = get_flies_ephys(all_folders);
        
    % setup summary table
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','prefHead','preJump','postJump','activity'};
     jump_summary = cell2table(cell(1,14),'VariableNames',headers); 
     count = 1; 

    % Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    for ff = 1:folderNum      
      % Get folder information
      folder = folders(ff).folder;
      cell_count = find(uniqueFlies == flies(ff));
             
        % Load in fictrac & ephys data
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        processedDir = fullfile(folder,'processedData');
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        load(fullfile(processedDir,'pro_trialData.mat'))
        bData = pro_behaviourData{1}; 
        tData = pro_trialData{1}; 

        % window (s) around jumps to average & plot
        window = 15;
        [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, window,window, 1000);
        % window (s) around jumps to evaluate rho
        [jump_array_rho, ~, ~] = detect_jumps_ephys(bData.frY, 30,30, 1000);
        % cumulative speed
        total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5);
        % remove indicies where fly didn't move
        no0vel_idx = find(total_mov_mm > 3);
        angle = bData.angle(no0vel_idx); 

        activity = tData.smooth_Vm;
        activity = activity(no0vel_idx);

        % bin membrane potential by HD, create estimate preferred HD of
        % cell 
        edges_angle = [-180:30:180];
        [angleBins, centers_angle, ~] = binData(activity, angle, edges_angle);
        % lightly smooth binned data to better estimate curve peak
        prefHead_smooth = smoothdata(angleBins,'gaussian',3); 
        prefHead = centers_angle(prefHead_smooth == max(prefHead_smooth));

        jumpIdx = [jump_array(1,1):jump_array(1,3)]';  
        jumpIdx_rho = [jump_array_rho(1,1):jump_array_rho(1,3)]';

        % discard jumps if window start idx is < 1
        if jumpIdx(1) < 1 || jumpIdx_rho(1) < 1
            jump_array(1,:) = []; 
            jump_array_rho(1,:) = []; 
        end

        % discard jumps if window end is > trial length
        for jump = 1:size(jump_array,1)
            jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
            if jumpIdx(end) > size(bData,1)
                jump_array(jump,:) = []; 
                jump_array_rho(jump,:) = []; 
            end
        end 

        % iterate through remaining jumps
        for jump = 1:size(jump_array,1)
            jumpSize = jump_array(jump,4);
            jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
            jump_vf = bData.vel_for(jumpIdx);
            jump_vy = bData.vel_yaw(jumpIdx);
            jump_vs = bData.vel_side(jumpIdx);
            jump_angle = bData.angle(jumpIdx);

            preJumpIdx = [jump_array(jump,1):jump_array(jump,2)-1];
            preJumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,2)-1];

            % amount of movement prior to jump
            timeMov = (sum(abs(bData.vel_for(preJumpIdx)) + abs(bData.vel_side(preJumpIdx)) + abs(deg2rad(bData.vel_yaw(preJumpIdx))*4.5) > 3))/1000; % seconds

            % average baseline Vm pre jump
            jump_Vm = tData.smoothVm(jumpIdx) - mean(tData.smoothVm(jump_array(jump,2) - (2*1000):jump_array(jump,2)-1)); 
            % baseline corrected pre jump Vm
            changeJump_pre =  mean(tData.smoothVm(jump_array(jump,2) - (1*1000):jump_array(jump,2)-1),'omitnan');
            % baseline corrected post jump Vm
            changeJump_post = mean(tData.smoothVm(jump_array(jump,2)+1:jump_array(jump,2) + (1*1000)),'omitnan');

            % calculate pre jump path straightness 
            [rho, ~] = CalculateAverageHeading_ephys(bData,1, preJumpIdx_rho);
            % calculate pre jump angle 
            [~, preJumpAngle] = CalculateAverageHeading_ephys(bData,1, preJumpIdx);

            % set up data arrays
            postJumpIdx = jump_array(jump,2)+1;
            sizePostJump = size(bData.angle(postJumpIdx:jump_array(jump,3)));
            preJump = ones(sizePostJump) * preJumpAngle;

            % calculate difference between pre jump & post jump cue angle
            cueDiff = angdiff(preJump,deg2rad(bData.angle(postJumpIdx:jump_array(jump,3)))); 

            % decide if jump was 'corrected for' or not by how close to pre
            % jump position the post jump cue was returned to
            if jumpSize == 180
                correct = sum(abs(cueDiff) < deg2rad(60)); 
            else
                correct = sum(abs(cueDiff) < deg2rad(30));
            end


            % collect relevent info 
            jump_summary.folder(count) = {folder};
            jump_summary.trial(count) = {1}; 
            jump_summary.jumpSize(count) = {jumpSize}; 
            jump_summary.jumpIdx(count) = {jumpIdx};
            jump_summary.velFor(count) = {jump_vf};
            jump_summary.velYaw(count) = {jump_vy};
            jump_summary.velSide(count) = {jump_vs}; 
            jump_summary.cueAngle(count) = {jump_angle};
            jump_summary.Vm(count) = {jump_Vm};
            jump_summary.timeMov(count) = {timeMov}; 
            jump_summary.preJump(count) = {bData.angle(jump_array(jump,2)-1)}; 
            jump_summary.postJump(count) = {bData.angle(jump_array(jump,2)+1)}; 
            jump_summary.prefHead(count) = {prefHead}; 
            jump_summary.activity(count) = {jump_Vm};
            jump_summary.rho(count) = {rho};
            jump_summary.correct(count) = {correct};
            jump_summary.changeJump_pre(count) = {changeJump_pre};
            jump_summary.changeJump_post(count) = {changeJump_post};
            jump_summary.cell(count) = {cell_count};

            % classify jumps 
            if  timeMov >= 3 && correct >= 5 && rho >= 0.88
                jump_summary.corrected(count) = {1}; 
            elseif timeMov >= 3 && rho < 0.8 && correct < 5
                jump_summary.corrected(count) = {0}; 
            elseif timeMov < 5
                jump_summary.corrected(count) = {2}; 
            else
                jump_summary.corrected(count) = {3}; 
            end

            count = count + 1; 

        end
    end 
    
    % separate jumps into relevent categories
    allCorr = jump_summary((cell2mat(jump_summary.corrected) == 1),:); 
    all180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180,:); 
    all90 = jump_summary(abs(cell2mat(jump_summary.jumpSize)) == 90 ,:);
    
    Sum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & (cell2mat(jump_summary.corrected) == 1),:); 
    Sum180_not = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 0,:); 
    Sum90 = jump_summary((cell2mat(jump_summary.jumpSize)) == 90  & (cell2mat(jump_summary.corrected) == 1),:);
    Sum90_not = jump_summary((cell2mat(jump_summary.jumpSize)) == 90 & cell2mat(jump_summary.corrected) == 0,:);
    Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90  & (cell2mat(jump_summary.corrected) == 1),:);
    Summ90_not = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 0,:);
    %% visualization plots
    
    % plot jump statistics 
    
    corrPreDiff = angdiff(deg2rad(cell2mat(Sum180.prefHead)),deg2rad(cell2mat(Sum180.preJump))); % preJump - prefHead + = clockwise - = counterclockwise
    corrPostDiff = angdiff(deg2rad(cell2mat(Sum180.prefHead)),deg2rad(cell2mat(Sum180.postJump))); % postJump - prefHead + = clockwise - = counterclockwise
    corrDiff = (rad2deg(angdiff(abs(corrPreDiff),abs(corrPostDiff)))); % PostDiff - PreDiff = 0 distance didn't change, >0 distance increased, <0 distance decreased
    
    ncorrPreDiff = angdiff(deg2rad(cell2mat(Sum180_not.prefHead)),deg2rad(cell2mat(Sum180_not.preJump))); % preJump - prefHead + = clockwise - = counterclockwise
    ncorrPostDiff = angdiff(deg2rad(cell2mat(Sum180_not.prefHead)),deg2rad(cell2mat(Sum180_not.postJump))); % postJump - prefHead + = clockwise - = counterclockwise
    ncorrDiff = (rad2deg(angdiff(abs(ncorrPreDiff),abs(ncorrPostDiff)))); % PostDiff - PreDiff = 0 distance didn't change, >0 distance increased, <0 distance decreased
    
    jump_change_all = [corrDiff;ncorrDiff];
    g1 = repmat({'High'},length(corrDiff),1);
    g2 = repmat({'Low'},length(ncorrDiff),1);

    base2 = ones(size(ncorrDiff))*1.4;  
    base1 = ones(size(corrDiff))*0.6;
    
    cells = unique([cell2mat(Sum180_not.cell);cell2mat(Sum180.cell)]);
    nCol = length(cells);    
    colormap = cbrewer2('RdBu',nCol);
    cell_num_corr = []; 
    cellColours_corr = []; 
     for jump = 1:size(Sum180,1)
        cell_num_corr(jump) = find(Sum180.cell{jump} == cells); 
        cellColours_corr(jump,:) = colormap(cell_num_corr(jump),:); 
     end
    
    cell_num_ncorr = []; 
    cellColours_ncorr = []; 
     for jump = 1:size(Sum180_not,1)
        cell_num_ncorr(jump) = find(Sum180_not.cell{jump} == cells); 
        cellColours_ncorr(jump,:) = colormap(cell_num_ncorr(jump),:); 
    end
    
    figure();hold on;
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    scatter(base1,corrDiff,[],cellColours_corr,'filled')
    scatter(base2,ncorrDiff,[],cellColours_ncorr,'filled')
    plot([base1(1),base2(1)],[mean(corrDiff),mean(ncorrDiff)],'o-k','MarkerFaceColor','k')
    xlim([0 2])
    [~,p] = ttest2(corrDiff,ncorrDiff); 
    text(1.6,10,['p = ',num2str(p)])
    

    [p,stats] = vartestn(jump_change_all,[g1;g2],'TestType','BrownForsythe');
    
    % 180 jump plots
    
    xlimMin = -2; 
    xlimMax = 2; 
    smoothWin    = 250;  
    smoothWin_v = 300;    
    jumpTime = [-window+(1/1000):1/1000:window+(1/1000)]; 
    ZSum = zeros(size(Sum180,1),size(Sum180.jumpIdx{1} ,1));
    vfSum = zeros(size(Sum180,1),size(Sum180.jumpIdx{1} ,1));
    vySum = zeros(size(Sum180,1),size(Sum180.jumpIdx{1} ,1));
    vsSum = zeros(size(Sum180,1),size(Sum180.jumpIdx{1} ,1));
    jump_change180 = [];
    figure();
    sgtitle('180 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    ax1 = subplot(3,1,1); 
    hold on
    ax2 = subplot(3,1,2);
    hold on
    ax3 = subplot(3,1,3); 
    hold on


    for jump = 1:size(Sum180,1)  
        if ~isempty(Sum180.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(abs(Sum180.Vm{jump}),'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum180.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum180.velYaw{jump})','loess',smoothWin_v);
            vsSum(jump,:) = smoothdata(Sum180.velSide{jump}','loess',smoothWin_v);
            jump_change180(jump) = (Sum180.changeJump_post{jump} - Sum180.changeJump_pre{jump});
        end
    end  
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    mid = (size(ZSum,2)-1)/2;
    
    meanZSum = mean(ZSum,1);   
    meanPreJump = mean(meanZSum(mid - 1000:mid - 1)); 
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum,'omitnan') -  meanPreJump+ SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum,'omitnan') -  meanPreJump - SEM_Z(keepIndex)];
    plot(ax1,jumpTime,mean(ZSum,'omitnan') -  meanPreJump,'b','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(ax1,'Vm')
    xlabel(ax1,'Time (s)')
    xlim(ax1,[xlimMin,xlimMax])
    
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [mean(vfSum(:,keepIndex),'omitnan') + SEM_vf(keepIndex)]; 
    SEMlow = [mean(vfSum(:,keepIndex),'omitnan') - SEM_vf(keepIndex)];
    plot(ax2,jumpTime,mean(vfSum,'omitnan'),'b','LineWidth',1.5)
    patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax2,[xlimMin,xlimMax])
    ylabel(ax2,'forward velocity (mm/s)')
    xlabel(ax2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(vySum(:,keepIndex),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(vySum(:,keepIndex),'omitnan') - SEM_vy(keepIndex)];
    plot(ax3,jumpTime,mean(vySum,'omitnan'),'b','LineWidth',1.5)
    patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax3,[xlimMin,xlimMax])
    ylabel(ax3,'rotational speed (deg/s)')
    xlabel(ax3,'Time (s)')
    
    ZSum = zeros(size(Sum180_not,1),size(Sum180.jumpIdx{1} ,1));
    vfSum = zeros(size(Sum180_not,1),size(Sum180.jumpIdx{1} ,1));
    vySum = zeros(size(Sum180_not,1),size(Sum180.jumpIdx{1} ,1));
    vsSum = zeros(size(Sum180_not,1),size(Sum180.jumpIdx{1} ,1));
    jump_change180not = [];

    for jump = 1:size(Sum180_not,1) 
        if ~isempty(Sum180_not.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(abs(Sum180_not.Vm{jump}),'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum180_not.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum180_not.velYaw{jump})','loess',smoothWin_v);
            vsSum(jump,:) = abs(Sum180_not.velSide{jump})';
            jump_change180not(jump) = (Sum180_not.changeJump_post{jump} - Sum180_not.changeJump_pre{jump});
        end
    end

    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    meanZSum = mean(ZSum,1);   
    meanPreJump =  mean(meanZSum(mid - 1000:mid - 1)); 
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum,'omitnan') -  meanPreJump + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum,'omitnan') -  meanPreJump - SEM_Z(keepIndex)];
    plot(ax1,jumpTime,mean(ZSum,'omitnan') -  meanPreJump,'k','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(ax1,'Vm')
    xlabel(ax1,'Time (s)')
    xlim(ax1,[xlimMin,xlimMax])     
 
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [mean(vfSum(:,keepIndex),'omitnan') + SEM_vf(keepIndex)]; 
    SEMlow = [mean(vfSum(:,keepIndex),'omitnan') - SEM_vf(keepIndex)];
    plot(ax2,jumpTime,mean(vfSum,'omitnan'),'k','LineWidth',1.5)
    patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax2,[xlimMin,xlimMax])
    ylabel(ax2,'forward velocity (mm/s)')
    xlabel(ax2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(vySum(:,keepIndex),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(vySum(:,keepIndex),'omitnan') - SEM_vy(keepIndex)];
    plot(ax3,jumpTime,mean(vySum,'omitnan'),'k','LineWidth',1.5)
    patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(ax3,[xlimMin,xlimMax])
    ylabel(ax3,'rotational speed (deg/s)')
    xlabel(ax3,'Time (s)')
    
    linkaxes([ax1,ax2,ax3],'x')
     
    jump_change_all = [jump_change180';jump_change180not'];
    g1 = repmat({'High'},length(jump_change180),1);
    g2 = repmat({'Low'},length(jump_change180not),1);

    base2 = ones(size(jump_change180not))*1.4;  
    base1 = ones(size(jump_change180))*0.6;
    
    cells = unique([cell2mat(Sum180_not.cell);cell2mat(Sum180.cell)]);
    nCol = length(cells);    
    colormap = cbrewer2('RdBu',nCol);
    cell_num_corr = []; 
    cellColours_corr = []; 
     for jump = 1:size(Sum180,1)
        cell_num_corr(jump) = find(Sum180.cell{jump} == cells); 
        cellColours_corr(jump,:) = colormap(cell_num_corr(jump),:); 
     end
    
    cell_num_ncorr = []; 
    cellColours_ncorr = []; 
     for jump = 1:size(Sum180_not,1)
        cell_num_ncorr(jump) = find(Sum180_not.cell{jump} == cells); 
        cellColours_ncorr(jump,:) = colormap(cell_num_ncorr(jump),:); 
    end
    
    figure();hold on;
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    scatter(base1,jump_change180,[],cellColours_corr,'filled')
    scatter(base2,jump_change180not,[],cellColours_ncorr,'filled')
    xlim([0 2])

     [p,stats] = vartestn(jump_change_all,[g1;g2],'TestType','BrownForsythe');

% 90 jump plots        
    vfSum = zeros(size(Sum90,1),size(Sum180.jumpIdx{1} ,1));
    vySum = zeros(size(Sum90,1),size(Sum180.jumpIdx{1} ,1));
    vsSum = zeros(size(Sum90,1),size(Sum180.jumpIdx{1} ,1));
    ZSum = zeros(size(Sum90,1),size(Sum180.jumpIdx{1} ,1));
    jump_change90 = [];
    figure();
    sgtitle('90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    bx1 = subplot(3,1,1);
    hold on
    bx2 = subplot(3,1,2);
    hold on
    bx3 = subplot(3,1,3); 
    hold on
    
    jump_change90 = [];

    for jump = 1:size(Sum90,1) 
        if ~isempty(Sum90.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(abs(Sum90.Vm{jump}),'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum90.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum90.velYaw{jump})','loess',smoothWin_v);
            vsSum(jump,:) = smoothdata(Sum90.velSide{jump}','loess',smoothWin_v);
            jump_change90(jump) = (Sum90.changeJump_post{jump} - Sum90.changeJump_pre{jump});
        end
    end
    
            startCount = jump;
    for jump = 1:size(Summ90,1)
        if ~isempty(Summ90.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(abs(Summ90.Vm{jump}),'loess',smoothWin);
            vfSum(startCount,:) = smoothdata(Summ90.velFor{jump}','loess',smoothWin_v); 
            vySum(startCount,:) = smoothdata(abs(Summ90.velYaw{jump})','loess',smoothWin_v);
            vsSum(startCount,:) = smoothdata(Summ90.velSide{jump}','loess',smoothWin_v);
            jump_change90(startCount) = (Summ90.changeJump_post{jump} - Summ90.changeJump_pre{jump});
            startCount = startCount + 1; 
        end
    end
    
    meanZSum = mean(ZSum,1);   
    meanPreJump = mean(meanZSum(mid - 1000:mid - 1)); 
    
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan')- meanPreJump + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') -  meanPreJump - SEM_Z(keepIndex)];
    plot(bx1,jumpTime,mean(ZSum,'omitnan') - meanPreJump,'b','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(bx1,'|Z|')
    xlabel(bx1,'Time (s)')
    xlim(bx1,[xlimMin,xlimMax])
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [mean(vfSum(:,keepIndex)) + SEM_vf(keepIndex)]; 
    SEMlow = [mean(vfSum(:,keepIndex)) - SEM_vf(keepIndex)];
    plot(bx2,jumpTime,mean(vfSum,'omitnan'),'b','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx2,[xlimMin,xlimMax])
    ylabel(bx2,'forward velocity (mm/s)')
    xlabel(bx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(vySum(:,keepIndex),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(vySum(:,keepIndex),'omitnan') - SEM_vy(keepIndex)];
    plot(bx3,jumpTime,mean(vySum,'omitnan'),'b','LineWidth',1.5)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx3,[xlimMin,xlimMax])
    ylabel(bx3,'rotational speed (deg/s)')
    xlabel(bx3,'Time (s)')

    vfSum = zeros(size(Sum90_not,1),size(Sum180.jumpIdx{1} ,1));
    vySum = zeros(size(Sum90_not,1),size(Sum180.jumpIdx{1} ,1));
    vsSum = zeros(size(Sum90_not,1),size(Sum180.jumpIdx{1} ,1));
    ZSum = zeros(size(Sum90_not,1),size(Sum180.jumpIdx{1} ,1));
    jump_change90not = [];

    for jump = 1:size(Sum90_not,1) 
        if ~isempty(Sum90_not.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(abs(Sum90_not.Vm{jump}),'loess',smoothWin);
            vfSum(jump,:) = smoothdata(Sum90_not.velFor{jump}','loess',smoothWin_v); 
            vySum(jump,:) = smoothdata(abs(Sum90_not.velYaw{jump})','loess',smoothWin_v);
            vsSum(jump,:) = smoothdata(Sum90_not.velSide{jump}','loess',smoothWin_v);
            jump_change90not(jump) = (Sum90_not.changeJump_post{jump} - Sum90_not.changeJump_pre{jump});
        end
    end
    
        startCount = jump;
    for jump = 1:size(Summ90_not,1)
        if ~isempty(Summ90_not.jumpIdx{jump})
            ZSum(jump,:) = smoothdata(abs(Summ90_not.Vm{jump}),'loess',smoothWin);
            vfSum(startCount,:) = smoothdata(Summ90_not.velFor{jump}','loess',smoothWin_v); 
            vySum(startCount,:) = smoothdata(abs(Summ90_not.velYaw{jump})','loess',smoothWin_v);
            vsSum(startCount,:) = smoothdata(Summ90_not.velSide{jump}','loess',smoothWin_v);
            jump_change90not(startCount) = (Summ90_not.changeJump_post{jump} - Summ90_not.changeJump_pre{jump});
            startCount = startCount + 1; 
        end
    end
    
    meanZSum = mean(ZSum,1);   
    meanPreJump =  mean(meanZSum(mid - 1000:mid - 1)); 
    
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') -  meanPreJump + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') -  meanPreJump - SEM_Z(keepIndex)];
    plot(bx1,jumpTime,mean(ZSum,'omitnan') -  meanPreJump,'k','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    ylabel(bx1,'|Z|')
    xlabel(bx1,'Time (s)')
    xlim(bx1,[xlimMin,xlimMax])
    
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [mean(vfSum(:,keepIndex)) + SEM_vf(keepIndex)]; 
    SEMlow = [mean(vfSum(:,keepIndex)) - SEM_vf(keepIndex)];
    plot(bx2,jumpTime,mean(vfSum,'omitnan'),'k','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx2,[xlimMin,xlimMax])
    ylabel(bx2,'forward velocity (mm/s)')
    xlabel(bx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(vySum(:,keepIndex),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(vySum(:,keepIndex),'omitnan') - SEM_vy(keepIndex)];
    plot(bx3,jumpTime,mean(vySum,'omitnan'),'k','LineWidth',1.5)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    xlim(bx3,[xlimMin,xlimMax])
    ylabel(bx3,'rotational speed (deg/s)')
    xlabel(bx3,'Time (s)')
    linkaxes([bx1,bx2,bx3],'x')
    
    jump_change_all = [jump_change90';jump_change90not'];
    g1 = repmat({'High'},length(jump_change90),1);
    g2 = repmat({'Low'},length(jump_change90not),1);

    base2 = ones(size(jump_change90not))*1.4;  
    base1 = ones(size(jump_change90))*0.6;
    
    figure();hold on;
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    scatter(base1,jump_change90)
    scatter(base2,jump_change90not)
    xlim([0 2])

    [p,stats] = vartestn(jump_change_all,[g1;g2],'TestType','BrownForsythe');
    
    
    % plot change in distance between cell's pref HD and the fly's HD
    % produced by each jump
    
    base2 = ones(size(jump_change90))*1.4;  
    base1 = ones(size(jump_change180))*0.6;
    figure();hold on;
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    scatter(base1,abs(jump_change180))
    scatter(base2,abs(jump_change90))
    plot([0.6,1.4],[mean(abs(jump_change180)),mean(abs(jump_change90))],'o-k','MarkerFaceColor','k')
    xlim([0 2])
    ylim([0,12])
    title('Corrected 180 vs 90 jumps')
    
    [~,p] = ttest2(abs(jump_change180), abs(jump_change90));
    text(1.6,10,['p = ',num2str(p)])
    
    base2 = ones(size(jump_change90not))*1.4;  
    base1 = ones(size(jump_change180not))*0.6;
    
    figure();hold on;
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    scatter(base1,abs(jump_change180not))
    scatter(base2,abs(jump_change90not))
    plot([0.6,1.4],[mean(abs(jump_change180not)),mean(abs(jump_change90not))],'o-k','MarkerFaceColor','k')
    ylim([0,12])
    xlim([0 2])
    title('Uncorrected 180 vs 90 jumps')
    
    [~,p] = ttest2(abs(jump_change180not), abs(jump_change90not));
    text(1.6,10,['p = ',num2str(p)])
    
    % plot individual membrane potential traces around cue jumps, smoothed
    % for visualization
    
    smoothWin = 1; 
    
    ZSum180 = []; 
    for jump = 1:size(Sum180,1)  
        if ~isempty(Sum180.jumpIdx{jump})
            ZSum180(jump,:) = smoothdata((Sum180.Vm{jump}),'loess',smoothWin);
        end
     end
     
     ZSum180not = []; 
     for jump = 1:size(Sum180_not,1) 
        if ~isempty(Sum180_not.jumpIdx{jump})
            ZSum180not(jump,:) = smoothdata((Sum180_not.Vm{jump}),'loess',smoothWin);
        end
     end
     
     figure();
     set(gcf,'renderer','painters','color','w')
     subplot(1,2,1)
     plot(jumpTime,ZSum180,'color',[0,0,1,0.2])
     xlim([-2,2])
     subplot(1,2,2)
     plot(jumpTime,ZSum180not,'color',[0,0,0,0.2])
     xlim([-2,2])
     sgtitle('180 jumps')
    
    for jump = 1:size(Sum90,1) 
        if ~isempty(Sum90.jumpIdx{jump})
            ZSum90(jump,:) = smoothdata((Sum90.Vm{jump}),'loess',smoothWin);
        end
    end
    
        startCount = jump;
    for jump = 1:size(Summ90,1)
        if ~isempty(Summ90.jumpIdx{jump})
            ZSum90(startCount,:) = smoothdata((Summ90.Vm{jump}),'loess',smoothWin);
            startCount = startCount + 1; 
        end
    end

    for jump = 1:size(Sum90_not,1) 
        if ~isempty(Sum90_not.jumpIdx{jump})
            ZSum90not(jump,:) = smoothdata((Sum90_not.Vm{jump}),'loess',smoothWin);
        end
    end
    
        startCount = jump;
    for jump = 1:size(Summ90_not,1)
        if ~isempty(Summ90_not.jumpIdx{jump})
            ZSum90not(startCount,:) = smoothdata((Summ90_not.Vm{jump}),'loess',smoothWin);
            startCount = startCount + 1; 
        end
    end
    
    figure();
    set(gcf,'renderer','painters','color','w')
    subplot(1,2,1)
    plot(jumpTime,ZSum90,'color',[0,0,1,0.2])
    xlim([-2,2])
    subplot(1,2,2)
    plot(jumpTime,ZSum90not,'color',[0,0,0,0.2])
    xlim([-2,2])
    sgtitle('90 jumps')
     
     
     
end