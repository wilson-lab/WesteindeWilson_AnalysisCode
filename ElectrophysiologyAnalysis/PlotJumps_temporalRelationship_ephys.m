function PlotJumps_temporalRelationship_ephys(rootDir)       
% Makes plots between neural activity & behaviour around cue jumps zoomed
% in on the response that follows directly after the cue jump to judge
% whether neural activity or the fly's velocity changes first. Selecting
% for jumps where fly's performed a fast correction such that the time of
% response is similar across individual examples. 
 
    % gather experiment folder
    folders = get_folders_ephys(rootDir);

    % setup summary table     
    headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','prefHead','preJump','postJump','jump_vm','cueDiff'};
    jump_summary = cell2table(cell(1,15),'VariableNames',headers); 
    count = 1; 

    % Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    for ff = 1:folderNum
      
      % Get folder information
      folder = folders(ff).folder;
      
       
        % Load in fictrac & ephys data
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end
        processedDir = fullfile(folder,'processedData');
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        load(fullfile(processedDir,'pro_trialData.mat'))
        bData = pro_behaviourData{1}; 
        tData = pro_trialData{1}; 
        
        % window (s) around jumps to evaluate rho
        jumpWin = 15;
        [jump_array_rho, ~, ~] = detect_jumps_ephys(bData.frY, 30,jumpWin, 1000);        
        % window (s) around jumps to average & plot
        [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, jumpWin,jumpWin, 1000);

        % z-score data better allow each cell to have equal influence on
        % slope of response regardless of average range of membrane
        % potential values
        activity = (tData.smooth_Vm - median(tData.smooth_Vm))/mad(tData.smooth_Vm);
        
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
            timeMov = (sum(abs(bData.vel_for(preJumpIdx_rho)) + abs(bData.vel_side(preJumpIdx_rho)) + abs(deg2rad(bData.vel_yaw(preJumpIdx_rho))*4.5) > 3))/1000; % seconds
            % average baseline Vm pre jump
            preJump_base = mean(activity(preJumpIdx));
            % baseline corrected jump Vm
            jump_vm = activity(jumpIdx) - preJump_base; 

            % calculate pre jump path straightness 
            [rho, ~] = CalculateAverageHeading_ephys(bData,1.5, preJumpIdx_rho);
            % calculate pre jump angle 
            [~, preJumpAngle] = CalculateAverageHeading_ephys(bData,0, preJumpIdx);
            
            % set up data arrays
            postJumpIdx = jump_array(jump,2)+1;
            sizePostJump = size(bData.angle(postJumpIdx:jump_array(jump,3)));
            preJump = ones(sizePostJump) * preJumpAngle;                    
            cueDiff = angdiff(preJump,deg2rad(bData.angle(postJumpIdx:jump_array(jump,3)))); 

            % decide if jump was 'corrected for' or not by how close to pre
            % jump position the post jump cue was returned to
            if jumpSize == 180
                correct = sum(abs(cueDiff) < deg2rad(45),'omitnan'); 
            else
                correct = sum(abs(cueDiff) < deg2rad(30),'omitnan');
            end 

            % collect relevent info 
            jump_summary.folder(count) = {folder}; 
            jump_summary.trial(count) = {t}; 
            jump_summary.jumpSize(count) = {jumpSize}; 
            jump_summary.jumpIdx(count) = {jumpIdx};
            jump_summary.velFor(count) = {jump_vf};
            jump_summary.velYaw(count) = {jump_vy};
            jump_summary.velSide(count) = {jump_vs}; 
            jump_summary.cueAngle(count) = {jump_angle};
            jump_summary.timeMov(count) = {timeMov}; 
            jump_summary.preJump(count) = {bData.angle(jump_array(jump,2)-1)}; 
            jump_summary.postJump(count) = {bData.angle(jump_array(jump,2)+1)}; 
            jump_summary.jump_vm(count) = {jump_vm};
            jump_summary.cueDiff(count) = {cueDiff}; 

            % classify jumps - only want ones where the fly moved a lot
            if  timeMov >= 15  && correct >= 10 && rho >= 0.7
                jump_summary.corrected(count) = {1};
            elseif correct >= 10 && rho < 0.7
                jump_summary.corrected(count) = {2}; 
            elseif correct <= 0 && timeMov >= 15
                jump_summary.corrected(count) = {0}; 
            else
                jump_summary.corrected(count) = {3}; 
            end

            count = count + 1; 

        end
    end
    
    % separate jumps into relevent categories
    Sum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 1,:); 
    nSum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 0,:); 
    Sum90 = jump_summary(cell2mat(jump_summary.jumpSize) == 90 & cell2mat(jump_summary.corrected) == 1,:);
    nSum90 = jump_summary(cell2mat(jump_summary.jumpSize) == 90 & cell2mat(jump_summary.corrected) == 0,:);
    Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 1,:);
    nSumm90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 0,:);
    %% visualization plots
    
    jumpTime = [-jumpWin+(1/1000):1/1000:jumpWin+(1/1000)];
    plotWin = 2; 
    vfSum = zeros(size(Sum180,1),size(jump_vf,1));
    vySum = zeros(size(Sum180,1),size(jump_vf,1));
    vsSum = zeros(size(Sum180,1),size(jump_vf,1));
    VmSum = zeros(size(Sum180,1),size(jump_vf,1));
    cueDiff_sum = zeros(size(Sum180,1),size(cueDiff,1));
    for jump = 1:size(Sum180,1)  
        VmSum(jump,:) = Sum180.jump_vm{jump};
        vfSum(jump,:) = Sum180.velFor{jump}'; 
        vySum(jump,:) = Sum180.velYaw{jump}';
        vsSum(jump,:) = Sum180.velSide{jump}';
        cueDiff_sum(jump,:) = Sum180.cueDiff{jump}';
    end
        
    SEM_Vm = std(VmSum,[],1,'omitnan') / sqrt(size(VmSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    figure; 
    sgtitle('180 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    yyaxis left
    hold on
    keepIndex = ~isnan(SEM_Vm);
    SEMhigh = [mean(abs(VmSum(:,keepIndex)),'omitnan') + SEM_Vm(keepIndex)]; 
    SEMlow = [mean(abs(VmSum(:,keepIndex)),'omitnan') - SEM_Vm(keepIndex)];
    patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,1],'FaceAlpha',0.1,'EdgeColor','none')
    plot(jumpTime,(mean(abs(VmSum),1)),'b','LineWidth',1.5)
    xline(0)
    ylim([min(SEMlow),max(SEMhigh)])
    ylabel('|Z(Vm)|')
    yyaxis right
    hold on
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(abs(vySum(:,keepIndex)),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(abs(vySum(:,keepIndex)),'omitnan') - SEM_vy(keepIndex)];
    patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[1,0,0],'FaceAlpha',0.1,'EdgeColor','none')
    plot(jumpTime,mean(abs(vySum),1),'r','LineWidth',1.5)
    ylim([min(SEMlow),max(SEMhigh)])
    ylabel('|Yaw|')
    xlim([-0.2,1])


    figure;
    set(gcf,'renderer','painters','color','w')
    subplot(2,1,1);
    sgtitle('180 jumps')
    plot(jumpTime,VmSum,'k')
    xlim([-plotWin,plotWin])
    box off
    subplot(2,1,2)
    plot(jumpTime,nVmSum,'r')
    xlim([-plotWin,plotWin])
    box off

        
    vfSum = zeros(size(Sum90,1),size(jump_vf,1));
    vySum = zeros(size(Sum90,1),size(jump_vf,1));
    vsSum = zeros(size(Sum90,1),size(jump_vf,1));
    VmSum = zeros(size(Sum90,1),size(jump_vf,1));
    cueDiff_sum = zeros(size(Sum90,1),size(cueDiff,1));
    
    for jump = 1:size(Sum90,1)  
        VmSum(jump,:) = Sum90.jump_vm{jump};
        vfSum(jump,:) = Sum90.velFor{jump}'; 
        vySum(jump,:) = Sum90.velYaw{jump}';
        vsSum(jump,:) = Sum90.velSide{jump}';
        cueDiff_sum(jump,:) = Sum90.cueDiff{jump}';
    end
    
    jumpNum = jump;
    for jump = 1:size(Summ90,1)
        VmSum(jump + jumpNum,:) = Summ90.jump_vm{jump};
        vfSum(jump + jumpNum,:) = Summ90.velFor{jump}'; 
        vySum(jump + jumpNum,:) = Summ90.velYaw{jump}';
        vsSum(jump + jumpNum,:) = Summ90.velSide{jump}';
        cueDiff_sum(jump + jumpNum,:) = Summ90.cueDiff{jump}';
    end
    
    SEM_Vm = std(VmSum,[],1,'omitnan') / sqrt(size(VmSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));

    figure; 
    sgtitle('90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    yyaxis left
    hold on
    keepIndex = ~isnan(SEM_Vm);
    SEMhigh = [mean(abs(VmSum(:,keepIndex)),'omitnan') + SEM_Vm(keepIndex)]; 
    SEMlow = [mean(abs(VmSum(:,keepIndex)),'omitnan') - SEM_Vm(keepIndex)];
    patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,1],'FaceAlpha',0.1,'EdgeColor','none')
    plot(jumpTime,(mean(abs(VmSum),1)),'b','LineWidth',1.5)
    xline(0)
    ylim([min(SEMlow),max(SEMhigh)])
    ylabel('|Z(Vm)|')
    yyaxis right
    hold on
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [mean(abs(vySum(:,keepIndex)),'omitnan') + SEM_vy(keepIndex)]; 
    SEMlow = [mean(abs(vySum(:,keepIndex)),'omitnan') - SEM_vy(keepIndex)];
    patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[1,0,0],'FaceAlpha',0.1,'EdgeColor','none')
    plot(jumpTime,mean(abs(vySum),1),'r','LineWidth',1.5)
    ylim([min(SEMlow),max(SEMhigh)])
    ylabel('|Yaw|')
    xlim([-0.2,1])




    figure;
    set(gcf,'renderer','painters','color','w')
    subplot(2,1,1);
    sgtitle('90 jumps')
    plot(jumpTime,VmSum,'k')
    xlim([-plotWin,plotWin])
    box off
    subplot(2,1,2)
    plot(jumpTime,nVmSum,'r')
    xlim([-plotWin,plotWin])
    box off

end