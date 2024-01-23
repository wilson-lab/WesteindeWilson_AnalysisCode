function activity_behaviour_correlation(rootDir)       
% calculate the correlation coefficient between neural activity and fly's
% movement at different lags. Look at jumps where fly responded fairly
% quickly after the cue jump 

    % get folders
     folders = get_folders_ephys(rootDir);
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','preJump','postJump','jump_vm','cueDiff'};
     jump_summary = cell2table(cell(1,14),'VariableNames',headers); 
     count = 1; 

    % Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);

    for ff = 1:folderNum
      
        % Get folder information
        folder = folders(ff).folder;
      
       
        % Load in data

        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        processedDir = fullfile(folder,'processedData');
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        load(fullfile(processedDir,'pro_trialData.mat'))
        bData = pro_behaviourData{1}; 
        tData = pro_trialData{1}; 

        % time (s) around jumps to look at
        jumpWin = 4;
        % calculate rho over a longer time window
        [jump_array_rho, ~, ~] = detect_jumps_ephys(bData.frY, 30,15, 1000);
        [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, jumpWin,jumpWin, 1000);

        activity = tData.smooth_Vm;
        
        % remove jumps where window extended beyond start & end of trial
        jumpIdx = [jump_array(1,1):jump_array(1,3)]';  
        jumpIdx_rho = [jump_array_rho(1,1):jump_array_rho(1,3)]';
        if jumpIdx(1) < 1 || jumpIdx_rho(1) < 1
            jump_array(1,:) = []; 
            jump_array_rho(1,:) = []; 
        end
        
        for jump = 1:size(jump_array,1)
            jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
            jumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,3)]';
            if jumpIdx(end) > size(bData,1) || jumpIdx_rho(end) > size(bData,1)
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
            [rho, ~] = CalculateAverageHeading_ephys(bData,1.5, preJumpIdx_rho);

            timeMov = (sum(abs(bData.vel_for(preJumpIdx_rho)) + abs(bData.vel_side(preJumpIdx_rho)) + abs(deg2rad(bData.vel_yaw(preJumpIdx_rho))*4.5) > 3))/1000; % seconds
            preJump_base = mean(activity(preJumpIdx));
            % subtract ave baseline, looking at change in membrane
            % potential
            jump_vm = activity(jumpIdx) - preJump_base;

            [~, preJumpAngle] = CalculateAverageHeading_ephys(bData,0, preJumpIdx);
            postJumpIdx = jump_array(jump,2)+ 100;
            sizePostJump = size(bData.angle(postJumpIdx:jump_array(jump,3)));
            preJump = ones(sizePostJump) * preJumpAngle;

            cueDiff = angdiff(preJump,deg2rad(bData.angle(postJumpIdx:jump_array(jump,3)))); 

            % determine if jump was corrected for or not
            if jumpSize == 180
                correct = sum(abs(cueDiff) < deg2rad(75),'omitnan'); 
            else
                correct = sum(abs(cueDiff) < deg2rad(40),'omitnan');
            end 

            % collect relevent info in summary table
            jump_summary.folder(count) = {folder}; 
            jump_summary.trial(count) = {t}; 
            jump_summary.jumpSize(count) = {jumpSize}; 
            jump_summary.rho(count) = {rho}; 
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
            jump_summary.correct(count) = {correct};

            % categorize jumps
            if  correct >= 10 && timeMov > 1 && rho > 0.88
                jump_summary.corrected(count) = {1};
            elseif correct < 10 && timeMov > 1
                jump_summary.corrected(count) = {0}; 
            else
                jump_summary.corrected(count) = {3};
            end

            count = count + 1; 

        end

    end 
    
    % split data into summary tables
    Sum90Corr = jump_summary(cell2mat(jump_summary.corrected) == 1 & abs(cell2mat(jump_summary.jumpSize)) == 90,:); 
    Sum180Corr = jump_summary(cell2mat(jump_summary.corrected) == 1 & cell2mat(jump_summary.jumpSize) == 180,:); 

    %% visualization plots 

    % iterate through jumps @ diff lags & calculate correlation 
    sampRate = 1000;
    lags = [-1:0.01:1];

    % 90 jumps
    corrWin = 1;
    corrValuesP = zeros(size(Sum90Corr,1),length(lags),2); 
    jumpI = jumpWin*sampRate + 1;

    for t = 1:length(lags)
        lag = round(lags(t)*1000); 
        for jump = 1:size(Sum90Corr,1)
            % + lag, comparing yaw values w/ activity values delayed by lag
            % - lag, comparing yaw values w/ activity values that preceded by lag
            activity = abs(Sum90Corr.jump_vm{jump}(jumpI + lag : jumpI + corrWin*sampRate + lag));          
            yaw = abs(Sum90Corr.velYaw{jump}(jumpI : jumpI + corrWin*sampRate));

            [R,p] = corr(activity, yaw,'Type','Pearson'); 

            corrValuesP(jump,t,1) = R; 
            corrValuesP(jump,t,2) = p;

        end
    end

    figure();
    meanCorr = mean(corrValuesP(:,:,1),1,'omitnan');

    SEM_corr = std(corrValuesP(:,:,1),[],1,'omitnan') / sqrt(size(corrValuesP(:,:,1),1));
    SEMhigh = [meanCorr + SEM_corr]; 
    SEMlow = [meanCorr - SEM_corr];
    set(gcf,'renderer','painters','color','w')
    patch([lags fliplr(lags)],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    hold on
    plot(lags,meanCorr,'k','LineWidth',1.5)
    xlabel('Lag (s)')
    ylabel('R')
    title('90 jumps Pearson')
    box off

    % 180 jumps
    corrWin = 1;
    corrValuesP = zeros(size(Sum180Corr,1),length(lags),2); 
    jumpI = jumpWin*sampRate + 1;

    for t = 1:length(lags)
        lag = round(lags(t)*1000); 
        for jump = 1:size(Sum180Corr,1)
            % + lag, comparing yaw values w/ activity values delayed by lag
            % - lag, comparing yaw values w/ activity values that preceded by lag
            activity = abs(Sum180Corr.jump_vm{jump}(jumpI + lag : jumpI + corrWin*sampRate + lag));          
            yaw = abs(Sum180Corr.velYaw{jump}(jumpI : jumpI + corrWin*sampRate));

            [R,p] = corr(activity, yaw,'Type','Pearson'); 

            corrValuesP(jump,t,1) = R; 
            corrValuesP(jump,t,2) = p;

        end
    end

    figure();
    meanCorr = mean(corrValuesP(:,:,1),1,'omitnan');

    SEM_corr = std(corrValuesP(:,:,1),[],1,'omitnan') / sqrt(size(corrValuesP(:,:,1),1));
    SEMhigh = [meanCorr + SEM_corr]; 
    SEMlow = [meanCorr - SEM_corr];
    set(gcf,'renderer','painters','color','w')
    patch([lags fliplr(lags)],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    hold on
    plot(lags,meanCorr,'k','LineWidth',1.5)
    xlabel('Lag (s)')
    ylabel('R')
    title('90 jumps Pearson')
    box off
    
end

    