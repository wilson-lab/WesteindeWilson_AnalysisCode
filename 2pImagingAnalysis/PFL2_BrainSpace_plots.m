function PFL2_BrainSpace_plots(rootDir)
% For each PB and FB PFL2 trial identify the estimated whole-trial goal, 
% bin data points from each ROI by 4 directional error bins relative to 
% that goal. Create plots of brain space (ROI) vs activity for each of the
% error bins. --> On average what did activity across the population look
% like when the fly was x degrees away from her goal HD and what was the
% bump amplitude difference when the fly was facing her goal HD vs when she
% was 180 degrees away?
% Create a scatter plot of bump amplitude difference vs rho and calcualte
% the strength of the relationship (Pearson's R)

% INPUT

% rootDir: directory to look for experimental folders

    foldersStruct = get_folders(rootDir,1,0);
    folders = string({foldersStruct.folder})'; 
    folders = [folders(contains(folders,'PB'));folders(contains(folders,'FB'))]; 

    % setup summary table
    headers = {'folder','angDiff','ampDiff_bump','ampVar_bump','ampDiff_z','ampVar_z','rho','estGoal','estAntiGoal','theta'};
    pattern_summary = cell2table(cell(1,10),'VariableNames',headers); 
    cnt = 1; 

    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    countFail = 1;
    for ff = 1:folderNum
        try
          %% Get folder information
          folder = folders(ff);
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

            threshold=3;    

            expID = get_expID(folder);
            expList = {expID};

            nTrial = 1; 

            [~,ftT, ~] = load_ft_data(expList, folder, 1, 0);

            processedData_dir = fullfile(folder,'processed_data');

            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
            end
            
            load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))
            
            [~, ~, rho_all, theta_all] = CalculateAverageHeading_wholeSegment(ftT,30,threshold,'all', 60);
            %%
            total_mov_mm = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);
            
            if sum(total_mov_mm > threshold)/60 > 60 % moved for at least 60 seconds
                %% calculate estimated goal HD based on whole trial activity
                
                no0vel_idx = find(total_mov_mm > threshold);
                
                % convert cue position to HD
                angle = -ftT.cueAngle{1}(no0vel_idx);

                window = 5; 
                edges_angle = [-180 - window/2:window:180 + window/2]; 

                Z = []; 
                for roi = 1:size(ZData,1)
                    Z(:,roi) = ZData.(3){roi}(no0vel_idx); 
                end

                % iterate through ROIs
                roiActivity = []; 
                for roi = 1:size(Z,2)
                    activity = Z(:,roi);
                    if isnan(activity) | isnan(angle)
                        activity = activity( ~isnan(activity) & ~isnan(angle)); 
                        angle = angle( ~isnan(activity) & ~isnan(angle)); 
                    end
                    % bin data within each roi across HDs
                    [binnedActivity, centers_angle] = binData(activity, angle, edges_angle);
                    roiActivity(roi,:) = binnedActivity; 
                end
                
                % HDs of -180 & 180 represent the same direction, combine
                roiActivity(:,1) = (roiActivity(:,1) + roiActivity(:,end))./2;
                roiActivity(:,end) = []; 
                centers_angle(:,1) = abs(centers_angle(:,1));
                centers_angle(:,end) = [];

                % fit sinusoid to binned activity for each HD
                [bump_params, ~, ~] = fit_sinusoid(roiActivity,[0,2*pi], 0);      
                
                % for each HD find the difference between the amplitude at
                % that HD and the HD + 180 degrees
                for b = 1:length(centers_angle)/2
                    bumpAmp_diff(b) = diff([bump_params.amp(b), bump_params.amp(b+length(centers_angle)/2)]); 
                end
                
                maxDiffidx = find(abs(bumpAmp_diff) == max(abs(bumpAmp_diff)));
                
                % find the HD & HD + 180 degrees that has the maximum
                % difference in bump amplitude & set the estimated goal as
                % the one with the smallest bump amplitude
                if bumpAmp_diff(maxDiffidx) > 0 
                    estGoal = centers_angle(maxDiffidx);
                    estAntiGoal = centers_angle(maxDiffidx + length(centers_angle)/2);
                else
                    estAntiGoal = centers_angle(maxDiffidx);
                    estGoal = centers_angle(maxDiffidx + length(centers_angle)/2);
                end

                %% Bin activity in quadrants relative to estimated goal HD
                angle = -ftT.cueAngle{1}(no0vel_idx);
                % shift HD relative to estimated goal HD
                angle = wrapTo180(wrapTo180(angle)-wrapTo180(estGoal));
                window = 90; 
                edges_angle = [-180 - window/2:window:180 + window/2];
                
                % iterate through ROIs
                roiActivity = []; 
                for roi = 1:size(Z,2)
                    activity = Z(:,roi);
                    if sum(isnan(activity)) | sum(isnan(angle))
                        activity = activity( ~isnan(activity) & ~isnan(angle)); 
                        angle = angle( ~isnan(activity) & ~isnan(angle)); 
                    end
                    [binnedActivity, centers_angle] = binData(activity, angle, edges_angle);
                    roiActivity(roi,:) = binnedActivity; 
                end
                
                % HDs of -180 & 180 represent the same direction, combine
                roiActivity(:,1) = (roiActivity(:,1) + roiActivity(:,end))./2;
                roiActivity(:,end) = []; 
                centers_angle(:,1) = abs(centers_angle(:,1));
                centers_angle(:,end) = [];
                
                % fit sinusoid to binned activity for each HD quadrant
                [bump_params, ~, ~] = fit_sinusoid(roiActivity,[0,2*pi], 0);

                % calculate the the bump amplitude difference between the
                % max & min bump amplitudes
                ampDiff_bump = max(bump_params.amp(1)) - min(bump_params.amp(3));
                Z_ranges = max(roiActivity,[],1)-min(roiActivity,[],1);
                ampDiff_Z = Z_ranges(1) - Z_ranges(3); 
                ampVar = var(bump_params.amp); 


                % collect summary data for each trial 
                pattern_summary.folder(cnt) = {folder};
                pattern_summary.ampDiff_bump(cnt) = {ampDiff_bump};
                pattern_summary.ampDiff_Z(cnt) = {ampDiff_Z};
                pattern_summary.ampVar(cnt) = {ampVar};
                pattern_summary.rho(cnt) = {rho_all}; 
                pattern_summary.estGoal(cnt) = {estGoal};
                pattern_summary.estAntiGoal(cnt) = {estAntiGoal};
                pattern_summary.theta(cnt) = {-rad2deg(theta_all)}; % rho is in cue pos coord not HD, correct for this

                cnt = cnt + 1; 

                % Plot activity pattern across brain space, 0 = HD with highest amp sinusoid across PFL2 neurons
                disp(folder)
                plot_brainSpace_activity_byHeading(folder, estGoal,rho_all)
                title([num2str(rho_all),' ',folder])
            end
        catch
             disp([folder, ' failed'])
        end  
    end
         for f = 1:size(pattern_summary,1)
         end

        % plots showing the max bump diff per trial against each trial's
        % average rho value
        figure();
        set(gcf,'color','w','renderer','painters')
        hold on
        scatter(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_bump))
        [p,~] = polyfit(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_bump),1); 
        x = linspace(0,max(cell2mat(pattern_summary.rho)),1000);
        fit1 = polyval(p,x);
        plot(x,fit1)
        [R,P] = corrcoef(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_bump));
        annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
        xlabel('rho')
        ylabel('max - min bump amp')
        
        figure();
        set(gcf,'color','w','renderer','painters')
        hold on
        scatter(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_Z))
        [p,~] = polyfit(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_Z),1); 
        x = linspace(0,max(cell2mat(pattern_summary.rho)),1000);
        fit1 = polyval(p,x);
        plot(x,fit1)
        [R,P] = corrcoef(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_Z));
        annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
        xlabel('rho')
        ylabel('max - min Z')

end





