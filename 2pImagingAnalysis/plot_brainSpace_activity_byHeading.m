function plot_brainSpace_activity_byHeading(folder, goal,rho) 
% plotting function used in PFL2_BrainSpace_plots to visualize the average
% pattern of activity across the PFL2 population in different directional
% error bins 

% INPUTS

% folder: folder name, string
% goal: estimated whole trial goal, float
% rho: estimated whole trial rho, float

    % bin width (degrees)
    window = 90; 

    threshold=3;

    % load in trial data
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

    % convert cue position to HD and remove indicies when the fly wasns't moving
    total_mov_mm = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);
    no0vel_idx = find(total_mov_mm > threshold);
    angle = -ftT.cueAngle{1};
    angle = angle(no0vel_idx); 
    angle = wrapTo180(wrapTo180(angle)-wrapTo180((goal)));
    edges_angle = [-180 - window/2:window:180 + window/2]; 

    Z = []; 
    for roi = 1:size(ZData,1)
        Z(:,roi) = ZData.(3){roi}(no0vel_idx); 
    end

    % iterate through the ROIs & bin by HD
    roiActivity = []; 
    for roi = 1:size(Z,2)
        activity = Z(:,roi);
        if isnan(activity) | isnan(angle)
            activity = activity( ~isnan(activity) & ~isnan(angle)); 
            angle = angle( ~isnan(activity) & ~isnan(angle)); 
        end
        [zscore, centers_angle] = binData(activity, angle, edges_angle);
        roiActivity(roi,:) = zscore; 
    end
    
    roiActivity(:,1) = (roiActivity(:,1) + roiActivity(:,end))./2;
    roiActivity(:,end) = []; 

    % plot binned activity 
    ncol = size(roiActivity,2);
    color = (cbrewer2('Greens', ncol));    
    color = circshift(color,8,1);
    
    brainSpace = linspace(180,-180,size(Z,2));
    legendNames = []; 
    legendNames(:,1)= centers_angle(1:end-1);
    legendNames = num2str(legendNames);
    
    figure();clf;
    set(gcf,'color','w','renderer','painters')
    hold on
    for h = 1:size(roiActivity,2)
        plot(brainSpace,roiActivity(:,h),'color',color(h ,:),'LineWidth',1.5)
    end
    legend({legendNames})
    xlabel('brain space')
    ylabel('z-scored df/f')
    legend boxoff
    title(num2str(rho))
    xlim([-180,180])


end




