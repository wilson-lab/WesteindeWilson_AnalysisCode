function makeRepTrace(trialDir)

    % load in trial data
    load(fullfile(trialDir,'pro_trialData.mat'));
    tData = processed_trialData{1};
    load(fullfile(trialDir,'pro_behaviourData.mat'));
    bData = processed_behaviourData{1};

    % open whole trial fig & choose plot time window
    figureDir = fullfile(trialDir,'figures','whole_trial1.fig');
    openfig(figureDir); 

    timeStart = input('start time: '); 
    timeEnd = input('end time: ');

    % smooth velocity data for visualization
    time = bData.time(bData.time > timeStart & bData.time < timeEnd);
    velFor = smoothdata(bData.vel_for(bData.time > timeStart & bData.time < timeEnd),'gaussian',900);
    velYaw = smoothdata(bData.vel_yaw(bData.time > timeStart & bData.time < timeEnd),'gaussian',900);
    stim = bData.stim(bData.time > timeStart & bData.time < timeEnd);
    stimOn = stim > 0 ; 

    activity = tData.scaledOutput(bData.time > timeStart & bData.time < timeEnd);
    
    red = [0.75, 0, 0, 0.5];

    % plot trace
    figure();
    hold on
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    h1 = subplot(5,1,1);
    h2 = subplot(5,1,2:3);
    h3 = subplot(5,1,4:5);
    hold([h1,h2,h3],'on') 

    plot(h1,time, activity, 'k')
    activityStim(~stimOn) = min(activity)-5; 
    activityStim(stimOn) = max(activity); 
    plot(h1,time, activityStim, 'color',red,'LineWidth',1.5)
    plot(h1,time(~stimOn), activityStim(~stimOn), 'color',[1,1,1],'LineWidth',3)
    ylim(h1,[min(activity)-10,max(activity)])
    plot(h2,time, velFor, 'k')
    ylim(h2,[min(velFor)-2,max(velFor)])
    fStim(~stimOn) = min(velFor)-1; 
    fStim(stimOn) = max(velFor); 
    plot(h2,time, fStim, 'color',red,'LineWidth',1.5)
    plot(h2,time(~stimOn), fStim(~stimOn), 'color',[1,1,1],'LineWidth',3)
    plot(h3,time, velYaw, 'k')
    ylim(h3,[min(velYaw)-10,max(velYaw)])
    yStim(~stimOn) = min(velYaw)-5; 
    yStim(stimOn) = max(velYaw); 
    plot(h3,time, yStim, 'color',red,'LineWidth',1.5)
    plot(h3,time(~stimOn), yStim(~stimOn), 'color',[1,1,1],'LineWidth',3)

end