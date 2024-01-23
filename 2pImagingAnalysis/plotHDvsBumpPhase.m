function plotHDvsBumpPhase(folder)
% calculates the change in HD and change in PFL2 bump phase at each
% timepoint then creates a scatter plot of these values and calculates the
% relationship between them 

% INPUT
% folder, trial folder to examine, string
% load data
load(fullfile(folder,'processed_data','bump_parameters_Trial001.mat'))
load(fullfile(folder,'processed_data','fictracData_Trial_001.mat'))

% cumulative velocity 
total_mov_mm = abs(ftT.velFor{1} + abs(ftT.velSide{1}) + abs(ftT.velYaw{1})*4.5);
no0vel_idx = find(total_mov_mm > 3);
% convert to agree with ROI & head direction units
bPhase = -wrapToPi(bump_params.pos_rad(no0vel_idx) + pi); 
rValue = bump_params.adj_rs(no0vel_idx);
% only look at indicies where the fit quality is positive
bPhase = bPhase(rValue > 0.1);
angle = ftT.cueAngle{1}(no0vel_idx);
angle = angle(rValue > 0.1);

changebPhase = [];
changeHD = []; 

% for different lags between the imaging & behaviour data calculate the
% bump phase and HD changes for each datapoint 
for lags = 0:0.2:1
    % s lag bump phase relative to HD
    lag = round(60 * lags); 
    % 5 sec windows
    window = 1.5*60; 
    count = 1;
    % move a window across through data to calculate phase and HD changes
    % from start of time window to end
    for i = 1+lag:window:length(bPhase)
        idx = i:i + window; 
        if (idx(end) <= length(bPhase)) && idx(1) >= 1 
            % collected lagged HD data, convert to radians 
            heading = -deg2rad(wrapTo180(angle(idx - lag))); 
            bPhaseStart = bPhase(idx(1));
            bPhaseEnd = bPhase(idx(end));
            changebPhase(count) = angdiff(bPhaseStart,bPhaseEnd); % 2nd - 1st + = clockwise change, - = counterclockwise in phase
            HDstart = heading(1); 
            HDend = heading(end); 
            changeHD(count) = angdiff(HDstart,HDend); % 2nd - 1st + = clockwise change, - = counterclockwise in HD
            count = count + 1; 
        end
    end
    
    % remove nans
    nanIdx = isnan(changeHD) | isnan(changebPhase);
    changebPhase(nanIdx) = []; 
    changeHD(nanIdx) = []; 
    
    figure()
    set(gcf,'color','w','renderer','painters')
    hold on
    scatter(rad2deg(changeHD),rad2deg(changebPhase))
    [p,~] = polyfit(changeHD,changebPhase,1); 
    x = linspace(-180,180,1000);
    fit1 = polyval(p,x);
    plot(x,fit1)
    ylim([-180,180])
    xlim([-180,180])
    [R,P] = corrcoef(changeHD,changebPhase);
    annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
    title(['lag: ',num2str(lags)])
    xlabel(folder)
end