function [triggerIdx, rho, Meno_chunks, not_Meno_chunks,bData_interp] = SegmentMenovsNotMeno_ephys(behaviourData, window, minVel,highThres,lowThres)
% Performs index by index assignment into menotaxing category based on
% calculated rho value assigned to that index. Rho is calculated over a
% window of set length centered on the index including only data points
% when the fly was moving. 

% INPUTS

    % behaviourData: trial's fictrac data
    % window: window to calculate idx rho over
    % minVel: minimum total velocity value, below = not moving
    % highThres: threshold when crossed with a + slope indicates start of meno bought
    % lowThres: threshold when crossed with a - slope indicates end of meno bought

%

% OUTPUTS

% triggerIdx: category of each index, boolean
% rho: rho values of each index, vector
% Meno_chunks: cell containing indices categorized as 'menotaxing', cell
% not_Meno_chunks: cell containing indices categorized as not 'menotaxing', cell
% bData_interp: downsampled behaviour data, used to assign correct category
% to 1000Hz sampled data in Group_Meno_notMeno_ephys.m

Meno_chunks = {}; 
not_Meno_chunks = {};

% ephys data had to be downsampled from 1000Hz to 60Hz due to speed & memory constraints
newSampRate = 60; 

vq = behaviourData.time(1):1/newSampRate:behaviourData.time(end);
fn = fieldnames(behaviourData);
if istable(behaviourData)
    fn = fn(1:end - 3); 
end
for field = 1:numel(fn)
    if isempty(regexp(fn{field},'disp'))
        bData_interp.(fn{field}) = interp1(behaviourData.time,behaviourData.(fn{field}),vq,'linear');
    end
end  

total_mov_mm = abs(bData_interp.vel_for) + abs(bData_interp.vel_side) + abs(deg2rad(bData_interp.vel_yaw))*4.5;        
no0vel_idx = find(total_mov_mm >= minVel); 


count = 1; 

% find jump idxes 
[jump_array, transitions, ~] = detect_jumps_ephys(bData_interp.frY, 5,5, 60);
jump_idx = [];
for jump = 1:size(jump_array,1)
    jump_idx = [jump_idx , jump_array(jump,2):jump_array(jump,3)];
end


% calculate & assign rho value to each datapoint
mean_headingVectors = [];

all_idx =  no0vel_idx;
if isempty(jump_idx)
    noJump_idx = all_idx;
else
    noJump_idx = all_idx(~ismember(all_idx, jump_idx));
end

keep_idx = intersect(noJump_idx, no0vel_idx);

% remove 0 velocity & jump idxes from data to evaluate
fn = fieldnames(bData_interp);
if istable(bData_interp)
    fn = fn(1:end - 3); 
end
for field = 1:numel(fn)
        bData_interp.(fn{field}) = bData_interp.(fn{field})(keep_idx);
end 

window = window * newSampRate; 
if window > length(keep_idx)
    window = length(keep_idx); 
end

% slide window over remaining data idx by idx to calculate rho
for i = 1 - window/2:1:length(bData_interp.angle) - window/2
    idx = i:i + window; 

    if idx(end) > length(bData_interp.angle) && idx(1) < 1 
        idx = 1:1:length(bData_interp.angle);
    elseif idx(end) > length(bData_interp.angle)
        idx = idx(1):1:length(bData_interp.angle);
    elseif idx(1) < 1 
        idx = 1:idx(end); 
    end

    angles_flyFor = bData_interp.angle(idx);
    if ~isempty(angles_flyFor)
        x = cosd(angles_flyFor); 
        y = sind(angles_flyFor); 
        idx_windows{count,1} = idx;
        mean_headingVectors(1,count)= sum(x)/length(x); 
        mean_headingVectors(2,count)= sum(y)/length(y);
        count = count + 1;
    else
        mean_headingVectors(1,count)= nan; 
        mean_headingVectors(2,count)= nan; 
        idx_windows{count,1} = idx;
        slope(count) = nan;
        count = count + 1;
    end
end

rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 

% assign idx 0 or 1 values based on assigned rho values & thresholds
triggerIdxRho = schmittTrigger(rho',highThres,lowThres);

triggerIdx = zeros(size(triggerIdxRho)); 
triggerIdx(triggerIdxRho == 1) = 1; 


% connect meno chunks if non meno is shorter than x seconds (optional)
pastMenoIdx = 1;
pastNotMenoIdx = 1;
for i = 1:length(triggerIdx)
iValue = triggerIdx(i);
if i > 1
    if iValue == 1
        timeFromMeno = i - pastMenoIdx; 
        if timeFromMeno < 0 * newSampRate % last transition from menotaxis occured less than x seconds before
            triggerIdx(pastMenoIdx:i) = 1; 
        end
        pastMenoIdx = i; 
    elseif iValue == 0
        timeFromNotMeno = i - pastNotMenoIdx;
        if timeFromNotMeno < 0.5 * newSampRate
            triggerIdx(pastNotMenoIdx:i) = 0; 
        end
        pastNotMenoIdx = i; 
    end
end
end

% Group meno & not meno data chunks in arrays

chunk = []; 
Mcount = 1; 
nMcount = 1; 
for idx = 1:length(triggerIdx)
    if isempty(chunk)
        if ~isnan(triggerIdx(idx))
            chunk(1) = idx;
            chunkValue = triggerIdx(idx);
        end
    else
        if (triggerIdx(idx) == chunkValue || isnan(triggerIdx(idx))) && idx ~= length(triggerIdx)
            chunk = [chunk, idx];
        else
            if chunkValue == 1
                Meno_chunks{Mcount} = chunk;
                Mcount = Mcount + 1; 
                chunk = [];
            else
                not_Meno_chunks{nMcount} = chunk;
                nMcount = nMcount + 1; 
                chunk = [];
            end
        end
    end
end

% Create 2D path trajectory meno = red not meno = black

notNan_idx = find(~isnan(bData_interp.angle) & ~isnan(bData_interp.vel_for) & ~isnan(bData_interp.vel_side));
yawAngPos = bData_interp.angle(notNan_idx);
fwdAngVel = bData_interp.vel_for(notNan_idx);
slideAngVel = bData_interp.vel_side(notNan_idx);
triggerIdx = triggerIdx(notNan_idx); 
transition = zeros(size(bData_interp.vel_side));

N = keep_idx';
V = find(transitions == 1);

A = repmat(N,[1 length(V)]);
[~,closestIndex] = min(abs(A-V'));

transition(closestIndex) = 1; 

transition = transition(notNan_idx); 

newSampRate = 60; 
% conversion factor between degrees and mm
circum = 9 * pi; % circumference of ball, in mm
mmPerDeg = circum / 360; % mm per degree of ball

% position incorporating heading - as if fly were walking on x-y plane,
%  x-y coordinates at each time point
% start with fly at (0,0) and facing 0 deg
zeroedYawAngPos = yawAngPos - yawAngPos(1); 

% movement in x (in degrees) at each time point
xChangePos = (fwdAngVel ./ newSampRate) .* sind(zeroedYawAngPos) + ...
    (slideAngVel ./ newSampRate) .* sind(zeroedYawAngPos + 90);  

% x position in mm (i.e. x-coordinate of fly's position at each time 
%  point), starts at 0
xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
minX = min(xPos);
maxX = max(xPos);

% movement in y (in degrees) at each time point
yChangePos = (fwdAngVel ./ newSampRate) .* cosd(zeroedYawAngPos) + ...
    (slideAngVel ./ newSampRate) .* cosd(zeroedYawAngPos + 90);

% y position in mm (i.e. y-coordinate of fly's position at each timepoint), starts at 0
yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
minY = min(yPos);
maxY = max(yPos);

% summary plots
time = bData_interp.time;
nMenoTime = time(triggerIdx == 0);
MenoTime = time(triggerIdx ==1); 
timeStart = time(1); 
timeEnd = max(time); 
nRho = rho(triggerIdx == 0); 

mAngle = -bData_interp.angle(triggerIdx ==1);
mRho = rho(triggerIdx == 1); 

figure(77);clf;
set(gcf,'color','w','renderer','painters')
h(1) =  subplot(2,1,1);
hold on
a = plot(behaviourData.time, -behaviourData.angle,'k');
b = plot(MenoTime(MenoTime > timeStart & MenoTime < timeEnd),mAngle(MenoTime > timeStart & MenoTime < timeEnd),'r');
try
    b.XData(abs(diff(b.XData)) > 20) = nan;
catch
end
ylabel('HD')
ylim([-180,180])
xlim([timeStart,timeEnd])
box off
h(2) = subplot(2,1,2);
hold on
c = plot(nMenoTime(nMenoTime > timeStart & nMenoTime < timeEnd),nRho(nMenoTime > timeStart & nMenoTime < timeEnd),'k');
try
    % remove lines connecting cue angle transitions from 0 to 180 & jumps
    c.XData(abs(diff(c.XData)) > 10) = nan;
catch
end
d = plot(MenoTime(MenoTime > timeStart & MenoTime < timeEnd),mRho(MenoTime > timeStart & MenoTime < timeEnd),'r');
try
    d.XData(abs(diff(d.XData)) > 10) = nan;
catch
end
ylabel('rho')
box off
xlim([timeStart,timeEnd])
linkaxes(h,'x')


% plot the 2D path trajectory
figure(88);clf;patch('XData',xPos(triggerIdx == 0),'YData',yPos(triggerIdx == 0),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
hold on
patch('XData',xPos(triggerIdx == 1),'YData',yPos(triggerIdx == 1),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
patch('XData',xPos(1),'YData',yPos(1),'EdgeColor','g','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',20);
patch('XData',xPos(transition == 1),'YData',yPos(transition == 1),'EdgeColor','b','FaceColor','none','LineStyle','none','Marker','o', 'MarkerSize',7);
set(gcf,'color','w');
xlabel('mm')
xlim([min(minX,minY),max(maxX, maxY)])
ylim([min(minX,minY),max(maxX, maxY)])

title(['window: ',num2str(window/newSampRate),' highThres: ', num2str(highThres),' lowThres: ', num2str(lowThres)], 'Interpreter', 'none')
end