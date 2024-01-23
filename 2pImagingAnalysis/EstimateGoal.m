function est_goal = EstimateGoal(ftT, i, window, jump_idx,sampRate)
% calculates goal HD associated with a specific data index

% INPUTS
% ftT: behaviour data, table
% i: index to calculate associated goal, int
% window: number data points around i to calculate goal from, float
% jump_idx: indices when cue jumps occured, float
% sampRate: sampling rate of data, float

    % cumulative speed
    speed = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1})*4.5;
    % mult by samp rate to get correct time window
    window = round(window * sampRate); 
    % data window centered on index of interest
    idx = i - round(window/2):i + round(window/2); 
    % make sure window doesn't extend past start or end of trial 
    if idx(end) > length(ftT.cueAngle{1})
        idx = idx(1):1:length(ftT.cueAngle{1});
    elseif idx(1) < 1 
        idx = 1:idx(end); 
    end
    % remove idx from window that were influenced by jumps
    if sum(ismember(idx, jump_idx))
        badIdx = ismember(idx, jump_idx);
        idx = idx(badIdx == 0);
    end
        angle_temp = ftT.cueAngle{1}(idx); 
        speed_temp = speed(idx); 
        % remove indicies when fly was stopped 
        angles_flyFor = angle_temp(speed_temp > 3); 
        if ~isempty(angles_flyFor) 
            % convert polar to cartesian coordinates (distance = 1)
            x = cosd(angles_flyFor); 
            y = sind(angles_flyFor); 
            % calculate relative distance traveled in each axis
            mean_headingVectors(1)= sum(x)/length(x); 
            mean_headingVectors(2)= sum(y)/length(y);
        else
            mean_headingVectors(1)= nan; 
            mean_headingVectors(2)= nan; 
        end
        
    % calculated average HD fly maintained around the index (goal)
    est_goal = atan2(mean_headingVectors(2),mean_headingVectors(1)); 
end
        
        