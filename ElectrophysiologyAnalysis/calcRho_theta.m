function [rho, theta] = calcRho_theta(bData, window, minVel, step)
% calculates rho and theta values associated with each data index. Due to
% processing time steps a window across the data and bulk assigns the same
% value to a set number of datapoints rather than a true sliding window.

% Inputs
% bData: behaviour data, table
% window: window length (seconds) to calculate rho & theta over, float
% minVel: movement threshold, float
% step: number of data points to step by, float

    % convert window to num data points
    window = window * 1000; 
    if window > length(bData.vel_for)
        window = length(bData.vel_for); 
    end
    count = 1;
    % find when cue jumps happened
    jump_idx = find(bData.jumps == 1);
   
    % step through data indicies 
    for i = step/2:step:length(bData.angle) - step/2
        
        % grab data window centered on the data index
        idx = i - window/2:i + window/2; 

        % make sure window does extend past the start or end of trial 
        if idx(end) > length(bData.angle)-step/2 && idx(1) < step/2 
            idx = step/2:1:length(bData.angle)-step/2;
        elseif idx(end) > length(bData.angle)-step/2
            idx = idx(1):1:length(bData.angle)-step/2;
        elseif idx(1) < step/2 
            idx = step/2:idx(end); 
        end
        
        % remove idx from window that were influenced by jumps
        if sum(ismember(idx, jump_idx))
            badIdx = ismember(idx, jump_idx);
            idx = idx(badIdx == 0);
        end

        angle_temp = bData.angle(idx); 
        % euclidian speed
        speed = sqrt(bData.vel_for.^2 + bData.vel_side.^2);
        speed_temp = speed(idx); 
        % remove indicies where fly didn't move
        angles_flyFor = angle_temp(speed_temp > minVel); 
        angles_flyFor = bData.angle(idx);
        if ~isempty(angles_flyFor)
            % convert from polar to cartesian coordinates
            x = cosd(angles_flyFor); 
            y = sind(angles_flyFor);
            idx_windows{count,1} = idx;
            % calculate relative distance travelled along both axes
            mean_headingVectors(1)= sum(x)/length(x); 
            mean_headingVectors(2)= sum(y)/length(y);
            count = count + 1;
        else
            mean_headingVectors(1)= nan; 
            mean_headingVectors(2)= nan; 
            idx_windows{count,1} = idx;
            slope(count) = nan;
            count = count + 1;
        end
            
        rho(i - step/2 + 1:i+step/2) = sqrt(mean_headingVectors(1)^2 + mean_headingVectors(2)^2); 
        theta(i - step/2 + 1:i+step/2) = atan2(mean_headingVectors(2),mean_headingVectors(1));  
    end
end