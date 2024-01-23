function [rho, theta] = CalculateAverageHeading_ephys(behaviourData,minVel, data_idx)
% calculates average rho and theta values across all indicies 

% Inputs
% minVel: movement threshold, float
% data_idx: indices to calculate rho & theta over, vector, string

    if strcmp(data_idx,'all')
        data_idx = 1:length(behaviourData.angle);
    end

    angle_temp = behaviourData.angle(data_idx); 
    % cumulative velocity
    speed = abs(behaviourData.vel_for(data_idx)) + abs(behaviourData.vel_side(data_idx)) + abs(deg2rad(behaviourData.vel_yaw(data_idx)))*4.5;    
    % remove indicies where fly didn't move
    angles_flyFor = angle_temp(speed > minVel); 
    if ~isempty(angles_flyFor) 
        % convert from polar to cartesian coordinates
        x = cosd(angles_flyFor); 
        y = sind(angles_flyFor);
        % calculate relative distance traveled along each axis
        mean_headingVectors(1)= sum(x)/length(x); 
        mean_headingVectors(2)= sum(y)/length(y);
    else
        mean_headingVectors(1)= nan; 
        mean_headingVectors(2)= nan; 
    end
    
    rho = sqrt(mean_headingVectors(1).^2 + mean_headingVectors(2).^2); 
    theta = atan2(mean_headingVectors(2),mean_headingVectors(1));
end