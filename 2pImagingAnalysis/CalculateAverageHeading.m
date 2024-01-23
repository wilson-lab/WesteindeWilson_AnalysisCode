function [rho, theta] = CalculateAverageHeading(ftT,minVel, data_idx)
% calculates the average rho and estimated goal HD over the entire data
% segement 

% INPUTS
% ftT: behaviour data, table
% minVel: minimum velocty value that counts as moving, float
% data_idx: data indices to calculate over, vector or string

    % if true calculate rho & goal over the entire trial
    if strcmp(data_idx,'all')
        data_idx = 1:length(ftT.cueAngle{1});
    end

    angle_temp = ftT.cueAngle{1}(data_idx); 
    % euclidian speed
    speed = sqrt(ftT.velFor{1}(data_idx).^2 + ftT.velSide{1}(data_idx).^2);
    % remove indicies when fly was stopped 
    angles_flyFor = angle_temp(speed > minVel); 
    % make sure the fly moved at somepoint
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
    
    % calculate & return straightness metric ( 0 to 1, rho) & goal HD (theta)
    rho = sqrt(mean_headingVectors(1).^2 + mean_headingVectors(2).^2); 
    theta = atan2(mean_headingVectors(2),mean_headingVectors(1));
end