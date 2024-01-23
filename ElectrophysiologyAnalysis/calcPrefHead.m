function [prefHead, bump_params] = calcPrefHead(bData, tData)
% estimates HD tuning curve, by binning neural activity by HD

% INPUTs
% bData: behaviour data, table
% tData, neural data, table

    total_mov_mm = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5);
    no0vel_idx = find(total_mov_mm >0);
    angle = bData.angle(no0vel_idx); 
    activity = tData.smoothVm;
    activity = activity(no0vel_idx);

    edges_angle = [-180:15:180];
    [angleBins, centers_angle, ~] = binData(activity(~isnan(activity)), angle(~isnan(angle)), edges_angle);
    % smooth a little to reduce influence of noise on peak detection
    prefHead_smooth = smoothdata(angleBins,'gaussian',6);
    
    figure(22);clf;
    plot(centers_angle,angleBins)
    hold on
    plot(centers_angle,prefHead_smooth)
    prefHead_smooth(isnan(prefHead_smooth)) = []; 
    prefHead = centers_angle(prefHead_smooth == max(prefHead_smooth));
    try
        bump_params = fit_sinusoid(prefHead_smooth, [0,2*pi], 1);
    catch
        bump_params.rs = 0 ; 
    end

end