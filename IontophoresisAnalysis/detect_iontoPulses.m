function [pulse_table] = detect_iontoPulses(bData, pulse_window, sampRate)
     % INPUTS
     
     % bData: trial behaviour data, table
     % pulse_window: length (s) of window around ionto pulses to collect, int
     % sampRate: sampling rate of data, int
     
     % OUTPUTS
     
     % pulse_table: table summarizing pulses
    
    
    stim = bData.stim; 

    count = 1; 
    prev_jump = 0; 
    transition = zeros(size(stim)); 

    % collect start & endpoints of iontophoresis stimuli 
    for idx = 1:length(stim)
        value = stim(idx);
        if count ~= 1 
            if value ~= prev_val && (count - prev_jump > 5) % last transition had to occur at least 5 datapoints earlier
                transition(count) = 1;
                prev_jump = count; 
            else
                transition(count) = 0; 
            end

            if ~(isnan(value))
                prev_val = value;
            else
                prev_val = prev_val;
            end
            
            count = count + 1; 
        else
            transition(count) = 0; 
            prev_val = value;
            count = count + 1;
        end
    end
 
    window_datapoints = round(pulse_window*sampRate); 
    pulse_idx = find(transition == 1);
    count = 1; 
    for idx = 1:length(pulse_idx)
        if rem(idx,2) ~= 0 
            pulse_StartStopidx(count,1) = pulse_idx(idx);
        else
            pulse_StartStopidx(count,2) = pulse_idx(idx); 
            count = count + 1; 
        end
    end
    
    pulse_StartStopidx(:,3) = pulse_StartStopidx(:,2) - pulse_StartStopidx(:,1); 

    count = 1; 
    pulse_array = zeros(size(pulse_StartStopidx,1),5); 
    for idx = 1:size(pulse_StartStopidx,1)
        pre_jump = pulse_StartStopidx(idx,1) - window_datapoints;
        post_jump = pulse_StartStopidx(idx,1) + window_datapoints;
        if post_jump > length(bData.stim)
            post_jump = length(bData.stim);
        end
        pulse_array(count,1) = pre_jump; 
        pulse_array(count,2) = pulse_StartStopidx(idx,1); 
        pulse_array(count,3) = pulse_StartStopidx(idx,2); 
        pulse_array(count,4) = post_jump; 
        pulse_array(count,5) = pulse_StartStopidx(idx,3); 
        count = count + 1;
    end
    
    headers = {'windowStart', 'pulseStart','pulseEnd','windowEnd','pulseLength'};
    pulse_table = array2table(pulse_array,'VariableNames',headers); 
    
end