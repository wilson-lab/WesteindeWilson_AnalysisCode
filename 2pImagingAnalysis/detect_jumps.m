function [jump_array_down, jump_array, transition] = detect_jumps(ftT, pre_jump_window, post_jump_window, volRate,down)
% Finds each index at which a cue jump occured, returns that index along 
% with the start and end index of windows of a set length around each jump  

% INPUTS

% ftT: fictrac data table, table
% pre_jump_window: time (s) prior to jump to collect, float/int
% post_jump_window: time (s) following jump to collect, float/int
% volRate: volume rate of imaging data (opt, [] if not using), float
% down: collect jump indices from downsampled behaviour data, boolean

% OUTPUTS

% jump_array_down: table summarizing downsampling jump data
% jump_array: table summarizing 60Hz jump data
% transition: jump indicies

    jump_array_down = [];
    jump_array = []; 
    %% original sampling rate
    if ~down
        % find jump timepoints & calculate jump size
        sampRate = 60; 
        yChannel = ftT.PanelsY{1}; 
        yChannel_temp = round(yChannel - [1 2 3 4]);
        yChannel_clean = zeros(size(yChannel));

        for row = 1:size(yChannel_temp,1)
            yChannel_idx = find(yChannel_temp(row,:) == 0);
            yChannel_clean(row,1) = yChannel_idx;
        end

        yChannel = yChannel_clean; 

        count = 1; 
        prev_jump = 0; 
        transition = zeros(size(yChannel)); 
        jumpSize = []; 
        jcount = 1; 

        for idx = 1:length(yChannel)
            value = yChannel(idx);
            if isnan(value)
                print ('nan approaching')
            end

            if count ~= 1 
                if value ~= prev_val && ~(isnan(value)) && (count - prev_jump > 5) 
                    transition(count) = 1;
                    prev_jump = count; 
                    jumpStep = yChannel_clean(idx+1) - yChannel_clean(idx-1); 
                    if jumpStep == 1
                        jumpSize(jcount) = 90; 
                    elseif jumpStep == -1
                        jumpSize(jcount) = -90; 
                    else
                        jumpSize(jcount) = 180; 
                    end
                    jcount = jcount + 1; 

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

        pre_window_datapoints = round(pre_jump_window*sampRate); 
        post_window_datapoints = round(post_jump_window*sampRate); 
        jump_idx = find(transition == 1); 
        jump_array = zeros(length(jump_idx),4); 

        count = 1; 

        for idx = 1:length(jump_idx)
            pre_jump = jump_idx(idx) - pre_window_datapoints;
            post_jump = jump_idx(idx) + post_window_datapoints;
            jump_array(count,1) = pre_jump; 
            jump_array(count,2) = jump_idx(idx); 
            jump_array(count,3) = post_jump; 
            jump_array(count,4) = jumpSize(count);
            count = count + 1;
        end


    end
    
    %% downsampled jump idx
    if down
        % find jump timepoints & calculate size
        yChannel = ftT.PanelsY{1}; 
        yChannel_temp = yChannel - [1 2 3 4];

        for row = 1:size(yChannel_temp,1)
            yChannel_idx = knnsearch(yChannel_temp(row,:)',0);
            yChannel_clean(row,1) = yChannel_idx;
        end

        yChannel = yChannel_clean; 

        count = 1; 
        prev_jump = 0; 
        jcount = 1; 
        transition = zeros(size(yChannel)); 

        for idx = 1:length(yChannel)
            value = yChannel(idx);
            if isnan(value)
                print ('nan approaching')
            end

            if count ~= 1 
                if value ~= prev_val && ~(isnan(value)) && (count - prev_jump > 5) 
                    transition(count) = 1;
                    prev_jump = count; 
                    jumpStep = yChannel_clean(idx+1) - yChannel_clean(idx-1); 
                    if jumpStep == 1
                        jumpSize(jcount) = 90; 
                    elseif jumpStep == -1
                        jumpSize(jcount) = -90; 
                    else
                        jumpSize(jcount) = 180; 
                    end
                    jcount = jcount + 1; 
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

        pre_window_datapoints = round(pre_jump_window*volRate); 
        post_window_datapoints = round(post_jump_window*volRate);  
        jump_idx = find(transition == 1); 
        jump_array = zeros(length(jump_idx),3); 

        count = 1; 

        for idx = 1:length(jump_idx)
            pre_jump = jump_idx(idx) - pre_window_datapoints;
            post_jump = jump_idx(idx) + post_window_datapoints;
            jump_array(count,1) = pre_jump; 
            jump_array(count,2) = jump_idx(idx); 
            jump_array(count,3) = post_jump; 
            jump_array(count,4) = jumpSize(count);
            count = count + 1;
        end

        jump_array_down = jump_array;
    end

        %Debugging plot 
    %     figure()
    %     plot(yChannel)
    %     hold on
    %     plot(transition,'ro')

    
end