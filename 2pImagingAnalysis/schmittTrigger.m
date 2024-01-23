function triggerIdx = schmittTrigger(xdata,highT,lowT)
% schmittTrigger
% implements a shmitt trigger and pulls indices accordingly. once x exceeds
% high threshold, trigger will remain active until x falls below low
% threshold
%
% INPUT
% xdata - dataset to apply trigger against
% highT - high threshold to cross in order to be active
% lowT - low threshold to cross in order to be inactive
% initialize
    
    if size(xdata,2) > size(xdata,1)
        xdata = xdata';
    end
    [x,n] = size(xdata);
    triggerIdx = zeros(x,n);

    for e = 1:n
        highIdx = xdata>=highT;
        lowIdx = xdata>=lowT;

        triggerIdx(:) = highIdx; %initialize w/all high points
        trig=0; %set trigger off
        for i=1:x
            % turn trigger ON if above HIGH
            if highIdx(i)==1
                trig=1;
            end
            % if above LOW add point
            if lowIdx(i)==1
                if trig==1
                    triggerIdx(i)=1;
                end
                % else, turn trigger OFF if below LOW
            elseif isnan(xdata(i)) % if nan set to prev value --> grace periods around jumps & not moving
                if i == 1
                    triggerIdx(i)=0;
                else
                    triggerIdx(i) = triggerIdx(i - 1);
                end
            else
                trig=0;
            end
        end
    end
end