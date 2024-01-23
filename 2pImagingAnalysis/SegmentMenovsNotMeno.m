function [triggerIdx, rho, Meno_chunks, not_Meno_chunks] = SegmentMenovsNotMeno(ftT, window, minVel,highThres,lowThres, folder, trial,volRate,toPlot)
% Performs index by index assignment into menotaxing category based on
% calculated rho value assigned to that index. Rho is calculated over a
% window of set length centered on the index including only data points
% when the fly was moving. 

% INPUTS

% ftT: behaviour data, table
% window: time window to calculate rho over for each datapoint (s), float
% minVel: velocity threshold total velocity must be over to count as moving, float
% highThres: rho value to be above to start goal directed segment, float
% low Thres: rho value when below that ends goal directed segment, float
% folder: experiment folder, string
% trial: trial number, int
% volRate: imaging acquisition rate, float
% toPlot: whether to plot, boolean

% OUTPUTS

% triggerIdx: category of each index, boolean
% rho: rho values of each index, vector
% Meno_chunks: cell containing indices categorized as 'menotaxing', cell
% not_Meno_chunks: cell containing indices categorized as not 'menotaxing', cell


    triggerIdx = [];
    rho = [];
    Meno_chunks = {}; 
    not_Meno_chunks = {};

    try

            sampRate = 60;     
            window = window * sampRate; 
            if window > length(ftT.velFor{1})
                window = length(ftT.velFor{1}); 
            end
            count = 1; 


            mean_headingVectors = [];
            speed = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);

            % gather indices around jumps to discount from rho calculation
            [~, jump_array, ~] = detect_jumps(ftT, 5, 5,volRate,0);
            jump_idx = [];
            for jump = 1:size(jump_array,1)
                jump_idx = [jump_idx , jump_array(jump,2):jump_array(jump,3)];
            end

            % iterate through every datapoint & calculate its associated rho
            for i = 1 - window/2:1:length(ftT.cueAngle{1}) - window/2
                idx = i:i + window; 
                if idx(end) > length(ftT.cueAngle{1})
                    idx = idx(1):1:length(ftT.cueAngle{1});
                elseif idx(1) < 1 
                    idx = 1:idx(end); 
                end

                if sum(ismember(idx, jump_idx)) % remove idx from window that were influenced by jumps
                    badIdx = ismember(idx, jump_idx);
                    idx = idx(badIdx == 0);
                end

                    angle_temp = ftT.cueAngle{1}(idx); 
                    speed_temp = speed(idx); 
                    % remove indicies when fly was stopped 
                    angles_flyFor = angle_temp(speed_temp > minVel); 
                    if ~isempty(angles_flyFor) 
                        % convert polar to cartesian coordinates (distance = 1)
                        x = cosd(angles_flyFor); 
                        y = sind(angles_flyFor); 
                        idx_windows{count,1} = idx;
                        % calculate relative distance traveled in each axis
                        mean_headingVectors(1,count)= sum(x)/length(x); 
                        mean_headingVectors(2,count)= sum(y)/length(y);
                        count = count + 1;
                    else
                        mean_headingVectors(1,count)= nan; 
                        mean_headingVectors(2,count)= nan; 
                        slope(count) = nan;
                        count = count + 1;
                    end
            end
            % calculate straightness metric ( 0 to 1, rho) for each index
            rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 

            triggerIdx = schmittTrigger(rho',highThres,lowThres);

            % connect chunks separated by transient switches into the other category
            pastMenoIdx = 1;
            pastNotMenoIdx = 1;
            for i = 1:length(triggerIdx)
                iValue = triggerIdx(i);
                if i > 1
                    if iValue == 1
                        timeFromMeno = i - pastMenoIdx; 
                        if timeFromMeno < 0.5 * sampRate % last transition from menotaxis occured less than 0.5 seconds before
                            triggerIdx(pastMenoIdx:i) = 1; 
                        end
                        pastMenoIdx = i; 
                    elseif iValue == 0
                        timeFromNotMeno = i - pastNotMenoIdx;
                        if timeFromNotMeno < 0.5 * sampRate
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
                        if triggerIdx(idx) == chunkValue || isnan(triggerIdx(idx))
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

    %% plot 2D path trajectory meno = red not meno = black
    if toPlot

                        yawAngPos = ftT.cueAngle{1};
                        fwdAngVel = ftT.velFor{1};
                        slideAngVel = ftT.velSide{1};

                        sampRate = 60;
                        % conversion factor between degrees and mm
                        circum = 9 * pi; % circumference of ball, in mm
                        mmPerDeg = circum / 360; % mm per degree of ball

                        % position incorporating heading - as if fly were walking on x-y plane,
                        %  x-y coordinates at each time point
                        % start with fly at (0,0) and facing 0 deg
                        zeroedYawAngPos = yawAngPos - yawAngPos(1); 

                        % movement in x (in degrees) at each time point
                        xChangePos = (fwdAngVel ./ sampRate) .* sind(zeroedYawAngPos) + ...
                            (slideAngVel ./ sampRate) .* sind(zeroedYawAngPos + 90);  

                        % x position in mm (i.e. x-coordinate of fly's position at each time 
                        %  point), starts at 0
                        xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
                        minX = min(xPos);
                        maxX = max(xPos);

                        % movement in y (in degrees) at each time point
                        yChangePos = (fwdAngVel ./ sampRate) .* cosd(zeroedYawAngPos) + ...
                            (slideAngVel ./ sampRate) .* cosd(zeroedYawAngPos + 90);

                        % y position in mm (i.e. y-coordinate of fly's position at each time 
                        %  point), starts at 0
                        yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
                        minY = min(yPos);
                        maxY = max(yPos);

                        time = seconds(ftT.trialTime{1});
                        nMenoTime = time(triggerIdx == 0);
                        MenoTime = time(triggerIdx ==1); 
                        timeStart = time(1); 
                        timeEnd = max(time); 
                        nxPos = xPos(triggerIdx == 0);
                        nyPos = yPos(triggerIdx == 0);
                        angle =  -ftT.cueAngle{1}; 
                        nAngle = angle(triggerIdx ==0); 
                        nRho = rho(triggerIdx == 0); 

                        mAngle = angle(triggerIdx ==1);
                        mxPos = xPos(triggerIdx == 1);
                        myPos = yPos(triggerIdx == 1);
                        mRho = rho(triggerIdx == 1); 

                        clf
                        figure(66);
                        patch('XData',nxPos(nMenoTime > timeStart & nMenoTime < timeEnd),'YData',nyPos(nMenoTime > timeStart & nMenoTime < timeEnd),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                        hold on
                        patch('XData',mxPos(MenoTime > timeStart &  MenoTime < timeEnd),'YData',myPos(MenoTime > timeStart &  MenoTime < timeEnd),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                        set(gcf,'color','w');
                        xlabel('mm')
                        xlim([min(minX,minY),max(maxX, maxY)])
                        ylim([min(minX,minY),max(maxX, maxY)])

                        saveas(gcf, fullfile(folder,'images',['menoVsnonMeno_Path_00',num2str(trial),'.fig']))

    end
    catch
        disp(['folder ',folder,'trial_',num2str(trial),' failed'])
    end
end