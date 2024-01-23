function [xPos, yPos] = plot2DTrajectory(rootDir, toPlot, nTrial, showJumps)
% plot the 2D path of the fly 

% INPUTS

% rootDir: experiment folder, string
% toPlot: plot path, boolean
% nTrial: experiment trial, 'all' or int
% showJumps: show jumps on path, boolean

% OUTPUTS

% xPos: x position of fly
% yPos: y position of fly

    folders = get_folders(rootDir, 1, 0);
    
    if isempty(folders)
        folders(1).folder = rootDir; 
    end
   
    %% Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    for ff = 1:folderNum
          try
          %% Get folder information
              folder = folders(ff).folder;


            %% Load in fictrac & ROI data
                if strcmp(folder(end),'.')
                    folder = folder(1:end-2); 
                end

                % Get data files
                expID = get_expID(folder);
                expList = {expID};

                % Load metadata 
                [~, trialMd] = load_metadata(expList, folder);

                % Load imaging data
                roiData = load_roi_data(expList, folder);

                % Load FicTrac data
                [~,ftData_DAQ,~] = load_ft_data(expList, folder, 1, 0);


                if all(nTrial == 'all')
                    nTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
                    numTrials = 1:nTrials;
                else 
                    numTrials = nTrial; 
                end
                %% path code

                for nTrial = numTrials
                    yawAngPos = ftData_DAQ.cueAngle{nTrial};
                    fwdAngVel = ftData_DAQ.velFor{nTrial};
                    slideAngVel = ftData_DAQ.velSide{nTrial};
                    
                    sampRate = 30; 
                    % conversion factor between degrees and mm
                    circum = 9 * pi; % circumference of ball, in mm
                    mmPerDeg = circum / 360; % mm per degree of ball

                    % position incorporating heading - as if fly were walking on x-y plane,x-y coordinates at each time point
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

                    % y position in mm (i.e. y-coordinate of fly's position at each time point), starts at 0
                    yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
                    minY = min(yPos);
                    maxY = max(yPos);

                    [~, ~, jumps] = detect_jumps(ftData_DAQ, 2, 2,[],0);
                    time = seconds(ftData_DAQ.trialTime{1});

                    if toPlot
                        figure();patch('xData',xPos,'yData',yPos,'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                        hold on
                        if showJumps
                            patch(xPos(jumps == 1),yPos(jumps == 1),time(jumps == 1),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',10);
                        end
                        set(gcf,'Renderer','painters')
                        set(gcf,'color','w');
                        xlabel('mm')
                        xlim([min(minX,minY),max(maxX, maxY)])
                        ylim([min(minX,minY),max(maxX, maxY)])

                    end
                end
          catch
              disp(['Folder: ',folder,' failed'])
          end
    end
end