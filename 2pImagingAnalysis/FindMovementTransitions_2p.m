function [startIdx, stopIdx] = FindMovementTransitions_2p(rootDir, toPlot)
% finds indicies of transitions between stopped and moving states

% INPUTS
% rootDir: directory containing experiment folders, string
% toPlot: whether or not to show debugging plot, boolean

    % gather folder & fly information
    folders = get_folders(rootDir,1,1); 

    all_folders = string(folders.folder);

    for f = 1:size(all_folders,1)
        [start,finish] = regexp(all_folders(f),'_fly', 'ignorecase');
        dateIdx = regexp(all_folders(f),'\');
        dateIdx = [dateIdx(3)+1:dateIdx(4)-1]; 
        Date = char(all_folders(f));
        Date = Date(dateIdx);
        if isempty(finish)
            [start,finish] = regexp(all_folders(f),'_Fly');
        end
        fly_temp = char(all_folders(f));
        fly = fly_temp(start:finish + 1);
        flyID = strcat(Date,fly);
        flies(f) = string(flyID);
    end
        flies = flies';
        uniqueFlies = unique(flies); 
        flyCount = zeros(size(all_folders));
        for fly = 1:length(uniqueFlies)
            flyCount(flies == uniqueFlies(fly)) = fly;
        end

    transitionThreshold = round(1.5,1); 
    % how close together can a movement transition occur
    bufferTime = 0.5;
    next_stopLength = bufferTime * 60; 
    prev_startLength = bufferTime * 60; 

    % iterate through flies
    for t = 1:length(uniqueFlies)

        fly = uniqueFlies(t); 
        fly_folders = all_folders(flies == fly);
        % for each fly iterate through trials
        for f = 1:size(fly_folders,1)
                folder = fly_folders(f);

                if strcmp(folder(end),'.')
                    folder = folder(1:end-2); 
                end

                 processedData_dir = fullfile(folder,'processed_data');

                % Get data files
                expID = get_expID(folder);
                expList = {expID};
                [~, trialMd] = load_metadata(expList, folder);
                roiData = load_roi_data(expList, folder);

                try
                    numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
                catch
                    numTrials = 1; 
                end
                
                for nTrial = 1:numTrials
                    % load processed data
                    clear ZData ftT 
                    data_filelist = dir(processedData_dir);
                    for files = 1:length(data_filelist)
                        if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                            load(fullfile(processedData_dir,data_filelist(files).name));
                        end
                    end
                    load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))

                    % calculate cumulative speed & smooth
                    total_mov_mm = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs((ftT.velYaw{1}))*4.5;
                    total_mov_smooth = smoothdata(total_mov_mm,'gaussian',60); 
                    start_transitions = zeros(size(total_mov_mm)); 
                    prev_startIdx = 1; 

                    % iterate through data indices and identify movement
                    % transitions
                    for i = 2:length(total_mov_mm)-1
                        % set a buffer in which a new start transition
                        % cannot occur if one occured right before it
                        prev_startWindow = (i - prev_startLength):i-1;
                        if prev_startWindow(1) < 1
                            prev_startWindow = 1:prev_startWindow(end);
                        end
                        
                        % a start transition occurs at the index that the fly's velocity
                        % passeds the minimum threshold, but only if no other
                        % start transitions happend directly before it
                        if  (total_mov_mm(i + 1) > transitionThreshold && total_mov_mm(i - 1) < transitionThreshold) && start_transitions(i-1) ~= 1  && (all(total_mov_mm(prev_startWindow) < transitionThreshold)) 
                            start_transitions(i) = 1;
                            prev_startIdx = i; 
                        end
                    end

                    stop_transitions = zeros(size(total_mov_mm)); 
                    next_stopIdx = 2; 


                    for i = 2:length(total_mov_mm)-1
                        % set a buffer in which a new stop transition
                        % cannot occur if one occurs right after it
                        next_stopWindow = i+1:i+next_stopLength;
                        if next_stopWindow(end) > length(total_mov_mm)
                            next_stopWindow = next_stopWindow(1):length(total_mov_mm); 
                        end

                        if  (total_mov_mm(i + 1) < transitionThreshold && total_mov_mm(i - 1) > transitionThreshold) && stop_transitions(i-1) ~= 1 && (all(total_mov_mm(next_stopWindow) < transitionThreshold)) 
                            stop_transitions(i) = 1; 
                        end
                    end
                    
                    % debugging plot to compare transition points to
                    % movement data
                    if toPlot
                        figure();
                        plot(ftT.trialTime{1},total_mov_mm)
                        hold on
                        plot(ftT.trialTime{1}(logical(stop_transitions)),total_mov_mm(logical(stop_transitions)),'ko')
                        plot(ftT.trialTime{1}(logical(start_transitions)),total_mov_mm(logical(start_transitions)),'ro')
                    end

                    startIdx = find(start_transitions == 1); 
                    stopIdx = find(stop_transitions == 1);
                end
        end
    end
end
