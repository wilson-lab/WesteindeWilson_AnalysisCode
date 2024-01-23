function HDtuning_binby_goalvsprefHead(rootDir, minVel) 
% bin menotaxis ephys data by the difference between the cell's preferred
% HD and the fly's current estimated goal HD. Calculate & plot the average tuning
% curve between neural activity & HD within each bin. 
    
    all_folders = struct2table(get_folders_ephys(rootDir));
    all_folders = all_folders.folder; 
    [~, flyCount, uniqueFlies] = get_flies_ephys(all_folders);
    MenoData = []; 
    %% calculate difference b/w preferred HD & goal HD    
    % go through all ephys data and calculate the difference between that cell's 
    % preferred HD and the fly's goal angle for each meno chunk

    for currentFly = 1:length(uniqueFlies) % goes through each fly
        flyTrials = folders(flyCount == currentFly)';
        for trial = 1:size(flyTrials,1) % goes through each trial for a fly
           folder = flyTrials{trial};
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end
            % load data
            processedDir = fullfile(folder,'processedData');
            clear bData tData MenoDataFly
            try
                load(fullfile(processedDir,'pro_behaviourData.mat'))
                bData = pro_behaviourData{1};
                load(fullfile(processedDir,'pro_trialData.mat'))
                tData = pro_trialData{1};
            catch
                load(fullfile(folder,'pro_behaviourData.mat'))
                bData = processed_behaviourData{1};
                load(fullfile(folder,'pro_trialData.mat'))
                tData = processed_trialData{1};
            end
            load(fullfile(folder,'MenoDataFly_trial_1.mat'));
            % find cell's preferred HD
            [prefHead, ~] = calcPrefHead(bData, tData);

            vf = bData.vel_for;            
            vs = bData.vel_side;            
            vy = bData.vel_yaw;             
            
            % gets rid of trials w/ no heading change --> indicates problem
            MenoDataTrial = MenoDataFly(~ismembertol(MenoDataFly.rho,1,10^-10),:); 

            if size(MenoDataTrial,1) > 0
                for r = 1:size(MenoDataTrial,1)
                    MenoDataTrial.fRate{r} = tData.fRate_sec(MenoDataTrial.Indices{r});
                end
            end
            
            % can get confused b/w 360 & 0
            if length(prefHead) > 1
                prefHead = prefHead(1); 
            end

            pHeads = ones(size(MenoDataTrial,1),1) * prefHead; 
            MenoDataTrial.prefHead = pHeads;
            % neg = goal cue pos is counterclockwise to prefhead cue pos
            MenoDataTrial.goalDiff = rad2deg(angdiff(-deg2rad(pHeads),-MenoDataTrial.Goal)); 

            % flip r PFL3 tunings to combine with l PFL3 data
            if (isempty(regexp(folder,'R')))
                MenoDataTrial.goalDiff = -MenoDataTrial.goalDiff;
            end
            MenoDataTrial.fly = ones(size(MenoDataTrial,1),1) * currentFly;
            if size(MenoDataTrial,1) > 0
                MenoData = [MenoData; MenoDataTrial];
            end        
        end
    end

    %% Assign data chunks to prefHead-goal diff bins

    step= 72; % degree difference between prefHead-goal bins
    goalDiffBins = [-180 + step/2:step:180 - step/2];
    step2 = 72; % degree difference between head direction bins
    binMeno = []; 
    binVMMeno = [];
    count = 1;

    % iterate through flies
    for currentFly = 1:length(uniqueFlies)
        disp(['Fly ',num2str(currentFly)])
        MenoDataFly = MenoData(MenoData.fly == currentFly,:);

        % Assign prefhead-goal diffs to correct bins
        for bin = 1:length(goalDiffBins)
            if bin == length(goalDiffBins)
                binMeno{currentFly,bin} = MenoDataFly(MenoDataFly.goalDiff >= goalDiffBins(bin) - step/2 & MenoDataFly.goalDiff <= goalDiffBins(bin) + step/2,:);
            else
                binMeno{currentFly,bin} = MenoDataFly(MenoDataFly.goalDiff >= goalDiffBins(bin) - step/2 & MenoDataFly.goalDiff < goalDiffBins(bin) + step/2,:);
            end
        end

     % iterate through each data chunk assigned to menotaxing bins 
        for bin = 1:size(binMeno,2) 
            if ~isempty(binMeno{currentFly,bin})
                binbData = [];
                bintData = [];
                goalDiffs = []; 
                for trial = 1:size(binMeno{currentFly,bin},1)   
                    folder = binMeno{currentFly,bin}.Folder(trial);
                    if strcmp(folder(end),'.')
                        folder = folder(1:end-2); 
                    end
                    clear bData tData 
                    processedDir = fullfile(folder,'processedData');
                    try
                        load(fullfile(processedDir,'pro_behaviourData.mat'))
                        bData = pro_behaviourData{1};
                        load(fullfile(processedDir,'pro_trialData.mat'))
                        tData = pro_trialData{1};
                        try
                            tData.smoothVm = (tData.smoothVm - vmMed(currentFly));
                        catch
                            tData.smooth_Vm = (tData.smooth_Vm - vmMed(currentFly));
                        end
                    catch
                        load(fullfile(folder,'pro_behaviourData.mat'))
                        bData = processed_behaviourData{1};
                        load(fullfile(folder,'pro_trialData.mat'))
                        tData = processed_trialData{1};
                        try
                           tData.smoothVm = (tData.smoothVm - vmMed(currentFly));
                        catch
                            tData.smooth_Vm = (tData.smooth_Vm - vmMed(currentFly));
                        end
                    end

                    tData.fRate_sec = (tData.fRate_sec - fRateMed(currentFly));

                    % set estimated preferred HD to 0, all HD relative to
                    % pref HD
                    bData.angle = wrapTo180(wrapTo180(bData.angle)-wrapTo180(binMeno{currentFly,bin}.prefHead(trial)));

                    binbData = [binbData; bData(binMeno{currentFly,bin}.Indices{trial},:)];
                    bintData = [bintData; tData(binMeno{currentFly,bin}.Indices{trial},:)];

                    if isempty(goalDiffs)
                        goalDiffs(1:length(binMeno{currentFly,bin}.Indices{trial})) = binMeno{currentFly,bin}.goalDiff(trial);
                    else
                        goalDiffs(end:end + length(binMeno{currentFly,bin}.Indices{trial})) = binMeno{currentFly,bin}.goalDiff(trial);
                    end                
                end

                vf = binbData.vel_for;            
                vs = binbData.vel_side;            
                vy = binbData.vel_yaw;             
                angle = binbData.angle; 
                total_mov_mm = abs(vf) + abs(vs) + abs(deg2rad(vy))*4.5;
                no0vel_idx = find(total_mov_mm >= minVel);
                angle = angle(no0vel_idx);
      
                % FR
                activity = bintData.fRate_sec;
                activity = activity(no0vel_idx);

                edges_angle = [-180:step2:180] ;       
                [angleBins, angleCenters, ~] = binData(activity, angle, edges_angle);

                % correct for binning edge error so tuning curve max value is in center
                if sum(isnan(angleBins)) <= 1
                    [~,maxIdx] = max(angleBins,[],'omitnan');
                    centerDiff =  3 - maxIdx;
                    binMeno{currentFly,bin} = circshift(angleBins, centerDiff); 
                    allTuningCurvesFR(count,:) = circshift(angleBins, centerDiff); 

                else 
                    binMeno{currentFly,bin} = [];
                    allTuningCurvesFR(count,1:length(angleCenters)) = 0;
                end

                % VM     
                activity = bintData.smoothVm;  


                activity = activity(no0vel_idx);

                edges_angle = [-180:step2:180] ;       
                [angleBins, angleCenters, ~] = binData(activity, angle, edges_angle);

                if sum(isnan(angleBins)) <= 1
                    [~,maxIdx] = max(angleBins,[],'omitnan');
                    centerDiff = 3 - maxIdx;
                    binVMMeno{currentFly,bin} = circshift(angleBins, centerDiff); 
                    allTuningCurvesVM(count,:) = circshift(angleBins, centerDiff); 

                else
                    binVMMeno{currentFly,bin} = [];
                    allTuningCurvesVM(count,1:length(angleCenters)) = 0;
                end  

                count = count + 1;             
            end
        end
    end

    %% Plot

    numCol = length(goalDiffBins);
    numRow = 2; 

    figure();
    set(gcf,'color','w','renderer','painters')
    clear ax
    for p = 1:numCol*2
        ax(p) = subplot(numRow,numCol,p);
        hold on
    end

    aveTrace = [];
    SEM_activity = []; 
    for e = 1:size(binMeno,2)
        Col = (binMeno(find(~cellfun(@isempty,binMeno(:,e))),e));
        numCells = size(Col,1);
        Col =cell2mat(Col);
        Col = reshape(Col,[360/step2,numCells]);
        % subtract minimum value to concentrate on amplitude difference
        Col = Col - min(Col,[],1);
        % calculate mean tuning curve within bin
        aveTrace(:,e) = mean(Col,2,'omitnan');
        SEM_activity(:,e) = std(Col,[],2,'omitnan')/sqrt(size(Col,2));
    end

    % uncomment to plot individual traces rather than SEM 
    % for c = 1:size(binMeno,1)
    %     for d = 1:size(binMeno,2)
    %         if ~isempty(binMeno{c,d}) %%&& isempty(regexp(uniqueFlies(c),'R'))
    %             plot(ax(d), angleCenters, binMeno{c,d} - min(binMeno{c,d},[],'omitnan'))
    %             
    %             %aveTrace = aveTrace + binMeno{c,d};
    %         end
    %         plot(ax(d), angleCenters, aveTrace(:,d),'k','LineWidth',1.5)
    %            ylim(ax(d),[0 45])
    % %     ylim(ax(d),[0 31])
    %     end
    % end



    for d = 1:size(binMeno,2)
        SEMhigh = aveTrace(:,d) + SEM_activity(:,d);
        SEMlow = aveTrace(:,d) - SEM_activity(:,d);
        patch(ax(d),[angleCenters fliplr(angleCenters)],[SEMhigh' fliplr(SEMlow')],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
        plot(ax(d), angleCenters, aveTrace(:,d),'k','LineWidth',1.5)
        title(ax(d),num2str(goalDiffBins(d)))
       % ylim(ax(d),[0 38])
        ylim(ax(d),[0 20])
    end

    for p = 1 + numCol:numCol * 2
        ax(p) = subplot(numRow,numCol,p);
        hold on
    end


    aveTrace = [];
    SEM_activity = []; 
    for e = 1:size(binVMMeno,2)
        Col = (binVMMeno(find(~cellfun(@isempty,binVMMeno(:,e))),e));
        numCells = size(Col,1);
        Col =cell2mat(Col);
        Col = reshape(Col,[360/step2,numCells]); 
        % subtract minimum value to concentrate on amplitude difference
        Col = Col - min(Col,[],1);
        % calculate mean tuning curve within bin
        aveTrace(:,e) = mean(Col,2,'omitnan');
        SEM_activity(:,e) = std(Col,[],2,'omitnan')/sqrt(size(Col,2));
    end

    % uncomment to plot individual traces rather than SEM 
    % for c = 1:size(binVMMeno,1)
    %     for d = 1:size(binVMMeno,2)    
    %         if ~isempty(binVMMeno{c,d}) %&& isempty(regexp(uniqueFlies(c),'R'))    
    %             plot(ax(d + numCol),angleCenters, binVMMeno{c,d} - min(binVMMeno{c,d},[],'omitnan'))     
    %         end
    %         plot(ax(d + numCol), angleCenters, aveTrace(:,d),'k','LineWidth',1.5)
    %         % %     ylim(ax(d + numCol),[-68,-52])
    %     ylim(ax(d + numCol),[0,18])
    %     end
    % end

    for d = 1:size(binVMMeno,2)
        SEMhigh = aveTrace(:,d) + SEM_activity(:,d); 
        SEMlow = aveTrace(:,d) - SEM_activity(:,d);
        patch(ax(d + numCol),[angleCenters fliplr(angleCenters)],[SEMhigh' fliplr(SEMlow')],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
        plot(ax(d + numCol), angleCenters, aveTrace(:,d),'k','LineWidth',1.5)
    %   ylim(ax(d + numCol),[-68,-52])
        ylim(ax(d + numCol),[0,10])
    end
end
