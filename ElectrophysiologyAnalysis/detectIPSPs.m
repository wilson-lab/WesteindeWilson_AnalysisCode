function jumpSum = detectIPSPs(expFolder)       
% First find all jumps where the fly didn't move before & after then detect
% IPSP frequency & amplitude before and after the cue jumps. Collect across
% all flies & jumps and plot the change in the IPSP parameters that
% happened before & after the jump against how much the jump changed the
% fly's HD from the cell's preferred HD and calculate the direction &
% strength of the relationship. 


    % creates a table summarizing which trials had jumps where the fly didn't
    % move before or after & which jump numbers those were.
    noMov_jumps = Find_0velJumps(expFolder);


    folders = noMov_jumps.folder;
    for f = 1:length(folders)
        folderTemp = char(folders(f));
        fliesAll(f) = string(folderTemp(1:end -7));
    end
    fliesAll = fliesAll';
    flies = unique(fliesAll); 

    jCount = 1;
    % iterate through each fly
    for fly = 1:size(flies,1)
        flyJumps = noMov_jumps(fliesAll == flies(fly),:);
        trials = unique(flyJumps.folder); 
        % iterate through each trial
        for t = 1:size(trials,1)   
            % collect the non moving jumps
            trialJumps = flyJumps(flyJumps.folder == trials(t),:);
            if size(trialJumps,1) > 1        
                folder = trialJumps.folder(1); 
                if strcmp(folder(end),'.')
                    folder = folder(1:end-2); 
                end

                % Load in data
                processedDir = fullfile(folder,'processedData');
                clear bData tData jump_array
                load(fullfile(processedDir,'pro_behaviourData.mat'))
                bData = pro_behaviourData{1};
                load(fullfile(processedDir,'pro_trialData.mat'))
                tData = pro_trialData{1};

                % find jump indices
                [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, 10, 10, 1000);
                % calculate the cell's preferred HD
                prefHead = calcPrefHead(bData, tData);

                % prep data, light smoothing improves IPSP detection
                subtractedOutput = smoothdata(tData.scaledOutput - medfilt1(tData.scaledOutput, 500),'loess',25);
                smooth_subVm = smoothdata(medfilt1(tData.scaledOutput, 25),'loess',20);
                Vm_derivative = gradient(smooth_subVm); 

                % for each jump detect IPSPs
                for j = 1:size(trialJumps,1)
                    jump = trialJumps.jump(j); 

                    jumpWindowIdx = [jump_array(jump,1):jump_array(jump,3)];
                    jumpIdx = jump_array(jump,2);

                    jump_subOut = subtractedOutput(jumpWindowIdx);
                    jump_smoothVm = smooth_subVm(jumpWindowIdx);
                    jump_der = Vm_derivative(jumpWindowIdx);

                    % empirically set values that perform well for IPSP
                    % detection
                    minPeakHeight = 2;
                    negThres = 0.25;
                    satisfied = 0; 
                    while satisfied == 0 
                        % find rapidly decreasing mV
                        [~, locs_neg] = findpeaks(jump_der .* -1 ,'MinPeakHeight', negThres,'MinPeakDistance', 20);
                        % find negative peaks of mV
                        [value_amp, locs_amp] = findpeaks(jump_subOut .* -1 ,'MinPeakHeight', minPeakHeight,'MinPeakDistance', 20);

                        % IPSPs occur when rapid - change in membrane potential
                        % occurs right before a negative peak, taking into account
                        % both improved accuracy 
                        locs = [];
                        ipsp_amp = []; 
                        for i = 1:length(locs_neg)
                            timePoint = locs_neg(i); 
                            tdiff = timePoint - locs_amp; 
                            if ~isempty(tdiff(abs(tdiff) <= 30))
                                amp_idx = find(abs(tdiff) == min(abs(tdiff)));
                                ipsp_amp = cat(1, ipsp_amp, value_amp(amp_idx));
                                locs = cat(1, locs, timePoint);
                            end
                        end

                        ipsp = zeros(size(jump_subOut));
                        ipsp_amp = - ipsp_amp;
                        ipsp(locs) = 1;

                        figure(77);clf;
                        a = subplot(4,1,2);
                        plot(jump_subOut)
                        ylabel('subOutput')
                        yline(0)
                        b = subplot(4,1,1);
                        plot(jump_smoothVm)
                        ylabel('smoothVm')
                        hold on
                        c = subplot(4,1,3);
                        plot(jump_der)
                        ylabel('derivative')
                        yline(0)
                        d = subplot(4,1,4);
                        plot(ipsp)
                        linkaxes([a,b,c,d],'x')

                        adjust = input('adjust spike detection?','s');

                        if strcmp(adjust,'n')
                            satisfied = 1;
                        else
                            minPeakHeight = input('Min peak height? ');
                            negThres = input('Negative Threshold? ');
                        end
                    end

                    timePeriod = 5*1000; %5sec
                    jumpTime = (length(jump_der)-1)/2 + 1 ;


                    preJumpAmps = ipsp_amp(locs < jumpTime & locs > jumpTime - timePeriod);
                    postJumpAmps = ipsp_amp(locs > jumpTime & locs < jumpTime + timePeriod);

                    % count number of IPSPs that occured in the 5 seconds
                    % before and after the cue jump & average amplitude of the
                    % IPSPs
                    preJumpFreq = sum(ipsp(jumpTime - timePeriod - 1:jumpTime - 1))/5;
                    preJumpAngle = sum(bData.angle(jumpIdx - 1010:jumpIdx -10))/1000; 
                    preJumpAmp = sum(preJumpAmps)/length(preJumpAmps);
                    postJumpFreq = sum(ipsp(jumpTime + 1:jumpTime + timePeriod + 1))/5;
                    postJumpAngle = sum(bData.angle(jumpIdx + 10: jumpIdx + 1010))/1000; 
                    postJumpAmp = sum(postJumpAmps)/length(postJumpAmps);

                    jumpSum.folder(jCount) = string(folder); 
                    jumpSum.fly(jCount) = fly; 
                    jumpSum.prefHead(jCount) = prefHead; 
                    jumpSum.jumpNum(jCount) = jump; 
                    jumpSum.preJumpAngle(jCount) = preJumpAngle; 
                    jumpSum.postJumpAngle(jCount) = postJumpAngle; 
                    jumpSum.preJumpFreq(jCount) = preJumpFreq; 
                    jumpSum.postJumpFreq(jCount) = postJumpFreq;
                    jumpSum.preJumpAmp(jCount) = preJumpAmp; 
                    jumpSum.postJumpAmp(jCount) = postJumpAmp; 
                    jCount = jCount + 1; 
                end
            end
        end
    end
    
    fn = fieldnames(jumpSum);
    for k=1:numel(fn)
        temp.(fn{k})= jumpSum.(fn{k})';
    end
    jumpSum = struct2table(temp); 

    % calculate change in ISPS frequency and amplitude along with change in
    % HD difference from cell's preferred HD
    for jump = 1:size(jumpSum,1)
        preAngleDiff(jump) = abs(angdiff(deg2rad(jumpSum.prefHead(jump)), deg2rad(jumpSum.preJumpAngle(jump))));
        postAngleDiff(jump) = abs(angdiff(deg2rad(jumpSum.prefHead(jump)), deg2rad(jumpSum.postJumpAngle(jump))));
        changeAngleDiff(jump) = rad2deg((postAngleDiff(jump)- preAngleDiff(jump))); % - = got closer, + = got farther
        
        changeFreq(jump) = jumpSum.postJumpFreq(jump) - jumpSum.preJumpFreq(jump); % - freq went down, + freq went up
        changeAmp(jump) =  abs(jumpSum.postJumpAmp(jump)) - abs(jumpSum.preJumpAmp(jump)) ;
    end
    
    cells = unique(jumpSum.fly);
    nCol = length(cells); 
    
    colormap = cbrewer2('RdBu',nCol);
    
    for jump = 1:size(jumpSum,1)
        cell_num(jump) = find(jumpSum.fly(jump) == cells); 
        cellColours(jump,:) = colormap(cell_num(jump),:); 
    end
    
    figure();
    set(gcf,'color','w')
    set(gcf,'renderer','painters');
    hold on
    scatter(changeAngleDiff, changeFreq,[],cellColours,'filled')
    [coeff, S] = polyfit(changeAngleDiff, changeFreq,1);
    xLine = linspace(min(changeAngleDiff),max(changeAngleDiff),1000);
    [fitLine, ~] = polyval(coeff,xLine,S);
    [R,p] = corrcoef(changeAngleDiff,changeFreq); 
    annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(p(1,2))], 'EdgeColor','none')
    plot(xLine,fitLine,'k')
    xlim([-190 , 190])
    xlabel('change in distance from preferred heading (deg)')
    ylabel('change in IPSP frequency (Hz)')
    
    figure();
    set(gcf,'color','w')
    set(gcf,'renderer','painters');
    hold on
    scatter(changeAngleDiff, changeAmp,[],cellColours,'filled')
    [coeff, S] = polyfit(changeAngleDiff(~isnan(changeAmp)), changeAmp(~isnan(changeAmp)),1);
    xLine = linspace(min(changeAngleDiff(~isnan(changeAmp))),max(changeAngleDiff(~isnan(changeAmp))),1000);
    [fitLine, ~] = polyval(coeff,xLine,S);
    [R,p] = corrcoef(changeAngleDiff(~isnan(changeAmp)),changeAmp(~isnan(changeAmp))); 
    annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(p(1,2))], 'EdgeColor','none')
    plot(xLine,fitLine,'k')
    xlim([-190 , 190])
    xlabel('change in distance from preferred heading (deg)')
    ylabel('change in IPSP amplitude (mV)')
  
    cells = jumpSum.fly;
    pFreq = anovan(changeFreq,{cells', changeAngleDiff},'model',2,'Varnames',{'cell','changeAngle'},'continuous',2); 
    pAmp = anovan(changeAmp,{cells', changeAngleDiff},'model',2,'Varnames',{'cell','changeAngle'},'continuous',2);  

end