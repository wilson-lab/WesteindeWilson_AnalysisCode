 function SummarizeIontoSetStimTrials_flyAveraged(rootDir, pulse_window, sampRate, savePlot)
 
    % INPUTS
    
    %   rootDir: Directory to look for trial folders in, string
    %   pulse_window: Time window to grab before & after pulse (s), int
    %   sampRate: The sampling rate of the data, int
    %   savePlot: whether or not to save the plots or not, boolean

    % get experiment folders
    files = dir(fullfile(rootDir,'*/**'));
    folders = files([files.isdir]);
    count = 1; 
    for ff = 1:length(folders)
        if ~isempty(regexp(folders(ff).name,'fly'))
            flies(count) = folders(ff);
            count = count + 1; 
        end
    end
    
    % setup arrays to collect data from pulse windows across flies &
    % trials, row num = num diff pulse lengths
    allFlies_ScOAve(1:4,1:length(flies)) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    allFlies_forAve(1:4,1:length(flies)) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    allFlies_yawAve(1:4,1:length(flies)) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    allFlies_sideAve(1:4,1:length(flies)) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    
    flyAve_ScOAve(1:4,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    flyAve_forAve(1:4,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    flyAve_yawAve(1:4,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    flyAve_sideAve(1:4,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
    
    grey = [0.75,0.75,0.75];
    red = [0.75, 0, 0, 0.25];
    fly_count = length(flies); 
       
    % Iterate through every fly
    for f = 1:length(flies) 
        count = 1;
        setStimTrials = {}; 
        flyDir = fullfile(flies(f).folder,flies(f).name); 
        alltrials = get_folders_ephys_behaviour(flyDir, 1);
        
        % for each fly collect all trials of the right experiment type & load data
        for t = 1:length(alltrials)
            if ~isempty(regexp(alltrials(t).folder, 'ionto_sCL')) && isempty(regexp(alltrials(t).folder, 'badTrials'))
                setStimTrials{count} = alltrials(t).folder;
                count = count + 1;
            end
        end

        if ~isempty(setStimTrials)
            % collect all trials from a fly
            for trial = 1:length(setStimTrials) 
                trialDir = setStimTrials{trial}(1:end-2); 
                load(fullfile(trialDir,'pro_trialData.mat'));
                all_trialData{trial} = processed_trialData{1};
                load(fullfile(trialDir,'pro_behaviourData.mat'));
                all_behaviourData{trial} = processed_behaviourData{1};

                % collect ionto pulses from trial
                [trialPulses] = detect_iontoPulses(all_behaviourData{trial}, pulse_window, sampRate);
                
                for t = size(trialPulses,1):-1:1
                    % get rid of pulses if the window starts before time = 0
                    if trialPulses.windowStart(t) < 0 
                        trialPulses(t,:) = []; 
                    end
                end
                pulseLengths = unique(trialPulses.pulseLength);
                pulse_table{trial} = trialPulses;
            end
                nPulses = length(pulseLengths);
                ScOAve_trial(1:nPulses,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
                forAve_trial(1:nPulses,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
                yawAve_trial(1:nPulses,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};
                sideAve_trial(1:nPulses,1) = {zeros(pulse_window* 2 *sampRate + 1,1)};

                % Collect pulses of each length 
                for p = 1:length(pulseLengths)
                    for trial = 1:length(setStimTrials) % collect average resp per trial
                        pulse = pulseLengths(p);
                        pulses = pulse_table{trial}(pulse_table{trial}.pulseLength == pulse,:);
                        numPulses = size(pulses,1);

                        if numPulses > 0
                            time = [-pulse_window:1/sampRate:pulse_window];
                            ScOAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                            forAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                            yawAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                            sideAve1 = zeros(pulse_window* 2 *sampRate + 1,1);
                             
                            try % sometimes last pulse window goes past trial time depending on trial & window length
                                for i = 1:numPulses
                                    ScOAve1 = ScOAve1 + all_trialData{trial}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i));
                                    forAve1 = forAve1 + all_behaviourData{trial}.vel_for(pulses.windowStart(i):pulses.windowEnd(i));
                                    yawAve1 = yawAve1 + abs(all_behaviourData{trial}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i)));
                                    sideAve1 = sideAve1 + abs(all_behaviourData{trial}.vel_side(pulses.windowStart(i):pulses.windowEnd(i)));
                                end
                            catch
                                for i = 1:numPulses -1
                                    ScOAve1 = ScOAve1 + all_trialData{trial}.scaledOutput(pulses.windowStart(i):pulses.windowEnd(i));
                                    forAve1 = forAve1 + all_behaviourData{trial}.vel_for(pulses.windowStart(i):pulses.windowEnd(i));
                                    yawAve1 = yawAve1 + abs(all_behaviourData{trial}.vel_yaw(pulses.windowStart(i):pulses.windowEnd(i)));
                                    sideAve1 = sideAve1 + abs(all_behaviourData{trial}.vel_side(pulses.windowStart(i):pulses.windowEnd(i)));
                                end
                            end
                            
                            % average pulse responses within each trial
                            ScOAve1 = ScOAve1/i;
                            forAve1 = forAve1/i;
                            yawAve1 = yawAve1/i;
                            sideAve1 = sideAve1/i;

                        end
                        ScOAve_trial{p} = ScOAve_trial{p}  + ScOAve1;
                        forAve_trial{p}  = forAve_trial{p}  + forAve1;
                        yawAve_trial{p}  = yawAve_trial{p}  + yawAve1;
                        sideAve_trial{p}  = sideAve_trial{p}  + sideAve1;
                    end
                    
                    % average pulse responses within each fly
                    ScOAve_trial{p} = ScOAve_trial{p}/trial;
                    forAve_trial{p}  = forAve_trial{p}/trial;
                    yawAve_trial{p}  = yawAve_trial{p}/trial;
                    sideAve_trial{p}  = sideAve_trial{p}/trial;
                end
            
            allFlies_ScOAve(:,f) = ScOAve_trial;
            allFlies_forAve(:,f) = forAve_trial;
            allFlies_yawAve(:,f) = yawAve_trial;
            allFlies_sideAve(:,f) = sideAve_trial;
            
            
            flyAve_ScOAve = cellfun(@plus, flyAve_ScOAve, allFlies_ScOAve(:,f),'Uni',0);
            flyAve_forAve = cellfun(@plus, flyAve_forAve, allFlies_forAve(:,f),'Uni',0);
            flyAve_yawAve = cellfun(@plus, flyAve_yawAve, allFlies_yawAve(:,f),'Uni',0);
            flyAve_sideAve = cellfun(@plus, flyAve_sideAve, allFlies_sideAve(:,f),'Uni',0);        
                    
        else
            fly_count = fly_count - 1;            
        end
    end
    
    % average pulse responses across flies
    fly_countCell(1:4,1) = {fly_count};
    flyAve_ScOAve = cellfun(@(x,y) x/y, flyAve_ScOAve,fly_countCell, 'UniformOutput',false);
    flyAve_forAve = cellfun(@(x,y) x/y, flyAve_forAve,fly_countCell, 'UniformOutput',false);
    flyAve_yawAve = cellfun(@(x,y) x/y, flyAve_yawAve,fly_countCell, 'UniformOutput',false);
    flyAve_sideAve = cellfun(@(x,y) x/y, flyAve_sideAve,fly_countCell, 'UniformOutput',false);
    
    
    % Create plot for each pulse length showing the average responses & save plot in rootDir
    for p = 1:length(pulseLengths)
        stdArray_ScOAve = [];
        stdArray_vf = []; 
        stdArray_vs = []; 
        stdArray_vy = [];
        fig = figure();
        h2 = subplot(10,1,2:4);
        h3 = subplot(10,1,5:7);
        h4 = subplot(10,1,8:10);
        h1 = subplot(10,1,1);
        hold([h2,h3,h4,h1],'on')       
        for f = 1:length(flies)
            if sum(allFlies_ScOAve{p,f}) ~= 0
                stdArray_ScOAve = [stdArray_ScOAve; allFlies_ScOAve{p,f}'];
                stdArray_vf = [stdArray_vf;allFlies_forAve{p,f}']; 
                stdArray_vs = [stdArray_vs;allFlies_sideAve{p,f}']; 
                stdArray_vy = [stdArray_vy;allFlies_yawAve{p,f}'];
   
            end
        end
        SEM_ScOAve = std(stdArray_ScOAve,[],1) / sqrt(size(stdArray_ScOAve,1));
        SEM_vf = std(stdArray_vf,[],1) / sqrt(size(stdArray_vf,1));
        SEM_vs = std(stdArray_vs,[],1) / sqrt(size(stdArray_vs,1));
        SEM_vy = std(stdArray_vy,[],1) / sqrt(size(stdArray_vy,1));
        
        keepIndex = ~isnan(SEM_ScOAve);
        SEMhigh = [flyAve_ScOAve{p}(keepIndex) + SEM_ScOAve(keepIndex)']'; 
        SEMlow = [flyAve_ScOAve{p}(keepIndex) - SEM_ScOAve(keepIndex)']';
        patch(h1,[time(keepIndex) fliplr(time(keepIndex))],[SEMhigh fliplr(SEMlow)],grey,'FaceAlpha',1,'EdgeColor','none')
        
        keepIndex = ~isnan(SEM_vf);
        SEMhigh = [flyAve_forAve{p}(keepIndex) + SEM_vf(keepIndex)']'; 
        SEMlow = [flyAve_forAve{p}(keepIndex) - SEM_vf(keepIndex)']';
        patch(h2,[time(keepIndex) fliplr(time(keepIndex))],[SEMhigh fliplr(SEMlow)],grey,'FaceAlpha',1,'EdgeColor','none')
        
        keepIndex = ~isnan(SEM_vy);
        SEMhigh = [flyAve_yawAve{p}(keepIndex) + SEM_vy(keepIndex)']'; 
        SEMlow = [flyAve_yawAve{p}(keepIndex) - SEM_vy(keepIndex)']';
        patch(h3,[time(keepIndex) fliplr(time(keepIndex))],[SEMhigh fliplr(SEMlow)],grey,'FaceAlpha',1,'EdgeColor','none')
        
        keepIndex = ~isnan(SEM_vs);
        SEMhigh = [flyAve_sideAve{p}(keepIndex) + SEM_vs(keepIndex)']'; 
        SEMlow = [flyAve_sideAve{p}(keepIndex) - SEM_vs(keepIndex)']';
        patch(h4,[time(keepIndex) fliplr(time(keepIndex))],[SEMhigh fliplr(SEMlow)],grey,'FaceAlpha',1,'EdgeColor','none')
        
        plot(h1,time, flyAve_ScOAve{p}, 'k','LineWidth',1.5)
        plot(h2,time, flyAve_forAve{p}, 'k','LineWidth',1.5)
        plot(h3,time, flyAve_yawAve{p}, 'k','LineWidth',1.5)
        plot(h4,time, flyAve_sideAve{p}, 'k','LineWidth',1.5)
        
        ylabel(h1,'mV')
        ylabel(h2, 'vf mm/s')
        ylabel(h3, '|vy| deg/s')
        ylabel(h4, '|vs| mm/s')
        xlabel(h4, 'Time (s)')
        
        ylim(h1,[floor(min(cellfun(@(x) min(x),allFlies_ScOAve),[],[1 2])), ceil(max(cellfun(@(x) max(x),allFlies_ScOAve),[],[1 2]))])
        ylim(h2,[floor(min(cellfun(@(x) min(x),allFlies_forAve),[],[1 2])), ceil(max(cellfun(@(x) max(x),allFlies_forAve),[],[1 2]))])
        ylim(h3,[floor(min(cellfun(@(x) min(x),allFlies_yawAve),[],[1 2])), ceil(max(cellfun(@(x) max(x),allFlies_yawAve),[],[1 2]))])
        ylim(h4, [floor(min(cellfun(@(x) min(x),allFlies_sideAve),[],[1 2])), ceil(max(cellfun(@(x) max(x),allFlies_sideAve),[],[1 2]))])       

        rectangle(h1,'Position', [0, min(ylim(h1)), pulseLengths(p)/1000, max(ylim(h1)) + abs(min(ylim(h1)))], 'FaceColor',red, 'EdgeColor', 'none');
        rectangle(h2,'Position', [0, min(ylim(h2)), pulseLengths(p)/1000, max(ylim(h2)) + abs(min(ylim(h2)))], 'FaceColor',red, 'EdgeColor', 'none');
        rectangle(h3,'Position', [0, min(ylim(h3)), pulseLengths(p)/1000, max(ylim(h3)) + abs(min(ylim(h3)))], 'FaceColor',red, 'EdgeColor', 'none');
        rectangle(h4,'Position', [0, min(ylim(h4)), pulseLengths(p)/1000, max(ylim(h4)) + abs(min(ylim(h4)))], 'FaceColor',red, 'EdgeColor', 'none');
        sgtitle([num2str(pulseLengths(p)),' ms Pulse'])
        linkaxes([h1,h2,h3,h4],'x')
        set(gcf,'color','w')
        hallbutlast = [h1,h2,h3];
        set(hallbutlast,'XTickLabel',[])
        set(h1, 'YColor','k')
        set(hallbutlast, 'XColor','none')
        box off
        set(fig,'units','normalized','position',[0.25 0.15 0.25 0.75])
        set(fig,'Renderer','painters')
        if savePlot
            saveas(fig,fullfile(rootDir,[num2str(pulseLengths(p)),'msPulse_expAve_SEM.svg']))
        end
    end
    

end