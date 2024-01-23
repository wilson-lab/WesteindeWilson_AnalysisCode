function noMov_jumps = Find_0velJumps(rootDir)
% identifies cue jumps where the fly didn't move in the 5 seconds before 
% or after the jump

    varNames = ["folder","expID","trialNum","jump"];
    varTypes = ["string","string","double","double"];
    noMov_jumps = table('Size',[1 4],'VariableNames', varNames,'VariableTypes',varTypes); 

    folders = get_folders_ephys(rootDir);    
    jcount = 1; 

    for f = 1:size(folders,1)

        folder = folders(f).folder; 
        if strcmp(folder(end),'.')
            folder = folder(1:end-2); 
        end

        % find file path & load data
        cellIdx = regexp(folder,'cell');
        idx = regexp(folder,'\');
        idff = idx - cellIdx;
        iIdx = find(abs(idff) == min(abs(idff)));
        expID = folder(idx(iIdx -1)+1:idx(iIdx)-1);
        trialNum = folder(end);

        processedDir = fullfile(folder,'processedData');
        clear bData tData jump_array
        
        load(fullfile(processedDir,'pro_behaviourData.mat'))
        bData = pro_behaviourData{1};
        load(fullfile(processedDir,'pro_trialData.mat'))
        tData = pro_trialData{1};

        % detect jumps
        [jump_array, ~, ~] = detect_jumps_ephys(bData.frY, 5,5, 1000);

        for jump = 1:size(jump_array,1)
            jumpWindowIdx = [jump_array(jump,1):jump_array(jump,3)];

            if sum(jumpWindowIdx < 1) == 0 && jumpWindowIdx(end) < size(bData,1)
                % cumulative velocity
                total_mov_mm = abs(bData.vel_for(jumpWindowIdx)) + abs(bData.vel_side(jumpWindowIdx)) + abs(deg2rad(bData.vel_yaw(jumpWindowIdx))*4.5);
                no0vel_idx = find(total_mov_mm > 3);

                total_mov_mm_all = abs(bData.vel_for) + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5);
                no0vel_idx_all = find(total_mov_mm_all > 3);
                timeMov = length(no0vel_idx_all)/1000;

                % if fly didn't move around the jump but did run for at
                % least 60 seconds elsewhere during the trial, collect.
                % Screen for general health
                if isempty(no0vel_idx) && timeMov > 60
                    noMov_jumps.folder(jcount) = string(folder);
                    noMov_jumps.expID(jcount) = string(expID);
                    noMov_jumps.trialNum(jcount) = str2double(trialNum);
                    noMov_jumps.jump(jcount) = jump;      
                    jcount = jcount + 1; 

                    figure(22);clf;
                    subplot(3,1,1)
                    plot(bData.vel_for(jumpWindowIdx))
                    subplot(3,1,2)
                    plot(bData.vel_yaw(jumpWindowIdx))
                    subplot(3,1,3)
                    plot(bData.vel_side(jumpWindowIdx))
                end
            end
        end
    end
end