function CalcBumpParams(rootDir)
% iterates through trials & fits a sinusoid to PFL2 activity across
% brain space, extracts various parameters of interest from the fitted
% sinusoid & saves to trial data file. 

% INPUT

% rootDir: directory to look for experimental folders

    folders = get_folders(rootDir, 1, 0);
    count = 1; 

    for f = 1:size(folders,1)
        folder = folders(f).folder; 
        if contains(folder,'FB') || contains(folder,'PB')
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

            try

                expID = get_expID(folder);
                expList = {expID};

                % Load metadata 
                [~, trialMd] = load_metadata(expList, folder);

                % Load imaging data
                roiData = load_roi_data(expList, folder);

                processedData_dir = fullfile(folder,'processed_data');
                try
                    numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum));
                catch
                    numTrials = 1; 
                end

                % iterate through trials
                for nTrial = 1:numTrials
                    % load in 
                    clear ZData
                    data_filelist = dir(processedData_dir);
                    for files = 1:length(data_filelist)
                        if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                            load(fullfile(processedData_dir,data_filelist(files).name));
                        end  
                    end
                    load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']));

                    x_range = [0,2*pi];
                    activity = [];
                    for roi = 1:size(ZData,1)
                        activity(roi,:) = ZData.Z{roi};
                    end

                    minZ = min(activity,[],[1, 2]);
                    if minZ < 0 % The Z score better compensates for the unequal baseline F across ROIs but the fit doesn't work with negative values
                        activity = activity + abs(minZ);
                    end
                    
                    % fit sinusoid to neural activity at each timepoint
                    bump_params = fit_sinusoid(activity, x_range, 1);

                    %% Offset calculation

                    %we will compute and define the offset as the
                    %circular distance between the bump's position and the fly's heading
                    fly_pos_rad = wrapToPi(deg2rad(ftT_down.cueAngle{1}));
                    offset = wrapTo180(rad2deg(circ_dist(bump_params.pos_rad,-fly_pos_rad)));                
                    bump_params.offset = offset;

                    save(fullfile(processedData_dir,['bump_parameters_down_Trial00',num2str(nTrial),'.mat']),'bump_params')

                    % option to plot bump parameters to check
                    %plotBumpParameters(folder,dffData, ftT, nTrial, bump_params)
                end
            catch
                failedFolders{count} = folder; 
                disp(['Folder failed: ',folder])
                count= count + 1; 
            end
        end
    end
end
