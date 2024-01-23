function [expID] = get_expID(folder)
    
    %% Default settings
    arguments
        folder char
    end
    
    % Scan file names in experiment directory to get the expID
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end   
    
    try
        imgDataFiles = dir(fullfile(folder, '*trial*.tif'));
        expID = regexprep(imgDataFiles(1).name,'_.*','');
    catch
        [~, currDir] = fileparts(folder);
        expID = regexprep(currDir,'_.*','');
    end
    
end