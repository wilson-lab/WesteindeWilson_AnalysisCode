function [flies, flyCount, uniqueFlies] = get_flies_ephys(all_folders)
% organize experimental flies
    folders_temp = [];    
    cnt = 0;
    for f = 1:size(all_folders,1)
        cnt = cnt + 1;
        folders_temp{cnt} = all_folders{f};
        slash = regexp(all_folders{f},'\');
        expIDidx = [slash(6)+1:slash(7)-1];
        expID = char(all_folders{f});
        expID = expID(expIDidx);
        flies(cnt) = string(expID);
    end
    flies = flies';
    uniqueFlies = unique(flies); 
    flyCount = zeros(size(all_folders));
    for fly = 1:length(uniqueFlies)
        flyCount(flies == uniqueFlies(fly)) = fly;
    end
end