function [mean_bin, centers, variance] = binData(activity, behaviour, edges)

% returns for each bin, specificed by the edges which bin each data index
% lies in (bin) as well as the number of data points in each bin (N)
[N, edges, bin] = histcounts(behaviour(:,1), edges);
% sums datapoints belonging to the same bin & saves the values in temp
temp = accumarray(bin+1, activity, [length(edges) 1]);
temp_var = accumarray(bin+1, activity, [length(edges) 1],@var);
variance = temp_var(2:end); 
temp_std = sqrt(temp_var); % std = sqrt(variance)
SEM = temp_std(2:end)./sqrt(N)'; % SEM = std(data)/sqrt(num data points)
% finds the average of each bin value
mean_bin = bsxfun(@rdivide, temp(2:end), N');
% finds the center value of each bin
centers = edges(1:end-1)+diff(edges)/2;

end