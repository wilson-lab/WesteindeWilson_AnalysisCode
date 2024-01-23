function [bump_params, model_data, x_grid] = fit_sinusoid(activity, x_range, to_plot)
%%% fit_von_Mises(x_range, dff_data, to_plot)
%%% based on fitVonMises (Mel Basnak)
%%%
%%% Args:
%%%    activity (struct): contains dF_F
%%%    x_range (1x2 double array): range of the angles covered (rad)
%%%    to_plot (logical): whether to plot for debugging or not
%%%
%%% Returns:
%%%    bump_params (struct): bump pos (rad, idx), bump magnitude, bump width, and
%%%    adjusted R^2
%%%
%%% Tatsuo Okubo
%%% 2022/03/23
%%% Elena Westeinde
%%% 2022/11 - modified to fit a sinusoid function
assert(x_range(1) == 0, 'x_range(1) needs to be 0')

dff_data = activity';
T = size(dff_data, 1);

n_centroid = size(dff_data, 2);
mid_dist_2pi_wrap = linspace(x_range(1), x_range(2), n_centroid + 1)';  % + 1 for wrapping around the full circle  
mid_dist_2pi = mid_dist_2pi_wrap(1:end-1);

%% set up optimization parameters 

fo = fitoptions(...
    'Method', 'NonlinearLeastSquares',...
    'Lower', [0,-inf, -inf],... % [a,b,c,k,u]s
    'Upper',[inf,inf,inf],...
    'StartPoint', [1, 0, 0]);

%ft = fittype('a*exp(k*cos(x-u))+c', 'options', fo);
ft = fittype('a*sin((x-u))+c', 'options', fo);
% a = amp
% u = horizontal pos
% c = vertical pos

%% initialize arrays
bump_pos = nan(T, 1);
bump_mag = nan(T, 1);
bump_width = nan(T, 1);
adj_rs = nan(T, 1);

%%
x_grid = linspace(x_range(1), x_range(2), 100);  % for plotting the best fit curve
XL = [0, ceil(x_range(2)/ (2 * pi)) * 2 * pi];
YL(1) = min(min(dff_data));
YL(2) = max(max(dff_data));       

f = waitbar(0, 'Analyzing...');

%%

for t = 1:T  % for each time point
    waitbar(t / T, f, sprintf('Progress: %d %%', floor((t / T) * 100)));
    data_to_fit = dff_data(t, :);
    [model_data, gof] = fit(mid_dist_2pi, data_to_fit', ft, 'MaxIter', 20000, 'MaxFunEvals', 20000);
    adj_rs(t) = gof.adjrsquare;
    rs(t) = gof.rsquare;
    
    % get the coefficient values.
    coeff_names = coeffnames(ft);
    coefficientValues = coeffvalues(model_data);
    a = coefficientValues(strcmp(coeff_names,'a'));
    u = coefficientValues(strcmp(coeff_names,'u'));
    c = coefficientValues(strcmp(coeff_names,'c'));
    
    % calculate bump parameters
    pos_temp = x_grid(feval(model_data, x_grid) == max(feval(model_data, x_grid)));
   if any(pos_temp == 0) && any(round(pos_temp,2) == 6.28)
        bump_pos(t) = 0; 
   elseif length(pos_temp) > 1
       bump_pos(t) = nan; 
   else
        bump_pos(t) = pos_temp;
   end
    bump_amp(t) = a;
    
    % plot results if necessary
    if to_plot == true
        figure(10); clf;
        set(gcf, 'position', [500, 400, 800, 600])
        subplot(3,1,1:2)
        plot(mid_dist_2pi, dff_data(t,:)', 'k.-')
        hold on    
        plot(x_grid, feval(model_data, x_grid), 'r-', 'linewidth', 2);         % plot the fitted curve
        plot(bump_pos(t), feval(model_data,bump_pos(t)), 'ro'); % add the bump position estimate
        xlabel('ROI location (rad)');
        ylabel('dF/F');
        title(['Frame #',num2str(t), ' Adj R^2 =', num2str(gof.adjrsquare)]);
        box off
        set(gca, 'tickdir', 'out', 'xtick', 0:pi:4 * pi, 'xticklabels', {'0', '\pi', '2\pi', '3 \pi', '4 \pi'})
        xlim(XL);
        ylim(YL)
        line([2 * pi, 2 * pi], YL, 'color', 'black', 'linestyle', ':')
        
        subplot(3,1,3)
        text(0,0.5,['Bump amplitude = ',num2str(bump_amp(t))]);
        hold on
        text(0.25,0,['Bump pos = ',num2str(bump_pos(t))]); axis off
    end

end

close(f)

%% save the variables in a structure array
bump_params.pos_rad = bump_pos;
bump_params.pos_idx = (n_centroid - 1) * bump_pos / x_range(2) + 1;  % bump location in terms of ROI indices
bump_params.amp = bump_amp;
bump_params.adj_rs = adj_rs;
bump_params.rs = rs;