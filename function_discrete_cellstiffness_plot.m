function function_discrete_cellstiffness_plot(t_hist_sim_run, x_hist_sim_run,k_hist_sim_run, L, ...
    filepath_save_figs,num_sims, plot_times, density_evaluation_x_dx)

%% Plot the cell stiffness snapshots which is averaged over realisations


density_evaluation_x = 0:density_evaluation_x_dx:L; %positions to calculate the density in the domain

%For each time determine the cell stiffness accross the domain
%For each simulation determine the cell stiffness

%Create a matrix to store the cell stiffness at each position specified and at each plot_time
k_plot_position = zeros(num_sims,size(density_evaluation_x,2));


%Calculate times for plot
k_plot_timeloop_index = zeros(length(length(plot_times)),1);
for jj=1:length(plot_times)
    [~, inda] = min(abs(plot_times(jj) - plot_times));
    k_plot_timeloop_index(jj) = inda;
end

%for each time calculate the mean and standard deviation of the cell stiffness and plot this


for k_plot_timeloop = 1:size(plot_times,2)
    if ismember(k_plot_timeloop, k_plot_timeloop_index) ==1
        for sim_run1=1:num_sims
            %find the closest time
            timing_diff = t_hist_sim_run{sim_run1}(1:end) - plot_times(k_plot_timeloop);
            
            
            timing_diff(timing_diff <0) = inf;
            if sum(timing_diff==inf) == size(timing_diff,2)
                time_index_1 = size(timing_diff,2);
            else
                [~, time_index_1] = min(timing_diff);
            end
            
            %find the positions and k
            if time_index_1 == 1
                x_den_lookup = x_hist_sim_run{sim_run1}{time_index_1};
                k_lookup = k_hist_sim_run{sim_run1}{time_index_1};
            else
                x_den_lookup = x_hist_sim_run{sim_run1}{time_index_1 -1};
                k_lookup = k_hist_sim_run{sim_run1}{time_index_1 -1};
            end
            
            for q_eval_loop=1:size(density_evaluation_x,2)
                
                %determine which cell
                cum_sum_pos_x_plot1_iso = x_den_lookup - density_evaluation_x(q_eval_loop);
                cum_sum_pos_x_plot1_iso(cum_sum_pos_x_plot1_iso < 0) = inf;
                
                [~, pos_index_1] = min(cum_sum_pos_x_plot1_iso);
                
                pos_index_1 = pos_index_1-1;
                
                %determine the cell stiffness at the position
                if pos_index_1 == 0
                    k_in_loop = k_lookup(1);
                elseif  pos_index_1 ==  size(x_den_lookup,1)
                    k_in_loop = k_lookup(end);
                else
                    k_in_loop = k_lookup(pos_index_1);
                end
                %save the cell stiffness
                k_plot_position(sim_run1, q_eval_loop) = k_in_loop;
            end
            
        end
        
        %Determine the mean for each time point
        k_plot_position_mean = mean(k_plot_position,1);
        
        %Determine the standard deviation for each time point
        k_plot_position_std = std(k_plot_position,0,1);
        
        %for a clear figure - only plot standard deviation for every 50 times.
        valToKeep = 1:50:length(k_plot_position_std); 
        for i=1:length(k_plot_position_std)
            if ~any(valToKeep == i)
                k_plot_position_std(i) = 0;
            end
        end
        
        %plot figure
        figure
        errorbar(density_evaluation_x,k_plot_position_mean,k_plot_position_std,'-s','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','blue')
        xlabel('Position')
        ylabel('Cell stiffness')
        title(['Cell stiffness @ ' num2str(plot_times(k_plot_timeloop))])
        
        %save figure
        %print(gcf,'-depsc2',[filepath_save_figs 'Cellstiffness_evolution_dis_' num2str(k_plot_timeloop) '_std.eps']);
        saveas(gcf,[filepath_save_figs 'Cellstiffness_evolution_dis_' num2str(k_plot_timeloop) '_std.fig'])
        saveas(gcf,[filepath_save_figs 'Cellstiffness_evolution_dis_' num2str(k_plot_timeloop) '_std.jpg'])
        
        %Close figure
        close all
    end
    
    
end


end
