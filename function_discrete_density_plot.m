function function_discrete_density_plot(t_hist_sim_run, x_hist_sim_run, L, ...
    filepath_save_figs,num_sims, plot_times, density_evaluation_x_dx)

%% Plot the density snapshots which is averaged over realisations

density_evaluation_x = 0:density_evaluation_x_dx:L; %positions to calculate the density in the domain
m_springs_per_cell=1; %default value for 1 spring per cell

%For each time determine the density accross the domain
%For each simulation determine the density

%Create a matrix to store the density at each position specified and at each plot_time
q_plot_position = zeros(num_sims,size(density_evaluation_x,2));


%Calculate times for plot
q_plot_timeloop_index = zeros(length(length(plot_times)),1);
for jj=1:length(plot_times)
    [~, inda] = min(abs(plot_times(jj) - plot_times));
    q_plot_timeloop_index(jj) = inda;
end

%for each time calculate the mean and standard deviation of the density and plot this

for q_plot_timeloop = 1:size(plot_times,2)
    if ismember(q_plot_timeloop, q_plot_timeloop_index) ==1
        for sim_run1=1:num_sims
            %find the closest time
            timing_diff = t_hist_sim_run{sim_run1}(1:end) - plot_times(q_plot_timeloop);
            
            timing_diff(timing_diff <0) = inf;
            if sum(timing_diff==inf) == size(timing_diff,2)
                time_index_1 = size(timing_diff,2);
            else
                [~, time_index_1] = min(timing_diff);
            end
            
            %find the positions
            if time_index_1 == 1
                x_den_lookup = x_hist_sim_run{sim_run1}{time_index_1};
            else
                x_den_lookup = x_hist_sim_run{sim_run1}{time_index_1 -1};
            end
            
            for q_eval_loop=1:size(density_evaluation_x,2)
                
                %determine which cell
                cum_sum_pos_x_plot1_iso = x_den_lookup - density_evaluation_x(q_eval_loop);
                cum_sum_pos_x_plot1_iso(cum_sum_pos_x_plot1_iso < 0) = inf;
                
                [~, pos_index_1] = min(cum_sum_pos_x_plot1_iso);
                
                %caclulate the density at the position
                if pos_index_1 == 1
                    den_in_loop = 1/((x_den_lookup(2) - x_den_lookup(1))*m_springs_per_cell);
                elseif  pos_index_1 ==  size(x_den_lookup,1)
                    den_in_loop = 1/(((x_den_lookup(size(x_den_lookup,1) ) - x_den_lookup(size(x_den_lookup,1) -1)))*m_springs_per_cell);
                else
                    den_in_loop = 1/(((x_den_lookup(pos_index_1) - x_den_lookup(pos_index_1-1)))*m_springs_per_cell);
                end
                
                %save the density
                q_plot_position(sim_run1, q_eval_loop) = den_in_loop;
            end
            
        end
        
        %Determine the mean for each time point
        q_plot_position_mean = mean(q_plot_position,1);
        
        %Determine the standard deviation for each time point
        q_plot_position_std = std(q_plot_position,0,1);
        
       
        %for a clear figure - only plot standard deviation for every 50 times.
        valToKeep = 1:50:length(q_plot_position_std); 
        for i=1:length(q_plot_position_std)
            if ~any(valToKeep == i)
                q_plot_position_std(i) = 0;
            end
        end
        
        %Plot figure
        figure
        errorbar(density_evaluation_x,q_plot_position_mean,q_plot_position_std,'-s','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','blue')
        xlabel('Position')
        ylabel('Density')
        title(['Density @ ' num2str(plot_times(q_plot_timeloop))])
        
        %Save figure
        %print(gcf,'-depsc2',[filepath_save_figs 'Density_evolution_dis_' num2str(q_plot_timeloop) '_std.eps']);
        saveas(gcf,[filepath_save_figs 'Density_evolution_dis_' num2str(q_plot_timeloop) '_std.fig'])
        saveas(gcf,[filepath_save_figs 'Density_evolution_dis_' num2str(q_plot_timeloop) '_std.jpg'])
        
        %Close image
        close all
    end
    
    
end


end
