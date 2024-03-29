function function_discrete_total_cell_number(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims, plot_times,mixing_pop_id_hist_sim_run,one_or_two_pop)

%% Plot total cell number mean and variance - averaged over realisations

m_springs_per_cell=1; %default value for number of springs per cell

%Create a matrix to store the N for each sim at the times
N_plot = zeros(num_sims,size(plot_times,2));

%For each time extract the number of cells
for N_plot_timeloop = 1:size(plot_times,2)
    for sim_run1=1:num_sims
        %find the closest time
        timing_diff = t_hist_sim_run{sim_run1}(1:end) - plot_times(N_plot_timeloop);
        timing_diff(timing_diff <0) = inf;
        if sum(timing_diff==inf) == size(timing_diff,2)
            time_index_1 = size(timing_diff,2);
        else
            [~, time_index_1] = min(timing_diff);
        end
        N2_plot(sim_run1,N_plot_timeloop) = sum(cell2mat(mixing_pop_id_hist_sim_run{sim_run1}(time_index_1)))/m_springs_per_cell;
        N_plot(sim_run1,N_plot_timeloop) =  ((size(cell2mat(x_hist_sim_run{sim_run1}(time_index_1)),1)-1)/m_springs_per_cell);
        N1_plot(sim_run1,N_plot_timeloop) =  ((size(cell2mat(x_hist_sim_run{sim_run1}(time_index_1)),1)-1)/m_springs_per_cell)   -N2_plot(sim_run1,N_plot_timeloop);
    end
end

%Determine the mean for each time point
N1_mean = mean(N1_plot,1);
N2_mean = mean(N2_plot,1);
N_mean = mean(N_plot,1);

%Determine the standard deviation for each time point
N1_std = std(N1_plot,0,1);
N2_std = std(N2_plot,0,1);
N_std = std(N_plot,0,1);

if one_or_two_pop ==1 %If one population
    %Plot figure
    figure
    hold on
    errorbar(plot_times,N_mean,N_std,'-s','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue','Color','blue','LineStyle','-')
    xlabel('Time')
    ylabel('Cell Number')
    title('Cell Number with Time')
    box on
    
    %Save figure
    print(gcf,'-depsc2',[filepath_save_figs  'N_evolution_dis.eps']);
    saveas(gcf, [filepath_save_figs  'N_evolution_dis.fig'])
    saveas(gcf, [filepath_save_figs  'N_evolution_dis.jpg'])
    
elseif one_or_two_pop==2 %If two populations
    
    %Plot figure
    figure
    hold on
    errorbar(plot_times,N1_mean,N1_std,'-s','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red','Color','red','LineStyle','-')
    errorbar(plot_times,N2_mean,N2_std,'-s','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue','Color','blue','LineStyle','-')
    xlabel('Time')
    ylabel('Cell Number')
    title('Cell Number with Time')
    legend('N1','N2')
    box on
    
    %Save figure
    %print(gcf,'-depsc2',[filepath_save_figs  'N1N2_evolution_dis.eps']);
    saveas(gcf, [filepath_save_figs  'N1N2_evolution_dis.fig'])
    saveas(gcf, [filepath_save_figs  'N1N2_evolution_dis.jpg'])
end

end