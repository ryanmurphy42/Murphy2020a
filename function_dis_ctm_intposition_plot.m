function function_dis_ctm_intposition_plot(t_hist,s_hist ,filepath_load_figs,filepath_save_figs, plot_times, ctm_only)

%% Plot the interface position from the continuum model

%determine the interface position at each plot_time
intposition_pde_pop1 = [];
for time_vector_plot_counter = 1:length(plot_times)
    tcomparing = plot_times(time_vector_plot_counter);
    [~, first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
    intposition_pde_pop1(time_vector_plot_counter) = s_hist(first_compare_index);
end

if ctm_only==0 %if discrete and continuum
    %open the discrete
    openfig([filepath_load_figs '\' 'Intposition_d_multipop.fig'])
    hold on
    %plot the continuum
    plot(plot_times,intposition_pde_pop1,'g','linewidth',3)
    legend('Dis','Ctm')
    
else
    %plot the continuum
    figure
    plot(plot_times,intposition_pde_pop1,'g','linewidth',3)
    legend('Ctm')
    
end

%save the figure
print(gcf,'-depsc2',[filepath_save_figs '\' 'Intposition_dis_ctm' '.eps'])
saveas(gcf,[filepath_save_figs '\' 'Intposition_dis_ctm' '.fig'])
saveas(gcf,[filepath_save_figs '\' 'Intposition_dis_ctm' '.jpg'])

end