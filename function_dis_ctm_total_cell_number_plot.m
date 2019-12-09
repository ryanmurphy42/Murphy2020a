function function_dis_ctm_total_cell_number_plot(t_hist, q_hist,s_hist,L,dx, filepath_load_figs,filepath_save_figs, plot_times, ctm_only,one_or_two_pop)

%% Plot the total cell number


if one_or_two_pop==2 %if two populations
    
    
    x_discretisation = 0:dx:L; %Spatial discretisation
    
    %Calculate the number of cells in each population
    cell_number_pde_pop1 = [];
    for time_vector_plot_counter = 1:length(plot_times)
        tcomparing = plot_times(time_vector_plot_counter);
        [~, first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        
        %find the closest node
        xcomparing = s_hist(first_compare_index);
        [~, first_compare_indexx]   = min(abs(x_discretisation-xcomparing));
        
        x_discretisation(first_compare_indexx);
        
        nodes_with_pop1= first_compare_indexx;
        
        cell_number_pde_pop1(time_vector_plot_counter) = trapz(0:dx:((first_compare_indexx-1)*dx),q_hist(1:(nodes_with_pop1),first_compare_index));
    end
    
    
    cell_number_pde_pop2 = [];
    for time_vector_plot_counter = 1:length(plot_times)
        tcomparing = plot_times(time_vector_plot_counter);
        [~, first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        
        
        %find the closest node
        xcomparing = s_hist(first_compare_index);
        [~, first_compare_indexx]   = min(abs(x_discretisation-xcomparing));
        
        x_discretisation(first_compare_indexx);
        
        nodes_with_pop2= first_compare_indexx;
        
        total_nodes = length(q_hist(:,first_compare_index));
        
        if total_nodes == nodes_with_pop2
            cell_number_pde_pop2(time_vector_plot_counter) = 0;
        else
            cell_number_pde_pop2(time_vector_plot_counter) = trapz(((nodes_with_pop2-1)*dx):dx:L,q_hist((nodes_with_pop2:end),first_compare_index));
        end
    end
    
    %% Plots
    
    
    if ctm_only==0 %if discrete and continuum 
        
        %open the discrete figure
        openfig([filepath_load_figs '\' 'N1N2_evolution_dis.fig'])
        hold on
        
        %plot the total cell number for population 1 from the continuum model
        plot(plot_times,cell_number_pde_pop1,'m','linewidth',3)
        %plot the total cell number for population 2 from the continuum model
        plot(plot_times,cell_number_pde_pop2,'c','linewidth',3)
        
        legend('Dis N1','Dis - N2','Ctm - N1','Ctm - N2')
        
        %save figure
        print(gcf,'-depsc2',[filepath_save_figs '\' 'N1N2_evolution_dis_ctm' '.eps'])
        saveas(gcf,[filepath_save_figs '\' 'N1N2_evolution_dis_ctm' '.fig'])
        saveas(gcf,[filepath_save_figs '\' 'N1N2_evolution_dis_ctm' '.jpg'])
        
        
    else
         
         
        figure
        %plot the total cell number for population 1 from the continuum model
        plot(plot_times,cell_number_pde_pop1,'g','linewidth',3)
        hold on
        %plot the total cell number for population 2 from the continuum model
        plot(plot_times,cell_number_pde_pop2,'g','linewidth',3)
        
        %Save figures
        print(gcf,'-depsc2',[filepath_save_figs '\' 'N1N2_evolution_dis' '.eps'])
        saveas(gcf,[filepath_save_figs '\' 'N1N2_evolution_dis' '.fig'])
        saveas(gcf,[filepath_save_figs '\' 'N1N2_evolution_dis' '.jpg'])
        
    end
    
elseif one_or_two_pop==1 %if one population
    
    %Calculate the total number of cells
    cell_number_pde_totalpop = [];
    for time_vector_plot_counter = 1:length(plot_times)
        tcomparing = plot_times(time_vector_plot_counter);
        [~, first_compare_index]   = min(abs(t_hist(1:end) -tcomparing));
        cell_number_pde_totalpop(time_vector_plot_counter) = trapz(0:dx:L,q_hist(:,first_compare_index));
    end
    
    
    
    if ctm_only==0 %if discrete and continuum
        
        %open the discrete
        openfig([filepath_load_figs '\' 'N_evolution_dis.fig'])
        hold on
        
        %plot the continuum
        plot(plot_times,cell_number_pde_totalpop,'g','linewidth',3)
        
        legend('Dis','Ctm')
        
        %save the figures
        print(gcf,'-depsc2',[filepath_save_figs '\' 'N_evolution_dis_ctm' '.eps'])
        saveas(gcf,[filepath_save_figs '\' 'N_evolution_dis_ctm' '.fig'])
        saveas(gcf,[filepath_save_figs '\' 'N_evolution_dis_ctm' '.jpg'])
        
        
    elseif ctm_only==1 %continuum only
        
        %plot the total cell number
        figure
        plot(plot_times,cell_number_pde_totalpop,'g','linewidth',3)
        legend('Dis','Ctm')
        
        %save the figure
        print(gcf,'-depsc2',[filepath_save_figs '\' 'N_evolution_ctm' '.eps'])
        saveas(gcf,[filepath_save_figs '\' 'N_evolution_ctm' '.fig'])
        saveas(gcf,[filepath_save_figs '\' 'N_evolution_ctm' '.jpg'])
        
    end
    
end



end