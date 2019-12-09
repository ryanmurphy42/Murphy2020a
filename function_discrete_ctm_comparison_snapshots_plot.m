function function_discrete_ctm_comparison_snapshots_plot(t_hist, q_hist,k_hist,...
    L,dx, filepath_save_figs,  colouring,ctm_only, plot_times,filepath_load_figs)

%% Plot snapshots from continuum model

%colouring: 1-density, 2- cell stiffness

if colouring == 1 %density
    
    for t_loop = 1:length(plot_times)
        
        if ctm_only == 0
            %open the discrete figure
            openfig([filepath_load_figs  'Density_evolution_dis_' num2str(t_loop) '_std.fig'])
        elseif ctm_only == 1
            figure
        end
        %find the nearest ctm plot time
        [~, first_compare_index]   = min(abs(t_hist(1:end) - plot_times(t_loop)));
        
        %overlay the ctm result
        hold on
        plot(0:dx:L , q_hist(:,first_compare_index),'LineWidth',2, 'color','g')
        
        %update figure properties
        if ctm_only  == 0
            legend('Dis', 'Ctm')
            
            %save figures
            %print(gcf,'-depsc2',[filepath_save_figs '\' 'Density_evolution_dis_ctm' num2str(t_loop) '_std.eps'])
            saveas(gcf,[filepath_save_figs '\' 'Density_evolution_dis_ctm_' num2str(t_loop) '_std.fig'])
            saveas(gcf,[filepath_save_figs '\' 'Density_evolution_dis_ctm_' num2str(t_loop) '_std.jpg'])
            
        elseif ctm_only == 1
            legend('Ctm')
            title(['Density plot ctm at t=' num2str(t_hist(first_compare_index))])
            
            %save figures
            %print(gcf,'-depsc2',[filepath_save_figs '\' 'Density_evolution_ctm_ ' num2str(t_loop) '.eps'])
            saveas(gcf,[filepath_save_figs '\' 'Density_evolution_ctm_' num2str(t_loop) '.fig'])
            saveas(gcf,[filepath_save_figs '\' 'Density_evolution_ctm_' num2str(t_loop) '.jpg'])
        end
        close all
    end
    
    
    
    
elseif colouring  == 2 %cell stiffness
    
    
    for t_loop = 1:length(plot_times)
        
        if ctm_only == 0
            %open the discrete figure
            openfig([filepath_load_figs '\' 'Cellstiffness_evolution_dis_' num2str(t_loop) '_std.fig'])
        elseif ctm_only == 1
            figure
        end
        %find the nearest ctm plot time
        [~, first_compare_index]   = min(abs(t_hist(1:end) -plot_times(t_loop)));
        
        %overlay the ctm result
        hold on
        plot(0:dx:L , k_hist(:,first_compare_index),'LineWidth',2, 'color','g')
        
        %update figure properties
        if ctm_only  == 0
            legend('Dis', 'Ctm')
            %save figures
            %print(gcf,'-depsc2',[filepath_save_figs '\' 'Cellstiffness_evolution_dis_ctm_' num2str(t_loop) '_std.eps'])
            saveas(gcf,[filepath_save_figs '\' 'Cellstiffness_evolution_dis_ctm_' num2str(t_loop) '_std.fig'])
            saveas(gcf,[filepath_save_figs '\' 'Cellstiffness_evolution_dis_ctm_' num2str(t_loop) '_std.jpg'])
            
        elseif ctm_only == 1
            legend('Ctm')
            title(['k plot ctm at t=' num2str(t_hist(first_compare_index))])
            %save figures
            %print(gcf,'-depsc2',[filepath_save_figs '\' 'Cellstiffness_evolution_ctm_' num2str(t_loop) '.eps'])
            saveas(gcf,[filepath_save_figs '\' 'Cellstiffness_evolution_ctm_' num2str(t_loop) '.fig'])
            saveas(gcf,[filepath_save_figs '\' 'Cellstiffness_evolution_ctm_' num2str(t_loop) '.jpg'])
        end
        
        close all
    end
    
    
end



end
