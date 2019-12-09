function function_characteristics_plot(t_hist, x_hist, k_hist,G_hist, mixing_pop_id_hist, a_hist, L,...
        filepath_save_figs,  colouring,  xticks_vec, yticks_vec,....
        cell_with_proliferation_event_hist,cell_with_death_event_hist,cell_event_hist,prolif_death_indicator_hist,sim_run,cell_death_lr_hist,m_springs_per_cell)

    %function_characteristcs_plot - plot characteristics for each realisation sequentially

%t_hist - times recorded
%x_hist - positions
%k_hist - cell stiffness
%a_hist - resting spring length
%G_hist - proliferation, death propensities
%mixing_pop_id_hist 
%L - length of domain
%filepath_save_figs
%colouring %1- density, 2 - cell stiffness, 3-resting cell length, 4-population_id colouring
%xticks_vec - x tickmarks
%yticks_vec- y tick marks


%% check if x_hist changes by more than one in any timestep

x_hist_double_change =[];
for t_step = 2:size( t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step == size_now -m_springs_per_cell
    elseif size_prev_step == size_now + m_springs_per_cell
    elseif size_prev_step == size_now
    else
        x_hist_double_change = [x_hist_double_change,t_hist(t_step) ];
    end
end


%% Calculate the proliferation times

prolif_times =[];
for t_step = 2:size( t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step < size_now
        prolif_times = [prolif_times,t_hist(t_step) ];
    end
end

%% Calculate the cell death times

death_times =[];
for t_step = 2:size(t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step > size_now
        death_times = [death_times,t_hist(t_step) ];
    end
end

%% Calculate the cell event times

event_times =[];
for t_step = 2:size(t_hist,2)
    if t_step == 2
        size_prev_step = size(x_hist{t_step},1);
    else
        size_prev_step = size_now;
    end
    
    size_now = size(x_hist{t_step},1);
    
    if size_prev_step > size_now %death event
        event_times = [event_times,t_hist(t_step) ];
    elseif size_prev_step < size_now %prolif event
        event_times = [event_times,t_hist(t_step) ];
    end
end

%% Calculate traceback history of particles

ch_plot_ord = {};

t_step_ch=1;
cell_prolif_count=0;
cell_death_count=0;
cell_death_index_count=0;

for t_step_ch = 1:(size(event_times,2)+1)
    
    if t_step_ch == 1 %if initial time label the cells
        ch_plot_ord{t_step_ch} = 1:size(x_hist{1},1);
    else %if not the initial time update the cell ordering
        
        %the cell ordering from the previous time step
        int_ch_plot_ord = ch_plot_ord{t_step_ch-1};
        
        %cell event
        cell_ev = cell_event_hist(t_step_ch-1);
        
        %determine whether prolif or death event
        if prolif_death_indicator_hist(t_step_ch-1) ==1 %if a proliferation event
            
            cell_prolif_count = cell_prolif_count + 1; %increase loop counter for cell prolif
            
            cell_prolif_ev = cell_with_proliferation_event_hist(cell_prolif_count);
            
            cell_prolif_ev_m = cell_prolif_ev*m_springs_per_cell;
            
            for int_ch_plot_ord_loop = 1:size(int_ch_plot_ord,2)
                
                if int_ch_plot_ord(int_ch_plot_ord_loop) <= cell_prolif_ev_m %if to the left index the same
                    int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop);
                    
                elseif int_ch_plot_ord(int_ch_plot_ord_loop) > cell_prolif_ev_m  %if to the right index increase
                    int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop) + 1;
                end
                
            end
            %add on the cell which proliferated
            int_ch_plot_ord = [int_ch_plot_ord,cell_prolif_ev+1];
            ch_plot_ord{t_step_ch} = int_ch_plot_ord;
            
        elseif prolif_death_indicator_hist(t_step_ch-1) == 2  %if a death event
            
            cell_death_count = cell_death_count + 1; %increase loop counter to cell death
            
            cell_death_ev = cell_with_death_event_hist(cell_death_count);
            
            for int_ch_plot_ord_loop = 1:size(int_ch_plot_ord,2)
                
                if int_ch_plot_ord(int_ch_plot_ord_loop) > 0
                    if int_ch_plot_ord(int_ch_plot_ord_loop) <= cell_death_ev %if to the left index the same
                        int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop);
                        
                    elseif int_ch_plot_ord(int_ch_plot_ord_loop) > cell_death_ev+1  %if to the right index decrease
                        int_ch_plot_ord(int_ch_plot_ord_loop) = int_ch_plot_ord(int_ch_plot_ord_loop) -1 ;
                        
                    elseif int_ch_plot_ord(int_ch_plot_ord_loop) == cell_death_ev+1  %if to the right index decrease
                        cell_death_index_count = cell_death_index_count -1;
                        int_ch_plot_ord(int_ch_plot_ord_loop) = cell_death_index_count ;
                    end
                end
            end
            ch_plot_ord{t_step_ch} = int_ch_plot_ord;
        end
        
    end
end

%%
% traceback from the final state to determine the cells at each timestep
% determine the starting point for each characteristic

x_order_hist = {};
traceback_Nfin_loop=1;
prolif_counter=0;
for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
    inter_hist = [];
    if traceback_Nfin_loop <= size(x_hist{1},1)
        for k_loop = 1:size(ch_plot_ord,2)
            inter_hist = [inter_hist,ch_plot_ord{k_loop}(traceback_Nfin_loop)];
        end
        x_order_hist{traceback_Nfin_loop} =   inter_hist;
    else
        %determine number of death events before cell boundary was generated then know the start time of this event
        
        %a new characteristic for each proliferation event
        %at what position did the proliferation event occur
        
        prolif_counter=prolif_counter+1;
        
        cum_sum_prolif_death_indicator_hist  =cumsum(prolif_death_indicator_hist==1);
        
        [~,time_begin_for_traj]= min(abs(cum_sum_prolif_death_indicator_hist-prolif_counter));
        
        for k_loop=1:(size(ch_plot_ord,2) -time_begin_for_traj)
            inter_hist = [inter_hist,ch_plot_ord{size(ch_plot_ord,2) - k_loop + 1}(traceback_Nfin_loop)];
        end
        x_order_hist{traceback_Nfin_loop} =   fliplr(inter_hist);
    end
end

%% determine the start times for each of the characteristic (based on the final positions)
%% and end times of each characteristics

%for each cell determine the positions
x_hist_traceback = {};
t_hist_traceback = {};
q_hist_traceback = {};
k_hist_traceback = {};
a_hist_traceback = {};
G_net_hist_traceback = {};
mixing_pop_id_hist_traceback = {};

for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2) %for each of the cells
    x_hist_follow_boundary =[];
    t_hist_follow_boundary =[];
    q_hist_follow_boundary =[];
    k_hist_follow_boundary =[];
    a_hist_follow_boundary =[];
    G_net_hist_follow_boundary = [];
    mixing_pop_id_hist_follow_boundary = [];
    
    for t_step = 1:size(t_hist,2) %for each timestep
        
        t_val = t_hist(t_step);
        
        %determine where t_val sits within event_times
        event_times_inf = event_times - t_val;
        event_times_inf(event_times_inf <= 0) = inf;
        
        if sum(event_times_inf== inf) == size(event_times_inf,2)
            event_time_pos = size(event_times,2)+1;
        else
            [~, event_time_pos] = min(event_times_inf);
        end
        
        %does the cell exist at this time
        if event_time_pos >= size(event_times,2) - size(x_order_hist{traceback_Nfin_loop},2) + 2 %cell exists
            
            event_time_pos_with_prolif = event_time_pos +  (size(x_order_hist{traceback_Nfin_loop},2)-size(x_order_hist{1},2));
            
            ch_plot_order_x = x_order_hist{traceback_Nfin_loop}(event_time_pos_with_prolif);
            %ch_plot_ord{prolif_time_pos}(traceback_Nfin_loop);
            
            if ch_plot_order_x > 0
                
                x_val = x_hist{t_step}(ch_plot_order_x);
                x_hist_follow_boundary = [x_hist_follow_boundary,x_val];
                t_hist_follow_boundary = [t_hist_follow_boundary,t_val];
                
                %density
                if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                    q_val = (1/(x_hist{t_step}(ch_plot_order_x+1) - x_hist{t_step}(ch_plot_order_x)));
                    q_hist_follow_boundary = [q_hist_follow_boundary,q_val];
                    
                elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                    q_val = (1/(x_hist{t_step}(ch_plot_order_x) - x_hist{t_step}(ch_plot_order_x-1)));
                    q_hist_follow_boundary = [q_hist_follow_boundary,q_val];
                    
                else %interior cell
                    q_val = 0.5*( (1/(x_hist{t_step}(ch_plot_order_x+1) - x_hist{t_step}(ch_plot_order_x))) + (1/(x_hist{t_step}(ch_plot_order_x) - x_hist{t_step}(ch_plot_order_x-1))));
                    q_hist_follow_boundary = [q_hist_follow_boundary,q_val];
                end
                
                %cell stiffness
                if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                    k_val = k_hist{t_step}(ch_plot_order_x);
                    k_hist_follow_boundary = [k_hist_follow_boundary,k_val];
                    
                elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                    k_val = k_hist{t_step}(ch_plot_order_x-1);
                    k_hist_follow_boundary = [k_hist_follow_boundary,k_val];
                    
                else %interior cell
                    k_val = 0.5*( k_hist{t_step}(ch_plot_order_x-1) + k_hist{t_step}(ch_plot_order_x));
                    k_hist_follow_boundary = [k_hist_follow_boundary,k_val];
                end
                
                %resting cell length
                if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                    a_val = a_hist{t_step}(ch_plot_order_x);
                    a_hist_follow_boundary = [a_hist_follow_boundary,a_val];
                    
                elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                    a_val = a_hist{t_step}(ch_plot_order_x-1);
                    a_hist_follow_boundary = [a_hist_follow_boundary,a_val];
                    
                else %interior cell
                    a_val = 0.5*( a_hist{t_step}(ch_plot_order_x-1) + a_hist{t_step}(ch_plot_order_x));
                    a_hist_follow_boundary = [a_hist_follow_boundary,a_val];
                end
                        
                %mixing populations id
                
                if iscell(mixing_pop_id_hist) == 0
                    mixing_pop_id_hist_follow_boundary=-1;
                else
                    if m_springs_per_cell == 1
                        if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x-1);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        else
                            mixing_pop_id_val = 0.5*( mixing_pop_id_hist{t_step}(ch_plot_order_x-1) + mixing_pop_id_hist{t_step}(ch_plot_order_x) );
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                        end
                    else
                        if  x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(1) %if at the first position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        elseif x_hist{t_step}(ch_plot_order_x) ==  x_hist{t_step}(end) %at the last position
                            mixing_pop_id_val = mixing_pop_id_hist{t_step}(ch_plot_order_x-1);
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                            
                        else
                            mixing_pop_id_val = 0.5*( mixing_pop_id_hist{t_step}(ch_plot_order_x-1) + mixing_pop_id_hist{t_step}(ch_plot_order_x) );
                            mixing_pop_id_hist_follow_boundary = [mixing_pop_id_hist_follow_boundary,mixing_pop_id_val];
                        end
                    end
                end
                
                
            end
        end
    end
    x_hist_traceback{traceback_Nfin_loop} = x_hist_follow_boundary;
    t_hist_traceback{traceback_Nfin_loop} = t_hist_follow_boundary;
    q_hist_traceback{traceback_Nfin_loop} = q_hist_follow_boundary;
    k_hist_traceback{traceback_Nfin_loop} = k_hist_follow_boundary;
    a_hist_traceback{traceback_Nfin_loop} = a_hist_follow_boundary;
    mixing_pop_id_hist_traceback{traceback_Nfin_loop} = mixing_pop_id_hist_follow_boundary;
end


%% Plot the characteristics

if colouring == 0 %each trajectory has its own colour
    
        figure
        for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
            x_plot =  x_hist_traceback{traceback_Nfin_loop};
            y_plot =  t_hist_traceback{traceback_Nfin_loop};
            scatter(x_plot,y_plot)
            hold on
        end
    
elseif    colouring ==1
    %use the density as the colouring
    %Add in an approximate particle density for colouring
    %determine the colours stepping throught the proliferation times
    
    figure
    for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
        x_plot = x_hist_traceback{traceback_Nfin_loop};
        y_plot =  t_hist_traceback{traceback_Nfin_loop};
        particle_density_plot  = q_hist_traceback{traceback_Nfin_loop};
        if size(x_plot,2) < 21
            scatter(x_plot,y_plot,10, particle_density_plot, 'fill')
        else
            scatter(x_plot([1:1:20,20:10:end]),y_plot([1:1:20,20:10:end]),10, particle_density_plot([1:1:20,20:10:end]), 'fill')
        end
        hold on
    end
    
    title(['Characteristics - Density Colouring - Sim Run ' num2str(sim_run)])
    
elseif colouring==2
    %use the cell stiffness for colouring
   
    figure
    for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
        x_plot = x_hist_traceback{traceback_Nfin_loop};
        y_plot =  t_hist_traceback{traceback_Nfin_loop};
        particle_density_plot  = k_hist_traceback{traceback_Nfin_loop};
        if size(x_plot,2) < 21
            scatter(x_plot,y_plot,10, particle_density_plot, 'fill')
        else
            scatter(x_plot([1:1:20,20:10:end]),y_plot([1:1:20,20:10:end]),10, particle_density_plot([1:1:20,20:10:end]), 'fill')
        end
        hold on
    end
    
    title(['Characteristics - Cell Stiffness Colouring - Sim Run ' num2str(sim_run)])
    
elseif colouring ==3
    %resting cell length
    
    figure
    for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
        x_plot = x_hist_traceback{traceback_Nfin_loop};
        y_plot =  t_hist_traceback{traceback_Nfin_loop};
        particle_density_plot  = a_hist_traceback{traceback_Nfin_loop};
        if size(x_plot,2) < 21
            scatter(x_plot,y_plot,10, particle_density_plot, 'fill')
        else
            scatter(x_plot([1:1:20,20:10:end]),y_plot([1:1:20,20:10:end]),10, particle_density_plot([1:1:20,20:10:end]), 'fill')
        end
        hold on
    end
      
    
elseif colouring  == 4
    %population id colouring
   
    tic
    figure
    hold on
    for traceback_Nfin_loop = 1:size(ch_plot_ord{end},2)
        x_plot =  x_hist_traceback{traceback_Nfin_loop};
        y_plot =  t_hist_traceback{traceback_Nfin_loop};
        
        x_plot_r = x_plot(mixing_pop_id_hist_traceback{traceback_Nfin_loop}==1);
        y_plot_r = y_plot(mixing_pop_id_hist_traceback{traceback_Nfin_loop}==1);
        
        x_plot_b = x_plot(mixing_pop_id_hist_traceback{traceback_Nfin_loop}==0);
        y_plot_b = y_plot(mixing_pop_id_hist_traceback{traceback_Nfin_loop}==0);
        
        x_plot_k = x_plot(mixing_pop_id_hist_traceback{traceback_Nfin_loop}==0.5);
        y_plot_k = y_plot(mixing_pop_id_hist_traceback{traceback_Nfin_loop}==0.5);
        
        if length(x_plot_r) > 0
            if size(x_plot_r,2) < 21
                scatter(x_plot_r,y_plot_r,10, 'b','fill')
            else
                scatter(x_plot_r([1:1:20,20:10:end]),y_plot_r([1:1:20,20:10:end]),10, 'b','fill')
            end
            %scatter(x_plot_r,y_plot_r,'b');
        end
        if length(x_plot_b) > 0
            if size(x_plot,2) < 21
                scatter(x_plot_b,y_plot_b,10,'r','fill')
            else
                scatter(x_plot_b([1:1:20,20:10:end]),y_plot_b([1:1:20,20:10:end]),10, 'r','fill')
            end
            %scatter(x_plot_b,y_plot_b,'r');
        end
        if length(x_plot_k) > 0
            if size(x_plot,2) < 21
                scatter(x_plot_k,y_plot_k,10,'k')
            else
                scatter(x_plot_k([1:1:20,20:10:end]),y_plot_k([1:1:20,20:10:end]),10,'k','fill')
            end
        end
    end
    toc
    
    title(['Characteristics - Population type Colouring - Sim Run ' num2str(sim_run)])
    
end


ylim([0,max([yticks_vec,t_hist])])
%flip the y axis axis
ax = gca;
ax.YDir = 'reverse';

xticks(xticks_vec);
yticks(yticks_vec);

%save the figure with the colorbar
if colouring ~= 4
    colorbar
end

%print(gcf,'-depsc2',  [filepath_save_figs 'Characteristics_with_colourbar_type_' num2str(colouring) '_sim_run_' num2str(sim_run) '2.eps']);
saveas(gcf, [filepath_save_figs 'Characteristics_with_colourbar_type_' num2str(colouring) '_sim_run_' num2str(sim_run) '.fig'])
saveas(gcf, [filepath_save_figs 'Characteristics_with_colourbar_type_' num2str(colouring) '_sim_run_' num2str(sim_run) '.jpg'])



end
