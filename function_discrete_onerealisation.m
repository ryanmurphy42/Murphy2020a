function [t_hist, x_hist_m, k_hist_m, a_hist_m, G_hist_m, cell_with_proliferation_event_hist,cell_with_death_event_hist,cell_event_hist, prolif_death_indicator_hist, cell_death_lr_hist,mixing_pop_id_m_hist] = ...
    function_discrete_onerealisation(L, N, k, a, eta, beta_prolif, gamma_death, x_current, prolif_law,death_law, ...
    timestop, time_record_vector, mixing_pop_id,characteristcs_include,dt,l_crit_death)


%% Initialise for storage

t_hist = [];
x_hist_m = {};
k_hist_m = {};
a_hist_m = {};
G_hist_m = {};
prolif_times=[];
death_times=[];
cell_with_proliferation_event_hist = [];
cell_with_death_event_hist = [];
cell_event_hist = [];
prolif_death_indicator_hist =[];
cell_death_lr_hist = [];

%% Initialise variables
t=0;
t_record_loop=1;

%% convert from cell properties to spring properties
N_cells=N;
m_springs_per_cell=1;
N_springs = N_cells*m_springs_per_cell;

%convert the positions to springs per cell
x_current_m=zeros(N_cells*m_springs_per_cell + 1,1);
k_m=zeros(N_cells*m_springs_per_cell,1);
a_m=zeros(N_cells*m_springs_per_cell,1);
beta_prolif_m=zeros(N_cells*m_springs_per_cell,1);
gamma_death_m=zeros(N_cells*m_springs_per_cell,1);
mixing_pop_id_m=zeros(N_cells*m_springs_per_cell,1);

%convert to spring properties
cell_spring_index = 0;
for ii=1:N_cells
    for jj=1:m_springs_per_cell
        cell_spring_index = cell_spring_index +1;
        x_current_m(cell_spring_index) = x_current(ii) + (x_current(ii+1) -x_current(ii))*((jj-1)/m_springs_per_cell);
        k_m(cell_spring_index) = k(ii)*m_springs_per_cell;
        a_m(cell_spring_index) = a(ii)/m_springs_per_cell;
        beta_prolif_m(cell_spring_index) = beta_prolif(ii)/m_springs_per_cell;
        gamma_death_m(cell_spring_index) = gamma_death(ii)/m_springs_per_cell;
        mixing_pop_id_m(cell_spring_index) = mixing_pop_id(ii);
    end
end
x_current_m(N_cells*m_springs_per_cell + 1) = L;
eta_m = eta/m_springs_per_cell;

%% calculate the proliferation and death propensities

[G_m,~,~] = function_proliferation_death_rates(N_cells,m_springs_per_cell, prolif_law,death_law,beta_prolif_m, gamma_death_m,x_current_m,l_crit_death);

%% save the spring properties
t_hist=[t_hist,t];
x_hist_m{t_record_loop} = x_current_m;
k_hist_m{t_record_loop} = k_m;
a_hist_m{t_record_loop} = a_m;
G_hist_m{t_record_loop} = G_m;
mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;



while (t<timestop)
    
    tau = timestop;
    
    if tau <= timestop + 1
        
        %Run the numerical simulation
        timesteps =1;
        
        if N_cells == 0
            
            t_record_loop=t_record_loop+1;
            t_hist=[t_hist,t];
            x_hist_m{t_record_loop} = x_current_m;
            k_hist_m{t_record_loop} = k_m/m_springs_per_cell;
            a_hist_m{t_record_loop} = a_m*m_springs_per_cell;
            G_hist_m{t_record_loop} = G_m*m_springs_per_cell;
            mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;
            t=timestop+1;
            
        else
            
            for timestep_loop = 1:timesteps
                
                %% Update position
                x_prev_m = x_current_m;
                for i=2:N_springs
                    x_current_m(i) = x_prev_m(i) + dt*(1/eta_m)*(    k_m(i-1)*(a_m(i-1) - (abs(x_prev_m(i) - x_prev_m(i-1))))*((x_prev_m(i) - x_prev_m(i-1))/(abs(x_prev_m(i) - x_prev_m(i-1)))) ...
                        + (k_m(i)*(a_m(i) - (abs(x_prev_m(i) - x_prev_m(i+1)) )))*((x_prev_m(i) - x_prev_m(i+1))/(abs(x_prev_m(i) - x_prev_m(i+1)))) );
                end
                
                
                t=t+dt;
                
                %calculate the absolute difference
                t_cum_sum_inf = abs(time_record_vector -t);
                
                %record data if close to timepoint
                if  min(t_cum_sum_inf) < 2.5*dt
                    t_record_loop=t_record_loop+1;
                    t_hist=[t_hist,t];
                    x_hist_m{t_record_loop} = x_current_m;
                    k_hist_m{t_record_loop} = k_m/m_springs_per_cell;
                    a_hist_m{t_record_loop} = a_m*m_springs_per_cell;
                    G_hist_m{t_record_loop} = G_m*m_springs_per_cell;
                    mixing_pop_id_m_hist{t_record_loop} = mixing_pop_id_m;
                end
                
            end
        end
        
        x_new_spring_positions = 0:(1/(2*m_springs_per_cell)):1;
        x_new_spring_positions_trunc = x_new_spring_positions(2:end-1)';
        
        %% update the propensities
        
        [G_m,~,~] = function_proliferation_death_rates(N_cells,m_springs_per_cell, prolif_law,death_law,beta_prolif_m, gamma_death_m,x_current_m,l_crit_death);
        
        %% calculate if there is a cell event
        r1 = rand(1);
        G_m_cumsum = cumsum(G_m)/sum(G_m);
        
        cell_event_happening =0;
        if  r1 < sum(G_m*dt)
            cell_event_happening =1;
            
            %% save data
            if characteristcs_include==1
                time_record_vector = sort([time_record_vector,t]);
            end
            
        end
        
        if  cell_event_happening ==1
            
            r2 = rand(1);
            %% Determine whether have cell proliferation or cell death
            
            G_m_cum_sum_inf = G_m_cumsum(1:end) -r2;
            
            G_m_cum_sum_inf(G_m_cum_sum_inf < 0) = inf;
            
            [~ , spring_event]   = min(G_m_cum_sum_inf);
            
            %determine cell event
            
            cell_event = floor((spring_event-1)/m_springs_per_cell)+ 1;
            
            if spring_event <=  N_springs %cell proliferation event
                
                if cell_event == 1
                    
                    %Update the positions and enforce the cell proliferation event by dividing the cell equally
                    %Relabel the x positions
                    
                    x_current_m = [x_current_m(1);
                        x_current_m(1) + x_new_spring_positions_trunc*(x_current_m(cell_event*m_springs_per_cell+1)-x_current_m(1));
                        x_current_m(cell_event*m_springs_per_cell+1:end)];
                    
                    
                elseif cell_event == N_cells
                    
                    %Update the positions and enforce the cell proliferation event by dividing the cell equally
                    %Relabel the x positions
                    
                    x_current_m = [x_current_m(1:(cell_event-1)*m_springs_per_cell +1);
                        x_current_m((cell_event-1)*m_springs_per_cell +1) + x_new_spring_positions_trunc*(x_current_m(end)-x_current_m((cell_event-1)*m_springs_per_cell +1)); ...
                        x_current_m(end)];
                    
                else
                    
                    %Update the positions and enforce the cell proliferation event by dividing the cell equally
                    %Relabel the x positions
                    
                    x_current_m = [x_current_m(1:(cell_event-1)*m_springs_per_cell +1);
                        x_current_m((cell_event-1)*m_springs_per_cell + 1) + x_new_spring_positions_trunc*(x_current_m(cell_event*m_springs_per_cell+1)-x_current_m((cell_event-1)*m_springs_per_cell+ 1));
                        x_current_m(cell_event*m_springs_per_cell+1:end)];
                    
                end
                
                %Update k
                k_m = [k_m(1:cell_event*m_springs_per_cell); k_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ;k_m(cell_event*m_springs_per_cell+1:end)];
                
                %Update a
                a_m = [a_m(1:cell_event*m_springs_per_cell); a_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ;a_m(cell_event*m_springs_per_cell+1:end)];
                
                %Update beta
                beta_prolif_m = [beta_prolif_m(1:cell_event*m_springs_per_cell); beta_prolif_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ...
                    ;beta_prolif_m(cell_event*m_springs_per_cell+1:end)];
                
                gamma_death_m = [gamma_death_m(1:cell_event*m_springs_per_cell); gamma_death_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ...
                    ;gamma_death_m(cell_event*m_springs_per_cell+1:end)];
                
                %Update mixing_pop_id
                mixing_pop_id_m = [mixing_pop_id_m(1:cell_event*m_springs_per_cell); mixing_pop_id_m(cell_event*m_springs_per_cell)*ones(m_springs_per_cell,1) ;mixing_pop_id_m((cell_event)*m_springs_per_cell + 1:end)];
                
                
                %Update the number of cells
                N_cells = N_cells + 1;
                
                N_springs = N_springs + m_springs_per_cell;
                
                cell_with_proliferation_event_hist = [cell_with_proliferation_event_hist,cell_event];
                
                cell_event_hist = [cell_event_hist,cell_event];
                
                prolif_death_indicator_hist = [prolif_death_indicator_hist, 1]; %1 correspond cell proliferation
                
                prolif_times = [prolif_times, t];
                
            else %cell death event
                %Update the positions and enforce the cell proliferation event by dividing the cell equally
                %Relabel the x positions
                cell_event_death = cell_event - N_cells;
                
                death_times = [death_times, t];
                
                x_new_spring_positions_death_coal = 0:(1/(m_springs_per_cell)):1;
                x_new_spring_positions_death_coal_trunc = x_new_spring_positions_death_coal(2:end-1)';
                
                if cell_event_death == 1 % first cell
                    
                    if m_springs_per_cell == 1
                        x_current_m = [x_current_m(cell_event_death);...
                            x_current_m(((cell_event_death+1)*m_springs_per_cell + 1):end)];
                    else
                        x_current_m = [x_current_m(cell_event_death);
                            x_current_m(cell_event_death) + x_new_spring_positions_death_coal_trunc*(x_current_m((cell_event_death+1)*m_springs_per_cell+1)-x_current_m(cell_event_death));
                            x_current_m(((cell_event_death+1)*m_springs_per_cell+1):end)];
                    end
                    
                    %Update k
                    k_m = [k_m(cell_event_death*m_springs_per_cell +1:end)];
                    
                    
                    %Update a
                    a_m = [a_m(cell_event_death*m_springs_per_cell +1:end)];
                    
                    %Update beta prolif
                    beta_prolif_m = [beta_prolif_m(cell_event_death*m_springs_per_cell +1:end)];
                    
                    %Update beta death
                    gamma_death_m = [gamma_death_m(cell_event_death*m_springs_per_cell +1:end)];
                    
                    %Update mixing_pop_id
                    mixing_pop_id_m = [mixing_pop_id_m(cell_event_death*m_springs_per_cell +1:end)];
                    
                    cell_death_lr_hist = [cell_death_lr_hist, 2 ];
                    
                    
                elseif cell_event_death == N_cells %last cell
                    
                    
                    if m_springs_per_cell == 1
                        x_current_m = [x_current_m(1:cell_event_death-1);
                            x_current_m(end)];
                    else
                        x_current_m = [x_current_m(1:((cell_event_death-2)*m_springs_per_cell+1));
                            x_current_m(((cell_event_death-2)*m_springs_per_cell+1)) + x_new_spring_positions_death_coal_trunc*(x_current_m(end)-x_current_m(((cell_event_death-2)*m_springs_per_cell+1)));
                            x_current_m(end)];
                    end
                    
                    %Update k
                    k_m = [k_m(1:(cell_event_death-1)*m_springs_per_cell) ;k_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    %Update a
                    a_m = [a_m(1:(cell_event_death-1)*m_springs_per_cell) ;a_m((cell_event_death)*m_springs_per_cell +1:end)];
                    
                    %Update beta
                    beta_prolif_m = [beta_prolif_m(1:(cell_event_death-1)*m_springs_per_cell) ;beta_prolif_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    gamma_death_m = [gamma_death_m(1:(cell_event_death-1)*m_springs_per_cell) ;gamma_death_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    %Update mixing_pop_id
                    mixing_pop_id_m = [mixing_pop_id_m(1:(cell_event_death-1)*m_springs_per_cell) ;mixing_pop_id_m((cell_event_death)*m_springs_per_cell +1:end)];
                    
                    cell_death_lr_hist = [cell_death_lr_hist, 1 ];
                    
                else
                    x_new_spring_positions_death_coal = 0:(1/(m_springs_per_cell)):1;
                    x_new_spring_positions_death_coal_trunc = x_new_spring_positions_death_coal(2:end-1)';
                    %% coalesce the two boundaries
                    
                    if m_springs_per_cell == 1
                        x_current_m = [x_current_m(1:cell_event_death-1);
                            (x_current_m(cell_event_death-1)+x_current_m(cell_event_death+2))/2;
                            x_current_m(cell_event_death+2:end)];
                    else
                        x_mid_death = (x_current_m(((cell_event_death-1)*m_springs_per_cell) +1) + x_current_m((cell_event_death)*m_springs_per_cell +1))/2;
                        
                        x_current_m = [x_current_m(1:(cell_event_death-2)*m_springs_per_cell +1);
                            x_current_m((cell_event_death-2)*m_springs_per_cell +1) + x_new_spring_positions_death_coal_trunc*(x_mid_death - ( x_current_m((cell_event_death-2)*m_springs_per_cell +1)        ));
                            x_mid_death;
                            x_mid_death + x_new_spring_positions_death_coal_trunc*(x_current_m(((cell_event_death+1)*m_springs_per_cell + 1))-x_mid_death);
                            x_current_m(((cell_event_death+1)*m_springs_per_cell + 1):end)];
                    end
                    
                    %Update k
                    k_m = [k_m(1:(cell_event_death-1)*m_springs_per_cell) ;k_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    %Update a
                    a_m = [a_m(1:(cell_event_death-1)*m_springs_per_cell) ;a_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    %Update beta
                    beta_prolif_m = [beta_prolif_m(1:(cell_event_death-1)*m_springs_per_cell) ;beta_prolif_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    gamma_death_m = [gamma_death_m(1:(cell_event_death-1)*m_springs_per_cell) ;gamma_death_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    %Update mixing_pop_id
                    mixing_pop_id_m = [mixing_pop_id_m(1:(cell_event_death-1)*m_springs_per_cell) ;mixing_pop_id_m((cell_event_death)*m_springs_per_cell + 1:end)];
                    
                    cell_death_lr_hist = [cell_death_lr_hist, 1 ];
                    
                    
                end
                
                %Update the number of cells
                N_cells = N_cells - 1;
                
                N_springs = N_springs - m_springs_per_cell;
                
                cell_with_death_event_hist = [cell_with_death_event_hist,cell_event_death];
                
                prolif_death_indicator_hist = [prolif_death_indicator_hist, 2]; %2 correspond to cell death
                
                cell_event_hist = [cell_event_hist,cell_event_death];
                
            end
        end
    else
        
    end
    
end

end


