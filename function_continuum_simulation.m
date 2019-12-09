function function_continuum_simulation(file_save_name,q_init_condition, ...
    k1,k2,a1,a2, beta_prolif_1, beta_prolif_2, gamma_death_1, gamma_death_2, l_crit_death, N1, N2, ...
    L,timestop,time_vector_record, folder_name,max_iters,err_tol,dx,dt,gaussian_n,...
    gaussian_sigma,gaussian_mu ,eta,k_hom,a_hom,beta_hom,gamma_hom,one_or_two_pop)

%Run continuum simulation

%% Initial conditions

%Number of nodes
nodesx = round(L/dx)+1.00;

%Spatial discretisation
x_discretisation = 0:dx:L;

% Density initial condition
if q_init_condition == 1
    %Uniform initial condition
    q0 = zeros(nodesx,1);
    q1=L/(N1+N2);
    for i=1:nodesx
        q0(i) = q1;
    end
elseif q_init_condition == 2
    %Gaussian initial condition
    function_q_init = @(x) ((gaussian_n)/(sqrt(2*pi*gaussian_sigma^2)))*(exp(-((x-gaussian_mu)^2)/(2*gaussian_sigma^2)));
    q0 = zeros(nodesx,1);
    for i=1:nodesx
        q0(i) = function_q_init((i-1)*dx);
    end
elseif  q_init_condition == 3
    %Mechanical relaxation initial condition - two populations
    interface_position = ((k1*a1/k2) + (L/N2) - a2)/( (k1/(k2*N1)) + (1/N2)   );
    q1 = N1/interface_position;
    q2 = N2/(L-interface_position);
    
    q0 = zeros(nodesx,1);
    for i=1:nodesx
        if (i-1)*dx <= interface_position
            q0(i) = q1;
        else
            q0(i) = q2;
        end
    end
end

%Interface position initial condition - if only one population dont have an interface.
if one_or_two_pop == 1
    interface_position = L+1;
end
s_prev = interface_position;


%% cell properties initial conditions
k = zeros(nodesx,1);
for i=1:nodesx
    if (i-1)*dx <= interface_position
        k(i) = k1;
    else
        k(i) = k2;
    end
end


a=zeros(nodesx,1);
for i=1:nodesx
    if (i-1)*dx <= interface_position
        a(i) = a1;
    else
        a(i) = a2;
    end
end


beta_prolif = zeros(nodesx,1);
for i=1:nodesx
    if (i-1)*dx <= interface_position
        beta_prolif(i) = beta_prolif_1;
    else
        beta_prolif(i) = beta_prolif_2;
    end
end


gamma_death = zeros(nodesx,1);
for i=1:nodesx
    if (i-1)*dx <= interface_position
        gamma_death(i) = gamma_death_1;
    else
        gamma_death(i) = gamma_death_2;
    end
end


%% Store the intial conditions

loop_count_stored = 1;
max_recorded_data_points = 10000; %value to initialise for below


q_hist = zeros(nodesx,max_recorded_data_points);
q_hist(:,loop_count_stored) = q0;

k_hist = zeros(nodesx, max_recorded_data_points);
k_hist(:,loop_count_stored) = k;

a_hist = zeros(nodesx, max_recorded_data_points);
a_hist(:,loop_count_stored) = a;

beta_prolif_hist = zeros(nodesx, max_recorded_data_points);
beta_prolif_hist(:,loop_count_stored) = beta_prolif;

gamma_death_hist = zeros(nodesx, max_recorded_data_points);
gamma_death_hist(:,loop_count_stored) = gamma_death;

node_velocity_hist=zeros(nodesx,max_recorded_data_points);

t_hist  = zeros(1,max_recorded_data_points);
s_hist  = zeros(1, max_recorded_data_points);
t_hist(1,loop_count_stored)=0;
s_hist(1,loop_count_stored)=s_prev;

% Store the variables to use in the temporal loop
q_prev =q0;
k_prev =k;
a_prev =a;
beta_prolif_prev = beta_prolif;
gamma_death_prev = gamma_death;

t=0; %initial time

bc_switch_made=0; %for two populations, 0-start assuming two populations.
%Changes to 1 when two populations reduce to 1.

%% Run the temporal loop.

while t < timestop
   
   
    t = t + dt; % Update time
    
    min_diff_to_time_record = min(abs(time_vector_record -t)); %how close is current time to times where we wish to record results
    
    if bc_switch_made ==0
        if one_or_two_pop == 2
            % If the number of cells in the right tissue is less than 1 then remove the tissue
            [~,snode] = min(abs(x_discretisation - s_prev));
            
            cellspop2  = trapz(((snode-1)*dx):dx:((nodesx-1)*dx),q_prev((snode:end)));
            
            if cellspop2 < 1
                k_prev = k1*ones(nodesx,1);
                a_prev = a1*ones(nodesx,1);
                beta_prolif_prev = beta_prolif_1*ones(nodesx,1);
                gamma_death_prev = gamma_death_1*ones(nodesx,1);
                bc_switch_made = 1;
                s_prev=((nodesx-1)*dx);
            end
            
            % If the number of cells in the left tissue is less than 1 then remove the tissue
            cellspop1  = trapz(0:dx:((snode-1)*dx),q_prev((1:snode)));
            if cellspop1 < 1
                k_prev = k2*ones(nodesx,1);
                a_prev = a2*ones(nodesx,1);
                beta_prolif_prev = beta_prolif_2*ones(nodesx,1);
                gamma_death_prev = gamma_death_2*ones(nodesx,1);
                bc_switch_made = 1;
                s_prev=0;
            end
        end
    end
    
    % Newton-Raphson method
    [q,k,a,beta_prolif,gamma_death,s,node_velocity,bc_switch_made,~,~] = function_continuum_newton_solver_linearpd(q_prev,k_prev,a_prev,beta_prolif_prev,gamma_death_prev,l_crit_death,s_prev,...
        nodesx,dt,dx,eta,max_iters,err_tol,k_hom, a_hom, beta_hom, gamma_hom,one_or_two_pop,bc_switch_made);
    
    
    %Update previous
    q_prev = q;
    k_prev = k;
    a_prev = a;
    beta_prolif_prev = beta_prolif;
    gamma_death_prev = gamma_death;
    s_prev = s;
    % Record results
    if  (min_diff_to_time_record < dt*10)
        loop_count_stored = loop_count_stored + 1;
        q_hist(:,loop_count_stored) = q;
        k_hist(:,loop_count_stored) = k;
        a_hist(:,loop_count_stored) = a;
        beta_prolif_hist(:,loop_count_stored) = beta_prolif;
        gamma_death_hist(:,loop_count_stored) = gamma_death;
        node_velocity_hist(:,loop_count_stored) = node_velocity;
        t_hist(loop_count_stored) = t;
        s_hist(loop_count_stored) = s;
    end
    
    %% Exit clauses
    if ( sum(q > 1000) > 1 )
        disp('q greater than 1000')
        q
        break
    elseif sum( isnan(q)) >0
        t
        disp('is nan q is true')
        break
    elseif max(abs(q)) > 1000
        t
        disp('max abs q > 1000')
        break
    elseif (t > 1.01*max(time_vector_record))
        t
        disp('1.01 times max time_vector')
        break
    end
    
    
end

%Save the data.
save([pwd '\' folder_name '\' file_save_name],'-v7.3');

