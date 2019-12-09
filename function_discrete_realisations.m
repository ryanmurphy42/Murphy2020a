function function_discrete_realisations(file_save_name, num_sims, k1, k2, a1, a2, prolif_law, ...
    beta_prolif_1, beta_prolif_2, death_law,gamma_death_1, gamma_death_2, characteristics_include, timestop, time_record_vector, N_init,...
    gaussian_n,gaussian_sigma, gaussian_mu, L,mixing_pop_id ,eta,dt,l_crit_death,q_init_condition,folder_name)

%% Initialise arrays for storage
t_hist_sim_run={};
x_hist_sim_run={};
k_hist_sim_run={};
a_hist_sim_run={};
G_hist_sim_run={};
cell_with_proliferation_event_hist_sim_run={};
cell_with_death_event_hist_sim_run = {};
cell_event_hist_sim_run = {};
prolif_death_indicator_hist_sim_run = {};
cell_death_lr_hist_sim_run= {};
mixing_pop_id_hist_sim_run={};

%% Initialise variables
N=N_init;
rng('default')

%% Cell stiffness, k, properties

k=zeros(N,1);
for ii =1:N
    if mixing_pop_id(ii) == 0
        k(ii) = k1;
    elseif mixing_pop_id(ii) == 1
        k(ii) = k2;
    end
end

%% Resting cell length, a, properties

a=zeros(N,1);
for ii =1:N
    if mixing_pop_id(ii) == 0
        a(ii) = a1;
    elseif mixing_pop_id(ii) == 1
        a(ii) = a2;
    end
end

%% beta prolif and beta death
beta_prolif=zeros(N,1);
gamma_death=zeros(N,1);
for ii =1:N
    if mixing_pop_id(ii) == 0
        beta_prolif(ii) = beta_prolif_1;
        gamma_death(ii) = gamma_death_1;
    elseif mixing_pop_id(ii) == 1
        beta_prolif(ii) = beta_prolif_2;
        gamma_death(ii) = gamma_death_2;
    end
end


%Initial density condition
x0 = linspace(0,L,N+1);
if q_init_condition == 1
    %Uniform initial condition
    x_current = x_init;
elseif q_init_condition ==2
    %Gaussian initial condition
    x_current= fsolve(@(x) function_p1d1_fsolve_initial_positions_gaussian_v1(N,L,gaussian_n,gaussian_sigma,gaussian_mu,x ),x0)';
elseif q_init_condition == 3
    %Mechanical relaxation initial condition
    x_current = function_density_mechrelax_initialcondition( N_init, L, k, a );
end


for sim_run=1:num_sims
    
    %% Discrete numerical solver for one realisation
    
    [t_hist, x_hist_m, k_hist_m, a_hist_m, G_hist_m, cell_with_proliferation_event_hist,cell_with_death_event_hist,cell_event_hist, prolif_death_indicator_hist, cell_death_lr_hist,mixing_pop_id_m_hist] = ...
        function_discrete_onerealisation(L, N, k, a, eta, beta_prolif, gamma_death, x_current, prolif_law,death_law, ...
        timestop, time_record_vector, mixing_pop_id,characteristics_include,dt,l_crit_death);
    
    %% Store variables
    
    t_hist_sim_run{sim_run} = t_hist;
    x_hist_sim_run{sim_run} = x_hist_m;
    k_hist_sim_run{sim_run} = k_hist_m;
    a_hist_sim_run{sim_run} = a_hist_m;
    G_hist_sim_run{sim_run} = G_hist_m;
    cell_with_proliferation_event_hist_sim_run{sim_run} = cell_with_proliferation_event_hist;
    cell_with_death_event_hist_sim_run{sim_run} = cell_with_death_event_hist;
    cell_event_hist_sim_run{sim_run} = cell_event_hist;
    prolif_death_indicator_hist_sim_run{sim_run} = prolif_death_indicator_hist;
    cell_death_lr_hist_sim_run{sim_run} = cell_death_lr_hist;
    mixing_pop_id_hist_sim_run{sim_run} = mixing_pop_id_m_hist;
    
    
end


%Save .mat file with outputs
save([pwd '\' folder_name '\' file_save_name],'-v7.3',...
    't_hist_sim_run',...
    'x_hist_sim_run',...
    'k_hist_sim_run',...
    'a_hist_sim_run',...
    'eta',...
    'L',...
    'num_sims',...
    'mixing_pop_id_hist_sim_run',...
    'cell_with_proliferation_event_hist_sim_run',...
    'cell_with_death_event_hist_sim_run',...
    'cell_event_hist_sim_run',...
    'prolif_death_indicator_hist_sim_run',...
    'cell_death_lr_hist_sim_run',...
    'G_hist_sim_run',...
    'timestop')
