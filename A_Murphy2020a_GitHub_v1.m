%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Murphy2020a - Mechanical cell competition in heterogeneous epithelial tissues.

% Key algorithms used to generate paper Figures 4 (with proliferation and death) and 6.
% Figure 4 - Homogeneous population - linear proliferation and death mechanism - shown for proliferation and death included
% Figure 6 - Heterogeneous population - linear proliferation and death mechanism

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Structure of code

%% 1) Figure 4c
% 1.1) Generate discrete characteristics
% 1.2) Generate discrete realizations
% 1.3) Generate continuum model results
% 1.4) Discrete plots - characteristics
% 1.5) Plot density snapshots and total cell number discrete
% 1.6) Overlay continuum density snapshots and total cell number over discrete

%% 2) Figure 6
% 2.1) Generate discrete characteristics
% 2.2) Generate discrete realizations
% 2.3) Generate continuum model results
% 2.4) Plot characteristic
% 2.5) Plot density snapshots, total cell number discrete, and interface position
% 2.6) Overlay continuum density snapshots, total cell number discrete, and interface position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Figure 4c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1) Generate discrete characteristics

simulation_id = 'Fig4c';

num_sims=1; %number of discrete realisations
k1=10; % cell stiffness population 1
k2=10; % cell stiffness population 2
a1=0; % resting cell length population 1
a2=0; % resting cell length population 2
prolif_law=202; %201 - constant, 202 - linear, 203 - logistic
beta_prolif_1=0.07; % beta population 1
beta_prolif_2=0.07; % beta population 2
death_law=202; %201 - constant, 202 - linear, 203 - logistic
gamma_death_1=0.35; % gamma population 1
gamma_death_2=0.35; % gamma population 2
characteristics_include=1; %1-include every proliferation event to reproduce characteristics, 0-fewer time points recorded
timestop=100; %final time
time_record_vector=0:0.1:timestop; %time points to record
gaussian_n=4*11.056818502502322; %scale factor if gaussian initial density condition
gaussian_sigma=sqrt(9); %standard deviation if gaussian initial density condition
gaussian_mu=5; %mean if gaussian initial condition
N_init=40; % initial number of cells
L=10; % length of the domain
mixing_pop_id = [0*ones(N_init,1)]; %mixing - here all zero so only one population
eta=1; %viscosity coefficient
dt = 0.0001; %timestep
l_crit_death =0.3; %if linear proliferation death law, the first value of length where death is zero
file_save_name = ['Results_Sim' simulation_id '_Discrete_Characteristics'];
q_init_condition=2; %1- uniform, 2-gaussia, 3- mechanical relaxation

%Create folder
folder_name = [simulation_id '_Discrete'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end

%Generate a .mat file for discrete realisations
function_discrete_realisations(file_save_name, num_sims, k1, k2, a1, a2, prolif_law, ...
    beta_prolif_1, beta_prolif_2, death_law, gamma_death_1, gamma_death_2, characteristics_include, timestop, time_record_vector, N_init,...
    gaussian_n,gaussian_sigma, gaussian_mu, L,mixing_pop_id ,eta,dt,l_crit_death,q_init_condition,folder_name)



%% 1.2) Generate discrete realizations
% Smaller number of time points recorded as do not need to reproduce characteristics with this data.

%Parameters which are different from 1.1)
num_sims=10;
characteristics_include=0; %0-recording fewer timesteps as not reproducing characteristics here
time_record_vector=0:5:timestop;

file_save_name = ['Results_Sim' simulation_id '_Discrete_R' num2str(num_sims)];


function_discrete_realisations(file_save_name, num_sims, k1, k2, a1, a2, prolif_law, ...
    beta_prolif_1, beta_prolif_2, death_law, gamma_death_1, gamma_death_2, characteristics_include, timestop, time_record_vector, N_init,...
    gaussian_n,gaussian_sigma, gaussian_mu, L,mixing_pop_id ,eta,dt,l_crit_death,q_init_condition,folder_name)


%% 1.3) Generate continuum model results
clear
clc

%This is for the linear proliferation and death mechanism %%
simulation_id = 'Fig4c';

q_init_condition = 2; %1-uniform, 2-gaussian, 3-mechanical relaxation
k1=10;
k2=10;
a1=0;
a2=0;
beta_prolif_1 = 0.07;
beta_prolif_2 = 0.07;
gamma_death_1 =0.35;
gamma_death_2 =0.35;
N1=20;
N2=20;
l_crit_death=0.3;
L=10;
gaussian_n=4*11.056818502502322;
gaussian_sigma=sqrt(9);
gaussian_mu=5;
dx=0.01; %continuum spatial step
dt =0.0001; %continum time step
eta=1;
k_hom=0; %0 - k heterogeneous, 1- k homogeneous (this is to improve speed for homogeneous properties)
a_hom=1; %0 - a heterogeneous, 1- a homogeneous (this is to improve speed for homogeneous properties)
beta_hom=1; %0 - beta heterogeneous, 1 - beta homogeneous (this is to improve speed for homogeneous properties)
gamma_hom=1; %0 - gamma heterogeneous, 1 - gamma homogeneous (this is to improve speed for homogeneous properties)
max_iters = 1000; %Newton-Raphson - maximum number of iterations per timestep
err_tol = 0.0001; %Newton-Rapshon - error tolerance
timestop=100;
time_vector_record=0:5:timestop;

file_save_name = ['Results_Sim' simulation_id '_Continuum'];

folder_name = [simulation_id '_Continuum'];
if ~exist(folder_name', 'dir')
    mkdir(folder_name);
end

one_or_two_pop=1; %1 - one population, 2 - two populations

tic
function_continuum_simulation(file_save_name,q_init_condition, ...
    k1,k2,a1,a2, beta_prolif_1, beta_prolif_2, gamma_death_1, gamma_death_2, l_crit_death, N1, N2, ...
    L,timestop,time_vector_record, folder_name,max_iters,err_tol,dx,dt,gaussian_n,gaussian_sigma,gaussian_mu ,eta,...
    k_hom,a_hom,beta_hom,gamma_hom,one_or_two_pop);
toc


%% 1.4) Discrete plots - characteristics
clear
clc

simulation_id = 'Fig4c';

filepath_save_figs = [pwd '\' simulation_id '_Discrete\'];

load([filepath_save_figs 'Results_Sim' simulation_id '_Discrete_Characteristics.mat']);

%check that the all the data is available to plot characteristics without error. If the below returns a result decreate dt.
m_springs_per_cell=1;% 1 spring per cell is standard - code should be amended for m>1 results (not included)
function_characteristic_plot_check(num_sims,t_hist_sim_run,x_hist_sim_run,m_springs_per_cell);


%define figure properties
xticks_vec=linspace(0,L,5); %figure x ticks
yticks_vec=[0,25,50,75,100]; %figure y ticks

%plot the characteristics
for sim_run = 1:num_sims
    
    %for each simulation extract the data
    cell_with_proliferation_event_hist=cell_with_proliferation_event_hist_sim_run{sim_run};
    cell_with_death_event_hist=cell_with_death_event_hist_sim_run{sim_run};
    cell_event_hist=cell_event_hist_sim_run{sim_run};
    prolif_death_indicator_hist=prolif_death_indicator_hist_sim_run{sim_run};
    cell_death_lr_hist = cell_death_lr_hist_sim_run{sim_run};
    
    t_hist=t_hist_sim_run{sim_run};
    x_hist=x_hist_sim_run{sim_run};
    k_hist =k_hist_sim_run{sim_run};
    a_hist =a_hist_sim_run{sim_run};
    G_hist =G_hist_sim_run{sim_run};
    mixing_pop_id_hist =mixing_pop_id_hist_sim_run{sim_run};
    
    %plot with density colouring
    colouring =1; %1- density, 2- cell stiffness, 3- resting cell length, 4-population id colouring
    function_characteristics_plot(t_hist, x_hist, k_hist,G_hist, mixing_pop_id_hist, a_hist, L,...
        filepath_save_figs, colouring, xticks_vec, yticks_vec,....
        cell_with_proliferation_event_hist,cell_with_death_event_hist,cell_event_hist,prolif_death_indicator_hist,sim_run,cell_death_lr_hist,m_springs_per_cell)
    
    close all
end


%% 1.5) Plot density snapshots and total cell number discrete
clear
clc

simulation_id = 'Fig4c';
num_sims=2000;

filepath_save_figs = [pwd '\' simulation_id '_Discrete\'];

load([filepath_save_figs 'Results_Sim' simulation_id '_Discrete_R' num2str(num_sims) '.mat']);

% 1.5.1) Density snapshots

density_evaluation_x_dx = 0.01; %spatial step to average the discrete realisations over
plot_times = 0:5:timestop;


function_discrete_density_plot(t_hist_sim_run, x_hist_sim_run, L, ...
    filepath_save_figs,num_sims, plot_times, density_evaluation_x_dx)

% 1.5.2) Total cell number

one_or_two_pop=1;
function_discrete_total_cell_number(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims,...
    plot_times,mixing_pop_id_hist_sim_run,one_or_two_pop)



%% 1.6) Overlay continuum density snapshots and total cell number over discrete
clear
clc

simulation_id = 'Fig4c';
filepath_save_figs = [pwd '\' simulation_id '_Continuum\'];
load([filepath_save_figs 'Results_Sim' simulation_id '_Continuum.mat']);

filepath_load_figs = [pwd '\' simulation_id '_Discrete\'];

% 1.6.1) Density snapshots

colouring=1;
ctm_only=0;
plot_times=0:5:timestop; %must match the save times for the discrete simulations (for accuracy)

function_discrete_ctm_comparison_snapshots_plot(t_hist, q_hist,k_hist,...
    L,dx, filepath_save_figs,  colouring,ctm_only, plot_times,filepath_load_figs)


% 1.6.2) Total cell number
ctm_only=0; %0 - discrete and continuum, 1- continuum only.
plot_times=0:5:timestop;
function_dis_ctm_total_cell_number_plot(t_hist, q_hist,s_hist,L,dx, filepath_load_figs,filepath_save_figs, plot_times, ctm_only,one_or_two_pop)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reset matlab

clear
clc

%% 2.1) Generate discrete characteristics


simulation_id = 'Fig6';

num_sims=1;
k1=10;
k2=20;
a1=0;
a2=0;
prolif_law=202;
beta_prolif_1=0.07;
beta_prolif_2=0.07;
death_law=202;
gamma_death_1=0.35;
gamma_death_2=0.35;
characteristics_include=1;
timestop=200;
time_record_vector=0:0.5:timestop;
gaussian_n=4*11.056818502502322;
gaussian_sigma=sqrt(9);
gaussian_mu=5;
N_init=40;
L=10;
mixing_pop_id = [0*ones(N_init/2,1); 1*ones(N_init/2,1)];
eta=1;
dt = 0.0001;
l_crit_death =0.3;
file_save_name = ['Results_Sim' simulation_id '_Discrete_Characteristics'];
q_init_condition=3; %start at steady state

folder_name = [simulation_id '_Discrete'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end


function_discrete_realisations(file_save_name, num_sims, k1, k2, a1, a2, prolif_law, ...
    beta_prolif_1, beta_prolif_2, death_law, gamma_death_1, gamma_death_2, characteristics_include, timestop, time_record_vector, N_init,...
    gaussian_n,gaussian_sigma, gaussian_mu, L,mixing_pop_id ,eta,dt,l_crit_death,q_init_condition,folder_name)


%% 2.2) Generate discrete realizations

%Smaller number of time points recorded as do not need to reproduce characteristics with this data.

%Same parameters as in 2.1 with the following exceptions
num_sims=10;
characteristics_include=0;
time_record_vector=0:5:timestop;

file_save_name = ['Results_Sim' simulation_id '_Discrete_R' num2str(num_sims)];

function_discrete_realisations(file_save_name, num_sims, k1, k2, a1, a2, prolif_law, ...
    beta_prolif_1, beta_prolif_2, death_law, gamma_death_1, gamma_death_2, characteristics_include, timestop, time_record_vector, N_init,...
    gaussian_n,gaussian_sigma, gaussian_mu, L,mixing_pop_id ,eta,dt,l_crit_death,q_init_condition,folder_name)



%% 2.3) Generate continuum model results
clear
clc


simulation_id = 'Fig6';
%This is for the linear proliferation and death mechanism %%

q_init_condition = 3; %mechanical equilibrium
k1=10;
k2=20;
a1=0;
a2=0;
beta_prolif_1 = 0.07;
beta_prolif_2 = 0.07;
gamma_death_1 =0.35;
gamma_death_2 =0.35;
N1=20;
N2=20;
l_crit_death=0.3;
L=10; %length of the domain
gaussian_n=4*11.056818502502322;
gaussian_sigma=sqrt(9);
gaussian_mu=5;
dx=0.01;
dt =0.0001;
eta=1;
k_hom=0;
a_hom=1;
beta_hom=1;
gamma_hom=1;
max_iters = 1000;
err_tol = 0.0001;
timestop=200;
time_vector_record=0:5:timestop;


file_save_name = ['Results_Sim' simulation_id '_Continuum'];

folder_name = [simulation_id '_Continuum'];
if ~exist([folder_name'], 'dir')
    mkdir([folder_name]);
end

one_or_two_pop=2;

tic
function_continuum_simulation(file_save_name,q_init_condition, ...
    k1,k2,a1,a2, beta_prolif_1, beta_prolif_2, gamma_death_1, gamma_death_2, l_crit_death, N1, N2, ...
    L,timestop,time_vector_record, folder_name,max_iters,err_tol,dx,dt,gaussian_n,gaussian_sigma,gaussian_mu ,eta,...
    k_hom,a_hom,beta_hom,gamma_hom,one_or_two_pop);
toc




%% 2.4) Plot characteristics - discrete
clear
clc

simulation_id = 'Fig6';

filepath_save_figs = [pwd '\' simulation_id '_Discrete\'];

load([filepath_save_figs 'Results_Sim' simulation_id '_Discrete_Characteristics.mat']);

%check that the all the data is available to plot characteristics without error. If the below returns a result decreate dt.
m_springs_per_cell=1;
function_characteristic_plot_check(num_sims,t_hist_sim_run,x_hist_sim_run,m_springs_per_cell);


%define figure properties
xticks_vec=linspace(0,L,5);
yticks_vec=[0,50,100,150,200];

%plot the characteristics
for sim_run = 1:num_sims
    
    cell_with_proliferation_event_hist=cell_with_proliferation_event_hist_sim_run{sim_run};
    cell_with_death_event_hist=cell_with_death_event_hist_sim_run{sim_run};
    cell_event_hist=cell_event_hist_sim_run{sim_run};
    prolif_death_indicator_hist=prolif_death_indicator_hist_sim_run{sim_run};
    cell_death_lr_hist = cell_death_lr_hist_sim_run{sim_run};
    
    t_hist=t_hist_sim_run{sim_run};
    x_hist=x_hist_sim_run{sim_run};
    k_hist =k_hist_sim_run{sim_run};
    a_hist =a_hist_sim_run{sim_run};
    G_hist =G_hist_sim_run{sim_run};
    mixing_pop_id_hist =mixing_pop_id_hist_sim_run{sim_run};
    
    %density colouring
    colouring =1;
    function_characteristics_plot(t_hist, x_hist, k_hist,G_hist, mixing_pop_id_hist, a_hist, L,...
        filepath_save_figs, colouring, xticks_vec, yticks_vec,....
        cell_with_proliferation_event_hist,cell_with_death_event_hist,cell_event_hist,prolif_death_indicator_hist,sim_run,cell_death_lr_hist,m_springs_per_cell)
    
    %population colouring
    colouring =4;
    function_characteristics_plot(t_hist, x_hist, k_hist,G_hist, mixing_pop_id_hist, a_hist, L,...
        filepath_save_figs, colouring, xticks_vec, yticks_vec,....
        cell_with_proliferation_event_hist,cell_with_death_event_hist,cell_event_hist,prolif_death_indicator_hist,sim_run,cell_death_lr_hist,m_springs_per_cell)
    
    close all
end



%% 2.5) Plot density snapshots, total cell number discrete, and interface position

clear
clc

simulation_id = 'Fig6';
num_sims=2000;

filepath_save_figs = [pwd '\' simulation_id '_Discrete\'];

load([filepath_save_figs 'Results_Sim' simulation_id '_Discrete_R' num2str(num_sims) '.mat']);


% 2.5.1) Density snapshots


density_evaluation_x_dx = 0.01;

plot_times=0:5:50;
function_discrete_density_plot(t_hist_sim_run, x_hist_sim_run, L, ...
    filepath_save_figs,num_sims, plot_times, density_evaluation_x_dx)

% 2.5.2) Cell stiffness snapshots
plot_times=0:5:50;
function_discrete_cellstiffness_plot(t_hist_sim_run, x_hist_sim_run,k_hist_sim_run, L, ...
    filepath_save_figs,num_sims, plot_times, density_evaluation_x_dx)

% 2.5.3) Total cell number
plot_times = 0:5:timestop;
one_or_two_pop=2;
function_discrete_total_cell_number(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims,...
    plot_times,mixing_pop_id_hist_sim_run,one_or_two_pop)

% 2.5.4) Interface position
plot_times = 0:5:timestop;
function_discrete_interfaceposition_plot(t_hist_sim_run, x_hist_sim_run, filepath_save_figs,num_sims, plot_times,mixing_pop_id_hist_sim_run)

%% 2.6) Overlay continuum density snapshots, total cell number discrete, and interface position
clear
clc

simulation_id = 'Fig6';
filepath_save_figs = [pwd '\' simulation_id '_Continuum\'];
load([filepath_save_figs 'Results_Sim' simulation_id '_Continuum.mat']);

filepath_load_figs = [pwd '\' simulation_id '_Discrete\'];


% 2.6.1) Density snapshots

colouring=1; %for density
ctm_only=0; %plot continuum overlayed on discrete
plot_times=0:5:50; %must match the plot times from the discrete times

function_discrete_ctm_comparison_snapshots_plot(t_hist, q_hist,k_hist,...
    L,dx, filepath_save_figs,  colouring,ctm_only, plot_times,filepath_load_figs)

% 2.6.2) Cell stiffness snapshots

colouring=2; %for cell stiffness
ctm_only=0; %plot continuum overlayed on discrete
plot_times=0:5:50; %must match the plot times from the discrete times

function_discrete_ctm_comparison_snapshots_plot(t_hist, q_hist,k_hist,...
    L,dx, filepath_save_figs,  colouring,ctm_only, plot_times,filepath_load_figs)

% 2.6.3) Total cell number
plot_times=0:5:timestop;
one_or_two_pop=2;
function_dis_ctm_total_cell_number_plot(t_hist, q_hist,s_hist,L,dx, filepath_load_figs,filepath_save_figs, plot_times, ctm_only,one_or_two_pop)

% 2.6.4) Interface position
plot_times=0:5:timestop;
%only for two populations
function_dis_ctm_intposition_plot(t_hist,s_hist, filepath_load_figs,filepath_save_figs, plot_times, ctm_only)
