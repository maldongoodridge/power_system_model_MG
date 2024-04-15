clear all
tic

load('IEEE_inputs_27_nodes.mat')
load('IEEE_initial27.mat')
% start
%INPUTS FOR STATISTICAL MODEL MCMC
start = 0;
number = 10; %numberof trials
tuning = 0.1;  %adjust this numberso that the 'rate' which is printed is around 0.15 - 0.4
tuning_attack =0.2;
network_skips = 3;
particle_skips = 6;
qq = 1;

%_POWER SYSTEMS MODEL INPUTS_________________________________________________________________________________________________________
bat_list = 0;                   %battery maximum power output list IN HUNDREDS OF MW
dynamic_or_static=1;                %Enter 0 for static attack; Enter 1 for dynamic attack 
nom_freq_range = [59.85, 60.15];     % this is the nominal frequency band for the frequency regulation MCMC

%____________________________________________________________________________________________________________________

for qq =  1:length(bat_list)
 
dist_location = logical(Local_Load+Local_Gens);
label = ['Full_Attack'];

initial_network_value = [
    0.54...     % vulnerable_ratio for loads
    4 ...       % time of day
    5 ...       % interval between attacks 
    ];
initial_attack = 5*rand(386,19);

NR = length(initial_network_value);     %NUMBER OF NODES WITH RENEWABLES DO NOT ALTER
Adjusted_Skipping_Algorithm_IEEE39_particles_cyclical(...
    qq,tuning,start,NR,NN,number,B,Local_Gens,Local_Load,M,...
    power,loads,areas,Gen_M,interconnectors,line_index,T,...
    initial, initial_network_value,nom_freq_range,...
    dist_location,maximum_loads_ratios,tuning_attack,initial_attack,network_skips,particle_skips) 

end 
toc