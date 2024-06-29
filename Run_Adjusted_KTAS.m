clear all
tic

load('Two_Area_Inputs_Cap_RLoad_RGen.mat')
load('Two_Area_initial.mat')

% start
%INPUTS FOR STATISTICAL MODEL MCMC
start = 0;
number = 50000;%numberof trials
network_tune = 0.1;
tuning = 0.1; %adjust this numberso that the 'rate' which is printed is around 0.15 - 0.4
skips = 3
network_skips =3;
qq = 1;

%_POWER SYSTEMS MODEL INPUTS_________________________________________________________________________________________________________
bat_list = 0;                   %battery maximum power output list IN HUNDREDS OF MW
dynamic_or_static=1;                %Enter 0 for static attack; Enter 1 for dynamic attack 
nom_freq_range = [59.85, 60.15];     % this is the nominal frequency band for the frequency regulation MCMC

%____________________________________________________________________________________________________________________

for qq =  1:length(bat_list)
 
dist_location = logical(Local_Load+Local_Gens);
label = [ 'Full_Attack' ];

initial_MCMC_value = [
    0.84...     % vulnerable_ratio for loads
    4,...   % time of day
    5 ...       % interval between attacks 
     ];
 initial_attack = -1+11*rand(13,6);
NR = length(initial_MCMC_value);     %NUMBER OF NODES WITH RENEWABLES DO NOT ALTER
% Adjusted_Skipping_Algorithm_KTAS(qq,tuning,start,NR,NN,number,Be,Be_Fault,Local_Gens,Local_Load,M,power,loads,...
%     initial_conditions,...
%     initial_MCMC_value,...
%     nom_freq_range,...
%     dist_location,...
%     maximum_load_ratios,skips) 
% 
Adjusted_Skipping_Algorithm_KTAS_particles(qq,tuning,start,NR,NN,number,Be,Be_Fault,Local_Gens,Local_Load,M,power,loads,...
    initial_conditions,...
    initial_MCMC_value,...
    nom_freq_range,...
    dist_location,...
    maximum_load_ratios,initial_attack,skips,network_skips,network_tune);
end 
toc