function [Network_failure_test,...         % BINARY: test for under frequency load shedding during simulation
          extreme_freq,...
          sequence,...          % 20 X 16 array of protection system activations
          nodalFailure,...      % 4 X 6 array of protection system ativations at each node
          Dynamics_TFVGBR,...   %list of Frequency, Voltage Governor Battery and ROCOF dynamics 
          Time_to_AGC,...       %this is the time until the next agc signal
          Initial_Battery,...    % this is the initial battery output
          cumulative_attack, ... % total effective cyber attack
          incremental_attack...     %each effective attack at each time
 ] = Updated_Two_Area_Model(...
                    renewables,...         % this is the disturbance vector from the renewables
                    be,...          % normal b matrix
                    Be_Fault,...      %faulted B-matrix in case of line trip
                    m,...             % generator/motor inertia
                    LG,...            % location of generators 
                    LL,...            %location of loads
                    pw,...            %pre-contingency power vector
                    lds,...           %loads
                    nn,...            %number of nodes
                    nr,...            %number of loads
                    qa,...            %initial conditions for dynamical system
                    batcap,...        %vector of battery maximum power outputs/inputs
                    nom_freq_range,...%nominal frequency range
                    r,...             %random values for agc time and intial battery output
                    dist_location,... %a 1 X 6 vector for the location of the disturbances: either the load or generators nodes
                    attack_list,...   %this is the load attack (scaled by initial load and 100MVA base
                    vr,...            %PERCENTAGE OF LOADS SUSCEPTIBLE TO LOAD ATTACKS
                     ci,...            %attack interval
                     mlr...           %max load ratio (signifies different times of the day)
                )    
% clear all 

 
global Time_List Rocof_List volt_list RateBatPow u power regime begin Local_Gens Local_Load gen_trip load_shed_event_tally ...
elec rocof maxfreq Rocof_Thres Ulim Llim FD FN UagcSOC LagcSOC UlimSOC LlimSOC...
numbat Uagc Lagc NN M D Be mem Pb agc B0 Pbat connector_power con_relay ...
bias line_trip nextAGCsignal AGCinterval batteryperformance load_shed_schd vdeadband...
minN maxN lse_tally connect_relay_delay rocof_movingaverage BatLoc freq_delay rocof_delay...
freq_list freq_movavg agc_bat_out emergency deadband droop_u droop_d  rad_to_hertz loads line_power_limit...
 volt_delay volt_dip_threshold uvls_trip ls_times vulnerable_ratio maxload_ratio cyberattack next_cyber_attack

 
%______________ PASSING VARIABLES FOR BATTERY FUNCTION, GOVERNOR DEADBAND, AVR DEADBAND,
% POWER FLOWS AND ROCOF_______________

%BATTERY FUNCTION SEE FUNCTION BELOW
Pb = @(Pbat) passing(Pbat);    % this passes the system state at each time moment to the battery model            

%GOVERNOR DEADBAND, SEE FUNCTION BELOW
deadband = @(gov,NN,FD) fdband(gov,NN,FD);  %limits the governor output due to deadband restraints

%VOLTAGE DEADBAND
vdeadband = @(gov,NN) vdband(gov, NN);  %limits the avr response due to deadband constraints

%POWER FLOWS IN THIRD ORDER MODEL
elec =  @(v,NN,b) (v(2*NN+1:3*NN).*(b.*sin(v(1:NN)-v(1:NN)'))*v(2*NN+1:3*NN));  %power flows at each node


%SWING EQUATION OF THIRD ORDER MODEL
 

rocof = @(gov,NN,M,D,u,power,elec,Be,Pb,mn,mx,FD,Local_Gens,loads,cyberattack,vulnerable_ratio,maxload_ratio,pbat) ...
     (M)\(-D*gov(NN+1:2*NN) + max(mn',min((Local_Gens'.*(gov(6*NN+1:7*NN)) + power'),mx'))+u'...
    - (1-vulnerable_ratio)*loads'- vulnerable_ratio*loads'.*max(0,min( maxload_ratio', 1 + cyberattack'))- elec(gov,NN,Be));
% + Pb(pbat)' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------------POWER SYSTEM MODEL------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%
j2 = @(gov,u,D,M,NN,X,Eref,Be,T,power,elec,...
    Pb,aq,govdroop,FD,GT,Te,Ta,Ts,Ke,Ka,Ks,mn,mx,shutoff,Local_Gens,loads,cyberattack,vulnerable_ratio,maxload_ratio,pbat) [...
    gov(NN+1:2*NN);                                                                                                 % PHASE ANGLES
    rocof(gov,NN,M,D,u,power,elec,Be,Pb,mn,mx,FD,Local_Gens,loads,cyberattack,vulnerable_ratio,maxload_ratio,pbat);                                     % RATE OF CHANGE OF FREQUENCY (SWING EQN)
    (1./T').*(shutoff'.*(Eref' - aq'.*gov(5*NN+1:6*NN)) - gov(2*NN+1:3*NN) + X.*(Be.*cos(gov(1:NN)-gov(1:NN)'))*gov(2*NN+1:3*NN)) ;  % TERMINAL VOLTAGE EQUATION 2*NN+1:3*NN
    vdeadband(gov,NN).*shutoff'.*aq'.*(1./Ts).*(Ks*(Eref' - gov(2*NN+1:3*NN)) - gov(3*NN+1:4*NN));                              % AVR: SENSOR 
    vdeadband(gov,NN).*shutoff'.*aq'.*(1./Ta).*(Ka*( -gov(3*NN+1:4*NN)) - gov(4*NN+1:5*NN));                     % AVR AMPLIFIER
    vdeadband(gov,NN).*shutoff'.*aq'.*(1./Te).*(Ke*gov(4*NN+1:5*NN) - gov(5*NN+1:6*NN));                              % AVR EXCITER
    deadband(gov,NN,FD).*shutoff'.*aq'.*(-govdroop'.*gov(NN+1:2*NN)); %GOVERNOR
    ]; 
%______________________________________________________________________________________________________________________________

%------------------------------------------------------------------INPUT VARIABLES-----------------------------------------------------
 
Be = be; M=m;Local_Gens = LG;Local_Load = LL;power = pw;NN = nn; loads = lds;
vulnerable_ratio = vr; maxload_ratio =mlr; cyber_interval = ci;
       
BatLoc = [5 6];
rad_to_hertz = 2*pi;            %CONVERTS RADIAN TO HERTZ MEASURE
pre_dist_loads = loads;

D = 2/(2*pi*60);                %dampening coefficients for each node
% fixed_D = D;
dM = diag(M(1:4,1:4));                             

Eref = ones(1,NN)*1;
X = (1.8-0.3)*(100/900)*(20/230)^2;                %resistance
T=8;

%DISTURBANCE VECTORS........................................................
% u = disturbance(renewables,NN,dist_location).*loads ;           %disturbance vector 
u = zeros(1,NN);% [0 0 0 0 0.1 0.01].*[0.1 0 0 0 power(5:6)];
load_attack = attack_list(1,:);


%INPUTS FOR PROTECTION SYSTEMS.............................................
connector_power = 4;
connect_relay_delay = 4;
line_power_limit = 5.1;   
maxfreq = 2;       
Rocof_Thres = 3;                                       %abs rocof threshold
freq_delay  = 2;
rocof_delay = 1;
volt_delay = 1.5;
volt_dip_threshold = 0.1;
maxN = Local_Gens*9;%+ Local_Load*0;
minN = 0;%Local_Gens*0 + Local_Load*-20;
load_shed_schd = [ -1 -1.5 -2 -2.5 -3 ];  
% 
shutoff = ones(1,NN);      %THIS IS A BINARY VARIABLE: INTIAL VALUE IS 1 AT ALL NODES. ...
                           % WHEN A GENERATOR DISCONNECTS DUE TO ROCOF OR
                           % OFGS, IT SETS FIELD VOLTAGE AT THAT NODE TO ZERO; GOVERNOR
                           % DIFFERENTIAL EQUATION AND OUTPUTS TO ZERO; AND
                           % THE AVR SENSOR, EXCITER AND AMPLIFIER OUTPUTS
                           % AND DIFFERENTIAL EQUATIONS TO ZERO.
                   
%.........................................................................

%governor model............................................................

aq = Local_Gens;           % places governors at generators
govdroop =  2;
GT = 0.25;

%...........................................................................

%AGC SETTINGS AND PARAMETERS..............................................

agclim = 0.5;                                           %abs limit for battery power to be provided by AGC
AGCinterval = 2;                                      %time interval for the impulse signal for AGC
mem = 0;                                                %required to memorize the previous time step
agc = 0;                                                %agc indicator: 0: AGC off; 1: AGC ON
nextAGCsignal = r(1)*AGCinterval;
Time_to_AGC = nextAGCsignal;

%BATTERY PARAMETERS AS INPUTS..............................................
RateBatPow = zeros(1,NN);
RateBatPow(BatLoc) = batcap;                                  %battery capabilility
B0 = agclim*(-0.05 + 0.1*r(2:end)).*RateBatPow;            %initial battery operating point
Initial_Battery = B0;
agc_bat_out = B0(BatLoc);
Ulim =   RateBatPow;                                          %upper limit of the capability of the battery                           
Llim = - RateBatPow;                                         %lower limit of the capability of the battery
numbat = sum(Local_Load);                                   %number of batteries in the network (required for the AGC)
FD = 0.05;                                                   %deadband frequency range
FN = 0.15;                                                   %maximum AGC limit
Fmax = 1;
Uagc =  agclim*RateBatPow;                                   %upper AGC power limit
Lagc = -agclim*RateBatPow;                                  %lower limit for battery power commanded by AGC
droop_u = (RateBatPow(BatLoc)-agc_bat_out)./(Fmax-FN);      %required to calculate battery droop response    
droop_d = (RateBatPow(BatLoc)+agc_bat_out)./(Fmax-FN);
Pbat = B0;                                                  %CORrrectedBattery Power following the AGC signal
batteryperformance = B0;
bias = 1.5;
emergency = zeros(1,NN);

%STATE OF CHARGE %VARIABLES................................................

% BatteryCapacity = 55;                                       %total battery capacity
% SOC = rand(1,NN)*BatteryCapacity;                           %initial state of charge of each battery
UagcSOC = Uagc(BatLoc);LagcSOC = Lagc(BatLoc);
UlimSOC = Ulim(BatLoc);LlimSOC = Llim(BatLoc);                %initializing the limits under State of Charge

%..........................................................................

%AVR MODEL PARAMETERS......................................................

Te = 1;
Ke = 10;
Ta = 0.1;
Ka = 10;
Ts = 0.05;
Ks = 1;

%INITIALISING LIST VARIABLES...............................................

gen_trip = ~Local_Gens;
line_trip = zeros(1,1); %there are 4 interconnectors to be tripped
Time_List = [];  
Rocof_List = [];
volt_list = [];
freq_list =[];
rocof_tally = zeros(1,NN);
OFGS_tally = zeros(1,NN);
UFLS_tally = zeros(1,NN);
uvls_trip = zeros(1,NN);

load_shed_event_tally = zeros(1,NN);
lse_tally = zeros(1,NN);
con_relay = zeros(1,1);
ls_times = zeros(1,NN);
load_shed_ufls = 0.15*loads;       %load shed amount at each node
load_shed_uvls = 0.1.*loads;

%INITIAL PARAMETERS FOR THE ODE SOLVER
simulation_duration = 60;
points2 =  60;
dt2 = simulation_duration/points2;
interval = 0:dt2:simulation_duration + 0.5;

cinterval = cyber_interval:cyber_interval:simulation_duration;
begin = 0;
next_cyber_attack = min(cinterval(cinterval>begin));      
 
attack_counter = 1;
current_vul_loads = 100*vulnerable_ratio*Local_Load.*loads.*(max(0,min( maxload_ratio, 1 + load_attack )));

effective_attack =  100*vulnerable_ratio*Local_Load.*loads.*(max(0,min( maxload_ratio, 1 + load_attack ))-1);

cumulative_attack = abs(effective_attack);
incremental_attack = [0 effective_attack]; % 0 for the time t

% commanded_incremental_attack = load_attack;
initial = qa;
record = [];
options = odeset('AbsTol',1e-5,'Events',@ProtectionTrip);
sequence = zeros(60,4*NN+length(line_trip)+1);
cyberattack = load_attack;

j=0;
while begin < simulation_duration
 
j = j+1;
regime = 1; 
[td,bat,~,~,failuretype] = ode45(@(td,bat) j2(bat,u,D,M,NN,X,Eref,Be,T,power,elec,Pb,aq,...
    govdroop,FD,GT,Te,Ta,Ts,Ke,Ka,Ks,minN,maxN,shutoff,Local_Gens,loads,cyberattack,...
    vulnerable_ratio,maxload_ratio,Pbat),...
    interval,initial,options);
 
[rocof_loc,gen_shed_loc,ufls_event,con_relay,uvls_event] = classify_failure(failuretype,NN);          % this functions classifies each type of failure and its location and saves the result in the ooutput variables
record = [record(1:end,:); td bat(:,[NN+1:3*NN,6*NN+1:7*NN])];
sequence(j,:) = [td(end), rocof_loc,gen_shed_loc,ufls_event,uvls_event,con_relay];
gen_trip = gen_trip + rocof_loc + gen_shed_loc; % this is a variable used to record the number of generator trips
line_trip  = line_trip + con_relay ;
rocof_tally = rocof_tally + rocof_loc;
OFGS_tally  = OFGS_tally  + gen_shed_loc;
UFLS_tally = UFLS_tally + ufls_event;
lse_tally = lse_tally + ufls_event;
uvls_trip = uvls_trip + uvls_event;
power(rocof_loc|gen_shed_loc) = 0;               %turns the power off at generators which experienced rocof or ofgs events
ls_times = max(ls_times,ufls_event*td(end)); %generates the most recent times of load shed events

if con_relay == 1
    Be = Be_Fault;
    d_M([1:3]) = sum(dM(1:2));      %if generator is shed AFTER a line trip, calculate new inertia for area 1
    d_M([4:6]) = sum(dM(3:4));      %if generator is shed AFTER a line trip, calculate new inertia for area 2
    M = diag(d_M);                  %M is now a 6 X 6 matrix
    D_mat(1:3) = 2*sum(loads(1:3))/sum(pre_dist_loads(1:3)); %D for area 1 after line trip
    D_mat(4:6) = 2*sum(loads(4:6))/sum(pre_dist_loads(4:6));  
    D = diag(D_mat);
end

if sum(rocof_loc+gen_shed_loc)>0        %tests for generation shedding event
    dM(rocof_loc|gen_shed_loc) = 0.04;  %set the inertia constant at the lost generator to approx 0;

    if line_trip > 0
        d_M([1:3]) = sum(dM(1:2));   %if generator is shed AFTER a line trip, calculate new inertia for area 1
        d_M([4:6]) = sum(dM(3:4));   %if generator is shed AFTER a line trip, calculate new inertia for area 2
        M = diag(d_M);               %M is now a 6 X 6 matrix
    else     
        M = sum(dM);                    % if there is no line trip, M is the sum of the 4 inertias
    end
end

if sum(uvls_event)>0
    loads = loads - load_shed_uvls.*uvls_event; %new load after under voltage load shed
end

if sum(ufls_event)>0
    loads = loads - ufls_event.*load_shed_ufls;   %new load after load shed event
end

if sum(ufls_event+uvls_event)>0
% loads = loads - load_shed_event.*load_shed;   %new load after load shed event

    if line_trip>0
         D_mat(1:3) = 2*sum(loads(1:3))/sum(pre_dist_loads(1:3)); %D for area 1 after line trip
         D_mat(4:6) = 2*sum(loads(4:6))/sum(pre_dist_loads(4:6));  
         D = diag(D_mat);
    else
         D = D*(sum(loads)/sum(pre_dist_loads)); 
    end
end

if any(failuretype==4*NN+2) %this is when the next attack is scheduled to occur
    attack_counter = attack_counter +1;
    load_attack = attack_list(attack_counter,:).*Local_Load;
    new_vul_loads = 100*vulnerable_ratio*Local_Load.*loads.*(max(0,min( maxload_ratio, 1 + load_attack )));

    cumulative_attack  = cumulative_attack + abs( new_vul_loads - current_vul_loads );
    incremental_attack = [incremental_attack; round(td(end),5) new_vul_loads-current_vul_loads]; %new - old = inc/dec
    cyberattack = load_attack;

    current_vul_loads = new_vul_loads;

end

shutoff(rocof_loc|gen_shed_loc) = 0;
initial = [bat(end,1:3*NN) [shutoff shutoff shutoff shutoff].*bat(end,3*NN+1:7*NN)];%bat(end,7*NN+1:8*NN)];% 
begin =  td(end) + 0.05;
interval = [begin,interval(interval>begin)];
next_cyber_attack = min(cinterval(cinterval>begin));
end

nodalFailure = [rocof_tally 0; OFGS_tally 0; UFLS_tally 0; uvls_trip 0; zeros(1,NN) line_trip];  
nodalFailure(end+1,:) = sum(nodalFailure,1);
sequence = sequence(1:end-1,:);
Network_failure_test = any(nodalFailure>0,'all');%sum(UFLS_tally,'all')>0;
Mov_RoCoF_TS = interp1(Time_List,rocof_movingaverage,record(:,1)); %rocof for time steps above (interpolated)
Mov_Freq_TS = interp1(Time_List,freq_movavg,record(:,1)); 
battery_TS = interp1(Time_List,batteryperformance,record(:,1));    %battery output for time steps above (interpolated)

Dynamics_TFVGBR  = [record(:,1),Mov_Freq_TS+60,record(:,NN+2:end),battery_TS, Mov_RoCoF_TS];
extreme_freq = [( (min(Dynamics_TFVGBR(:,2:NN+1),[],'all')<nom_freq_range(1)) | (max(Dynamics_TFVGBR(:,2:NN+1),[],'all')>nom_freq_range(2))  )>0, max(Dynamics_TFVGBR(:,2:NN+1),[],'all'), min(Dynamics_TFVGBR(:,2:NN+1),[],'all')];

end

function[position,term,direction]= ProtectionTrip(t,sol)

global NN gen_trip Local_Gens  power_change line_power_limit...
    maxfreq Rocof_Thres lse_tally load_shed_schd line_trip ...
    cyber_attack  Local_Load volt_dip_threshold uvls_trip

gen_shed =  @(s) (s.*Local_Gens)<maxfreq;
load_shed2 = @(s,mf) (Local_Load.*s)>mf;      % for load frequency activated UFLS
uvload_shed = @(s) (s<volt_dip_threshold);   
Construct_List(t,sol);         %launches the construct_list function which lists the time and rocof values from each step

[adjma, ma_freq, ma_volt] = MovAvg();           %calculates the 1s moving average for ROCOF
mf = load_shed_schd(lse_tally+1);

line_freq = sum(abs(ma_freq(5:6))>2);
ICT = power_change<line_power_limit & line_freq<1; %test for the line tripping

            %line power             OR     %line frequencies
ICT = power_change<line_power_limit | line_freq;                %test for the line tripping.this must be an OR because it should switch to false only when both conditions are false

position = [((adjma)<Rocof_Thres)+gen_trip,...
            gen_shed(ma_freq)+gen_trip,...
            load_shed2(ma_freq,mf)+(lse_tally>=4),...
            uvload_shed(ma_volt)+uvls_trip,...
            ICT+line_trip,...
            cyber_attack ]';

% find(~position)
term     = [ones(1,NN),ones(1,NN),ones(1,NN) ones(1,NN), 1 1]';
direction = [];   
 
end

function [AdjMA,Adj_FMA,Adj_VMA_dev] = MovAvg()

   global begin Local_Gens Local_Load  Time_List Rocof_List rocof_delay...
          freq_list freq_delay rocof_movingaverage freq_movavg ...
          volt_list volt_delay volt_mov_avg ls_times
   

window_rocof = Time_List(end)>=rocof_delay;
% window_freq = Time_List(end) >=freq_delay;
window_volt = Time_List(end) >=volt_delay;

rocof_movingaverage = movmean(Rocof_List,[rocof_delay,0],1,'SamplePoints',Time_List,'Endpoints','shrink');
volt_mov_avg = movmean(volt_list, [volt_delay,0],1,'SamplePoints',Time_List,'Endpoints','shrink');
freq_movavg = movmean(freq_list,[freq_delay,0],1,'SamplePoints',Time_List,'Endpoints','shrink');


AdjMA   = (rocof_movingaverage(end,:).*Local_Gens)*window_rocof*(Time_List(end)>begin); %captures the moving average of rocof
%if there is an activation of ufls, the counter resets for another delay
%before triggering again
Adj_FMA = freq_movavg(end,:).*(Time_List(end)>(ls_times+freq_delay)); %captures the moving average of rocof

Adj_VMA_dev = (1-max(0,volt_mov_avg(end,:))).*Local_Load*(Time_List(end)>begin)*window_volt;
 
end

function Construct_List(t,sol)
global Time_List Rocof_List u power regime rocof...
    elec Pb NN M D Be minN maxN freq_list Local_Gens connector_power connect_relay_delay...
    FD power_change con_relay volt_list rad_to_hertz loads cyberattack...
    vulnerable_ratio maxload_ratio batteryperformance Pbat ...
    cyber_attack next_cyber_attack

if t > 0 && t <= Time_List(end) 
   regime = 0;   
elseif regime == 0 && (t-Time_List(end))>0.1
    regime = 1;
end

if regime == 0 
    
   Time_List(end) = t;
   Rocof_List(end,:) = abs(rocof(sol,NN,M,D,u,power,elec,Be,Pb,minN,maxN,FD,Local_Gens,loads,cyberattack,vulnerable_ratio,maxload_ratio,Pbat))'/rad_to_hertz;
   freq_list(end,:) = sol(NN+1:2*NN)'/rad_to_hertz;
   volt_list(end,:) = sol(2*NN+1:3*NN);
   connector_power(end) = abs(sol(2*NN+5)*sol(3*NN)*Be(5,6)*sin(sol(5)-sol(6)));

   power_change = (connector_power(end,:) - connector_power(1,:))*(t>connect_relay_delay);
   batteryperformance(end,:)= Pbat;
   
else
   Time_List  = [Time_List; t]; 
   Rocof_List = [Rocof_List; abs(rocof(sol,NN,M,D,u,power,elec,Be,Pb,minN,maxN,FD,Local_Gens,loads,cyberattack,vulnerable_ratio,maxload_ratio,Pbat))'/rad_to_hertz];

   freq_list =  [freq_list; sol(NN+1:2*NN)'/rad_to_hertz];
%  
   volt_list = [volt_list;sol(2*NN+1:3*NN)'];

   
   if length(Time_List)>1
   BatteryModel(sol);                         %passes t and sol to the AGC function
   batteryperformance = [batteryperformance; Pbat];
   end

   if con_relay<1
   connector_power = [connector_power; abs(sol(2*NN+5)*sol(3*NN)*Be(5,6)*sin(sol(5)-sol(6)))];
   
   if t>connect_relay_delay
       connector_power(1) = [];
   end
       power_change = (connector_power(end) - connector_power(1))*(t>connect_relay_delay);
   end
end

if t>next_cyber_attack
    cyber_attack = 0;
else
    cyber_attack = 1;
end

end

function [rocof_loc,gen_shed_loc,load_shed_event,conn_relay_trip,uvls_load_shed_event] = classify_failure(failuretype,NN)

%INITIALIZE THE FAILURE VARIABLES
rocof_loc = zeros(1,NN); conn_relay_trip = 0; gen_shed_loc = zeros(1,NN);
load_shed_event=zeros(1,NN);uvls_load_shed_event= zeros(1,NN);

rocof_loc(failuretype(failuretype<=NN)) = 1;    %this calculates where rocof failures occur
gen_shed_loc(failuretype(failuretype >= (NN+1) & failuretype<= 2*NN)-NN) = 1; % this calcualtes where ofgs occurs
load_shed_event(failuretype(failuretype>2*NN & failuretype<=3*NN)-2*NN)  = 1;  % this calculates where ufls occurs
uvls_load_shed_event(failuretype(failuretype>3*NN & failuretype<=4*NN)-3*NN)  = 1;  % this calculates where uvls occurs
if failuretype == 4*NN+1
    conn_relay_trip = 1;  %this calculates where line trips occurs
end

if sum(rocof_loc+gen_shed_loc>1)>0
    rocof_loc((rocof_loc+gen_shed_loc)>1)=0;        %where both ofgs and rocof failures occur at the sme time, record as ofgs by setting rocof failure as 0
end

end


function BatteryModel(sol)
global FD FN B0 nextAGCsignal AGCinterval E_FCAS rad_to_hertz...
       NN Time_List agc Pbat BatLoc agc_bat_out R_FCAS 

freq = sol(NN+1:2*NN)'/rad_to_hertz;         % frequency deviation at time t
abs_freq = abs(freq);           % absolute value of the frequency deviation

E_FCAS = abs_freq(BatLoc)>=FN;  % test to see if frequency deviation is in the emergency band
R_FCAS = abs_freq(BatLoc)<FN  ; % indicates whethere frequency deviation at battery nodes are within agc command range (1)

%AGC Signal and AGC Command
if Time_List(end) >= nextAGCsignal
    nextAGCsignal = nextAGCsignal + AGCinterval; %sets the next time agc signal    
    agc = mean(freq)>FD;  %measures whether the average frequency across the network requires agc action - (1). Else 0
    agc_ref_freq = mean(freq);    % signal to be sent to ALL batteries 
    
   if agc>0 %frequency deviation is in the response range
       agc_bat_out = AGCBat(agc_ref_freq); %battery power under the agc command 
   else
       agc_bat_out = B0(BatLoc); % battery power is set to pre disturbance power when agc is 0. 
   end   
end

Pbat(BatLoc) =  EmergencyBat(freq,E_FCAS) + agc*(~E_FCAS).*agc_bat_out + (~agc)*(~E_FCAS).*B0(BatLoc);
% batteryperformance = [batteryperformance;Time_List(end) Pbat];

end

function Ebat = EmergencyBat(f,emergency) 

global UlimSOC LlimSOC droop_u droop_d B0 BatLoc FN
    
    local_bat_freq = f(BatLoc);
    ad_f = abs(local_bat_freq) - FN; %adjusted frequency by removing FN
    n_f = ad_f.*sign(local_bat_freq); %adjusted frequency with original sign
    ud = sign(local_bat_freq)>0;
    droop = droop_u.*(~ud) + droop_d.*(ud);
    Eb = (-droop.*(n_f) + B0(BatLoc)).*emergency; %emergency battery response on top of initial battery point
    Ebat = Eb.*(Eb>=LlimSOC & Eb<=UlimSOC) + LlimSOC.*(Eb<LlimSOC) + UlimSOC.*(Eb>UlimSOC); %with limits

end

function Abat = AGCBat(f)  
   global UagcSOC LagcSOC bias B0 BatLoc
   
   Agcbat = -bias*f + B0(BatLoc);
   Abat =  Agcbat.*(Agcbat>=LagcSOC & Agcbat<=UagcSOC) + LagcSOC.*(Agcbat<LagcSOC) + UagcSOC.*(Agcbat>UagcSOC);
end

 function w = fdband(gov,NN,FD)
global rad_to_hertz
w = abs(gov(NN+1:2*NN)/rad_to_hertz)>0.01;
end

function vd = vdband(gov,NN)
vd = abs(1 - gov(2*NN+1:3*NN)) > 0.05;
end

function v = passing(Pbat)
% global Pbat
v = Pbat;
end 

function u = disturbance(dist,NN,dist_loc)

u = zeros(1,NN);
u(dist_loc) = dist;
end

function B = line_disconnection(con_relay,B,position)
        matches = {[1],[2 3 4], [5], [6 7]}; %matches line trip code to loction in network
        line_fail_ind = [matches{find(con_relay)}];

        B(position(line_fail_ind,:))=0;
%         
%         if con_relay(6)==1
%             B(position(7,:)) = 0;
%         end
end

function [M,D] = islanding(line_trip,all_M,fixed_D,loads,areas,pre_dist_loads,Local_Gens)
    test = line_trip.*[3 5 7 11];%just prime numbers for unique additive identities later on
    
    if sum(test)==15 %[1 2 3] area 1 only islanded
        all_D(areas==1) = fixed_D*sum(loads(areas==1))/sum(pre_dist_loads);
        all_D(areas~=1) = fixed_D*sum(loads(areas~=1))/sum(pre_dist_loads);

        all_M(areas==1) = sum(all_M(Local_Gens & areas==1));
        all_M(areas~=1) = sum(all_M(Local_Gens & areas~=1));
        M = diag(all_M);
        D = diag(all_D);
    elseif sum(test)==14 %[1 4]area 2 only islanded

        all_D(areas==2) = fixed_D*sum(loads(areas==2))/sum(pre_dist_loads);
        all_D(areas~=2) = fixed_D*sum(loads(areas~=2))/sum(pre_dist_loads);

        all_M(areas==2) = sum(all_M(areas==2 & Local_Gens ));
        all_M(areas~=2) = sum(all_M(areas~=2 & Local_Gens));

        M = diag(all_M);
        D = diag(all_D);
    elseif sum(test)== 23 %[2 3 4] area 3 only islanded

        all_D(areas==3) =  fixed_D*sum(loads(areas==3))/sum(pre_dist_loads);
        all_D(areas~=3) =  fixed_D*sum(loads(areas~=3))/sum(pre_dist_loads);

        all_M(areas==3) = sum(all_M(areas==3 & Local_Gens));
        all_M(areas~=3) = sum(all_M(areas~=3 & Local_Gens));
        M = diag(all_M);
        D = diag(all_D);
    elseif sum(test)==sum([3 5 7 11]) % all areas islanded

        all_D(areas==1) =  fixed_D*sum(loads(areas==1))/sum(pre_dist_loads);
        all_D(areas==2) =  fixed_D*sum(loads(areas==2))/sum(pre_dist_loads);
        all_D(areas==3) =  fixed_D*sum(loads(areas==3))/sum(pre_dist_loads);
    
        all_M(areas==1) = sum(all_M(areas==1 & Local_Gens));
        all_M(areas==2) = sum(all_M(areas==2 & Local_Gens));
        all_M(areas==3) = sum(all_M(areas==3 & Local_Gens));
        M = diag(all_M);
        D = diag(all_D);
    else
        M = sum(all_M);
        D = fixed_D*sum(loads)/sum(pre_dist_loads);
    end
end