function Adjusted_Skipping_Algorithm_KTAS_particles(q,tune,start,NR,NN,number,Be,Be_Fault,Local_Gens,Local_Load,M,power,loads,...
initial_conditions,initial_value,nom_freq_range,...
dist_location,maximum_load_ratios,initial_attack,k,network_counter,network_tune)       %vector of battery location in the network)

dimension = length(initial_value);           
if start ==1

    n = num2str((q-1)*100);
    y = load(['Cyber_Attack',n,'MW_BatLogN4.mat'],'DATA');
    
    network_attacks = y.DATA.(['D',n]);
    current_network_attack = network_attacks{end}{1};
    ROU_Sequence= y.DATA.(['BatSequence',n]);
    Nodal_Activation= y.DATA.(['BatNodeAct',n]);
    AGC_Time = y.DATA.(['AGC_Initial_Time',n]);
    Pre_Contin_Battery = y.DATA.(['Initial_Bat',n]);
    Dynamics = y.DATA.(['Bat_Dynamics',n]);
    Renewables_Fluctuations= y.DATA.(['renewables',n]);
    Extreme_Frequencies = y.DATA.(['extreme_frequencies',n]);
    Total_Effective_Cyber_Attack = y.DATA.(['effective_attack',n]);
    Cyber_Attack_Sequence = y.DATA.(['attack_sequence',n]);
    Commanded_Attack_Sequence = y.DATA.(['commanded_attack_sequence',n]);
    time_of_day = y.DATA.(['Time_of_Day',n]);
    current_cumulative_attack = Total_Effective_Cyber_Attack(end,:) ;
else
    current_network_attack = initial_value;
    current_cumulative_attack = 5190*ones(1,NN); %choose a large number with low density in log normal
    network_attacks={};%zeros(1,dimension);
    Nodal_Activation= [];
    ROU_Sequence= {}; 
    AGC_Time = [];
    Pre_Contin_Battery = [];
    Dynamics = {};
    Renewables_Fluctuations = [];
    Extreme_Frequencies = [];
    Total_Effective_Cyber_Attack =[];
    Cyber_Attack_Sequence = {};
    Commanded_Attack_Sequence = {};
    time_of_day =[];
    current_particles = initial_attack;
end

jj = 0;

for i = 1:number

    proposal = current_network_attack;%[rand(),randi([1,4],1),1+59*rand()];
    
    network_trajectory= randn(1,NR);
    network_direction = network_trajectory./vecnorm(network_trajectory,2,2);
    

    counter = 1;
    
    while counter<=network_counter

        particles = current_particles;
        old_no_attacks = size(particles,1);

        network_distance = exprnd(network_tune,1);
        proposal = proposal + network_direction*network_tune;
        proposal = cyclical_box(proposal,[0 0 1/60],[1,4,1],"network");
        attack_interval = 60*proposal(3);

        diurnal_cycle = ceil(proposal(2));
        no_of_attacks = ceil(60/attack_interval);
        
        if old_no_attacks < no_of_attacks
            extra = no_of_attacks-old_no_attacks;
            repeat_index = randi([1 old_no_attacks],extra,1);
            particles = [particles;particles(repeat_index,:)];
        else
            particles = particles(1:no_of_attacks,:);
        end

        %     
        %       particles = -1+11*rand(no_of_attacks,6);
        
        %     k = 5;
        trajectories = randn(no_of_attacks,NN);
        particle_directions = trajectories./vecnorm(trajectories,2,2);
        subcounter=1;
        
        while subcounter<=k
        
        
            particle_jump =  exprnd(tune,no_of_attacks,1);%tuning_attack*rand(no_of_attacks,1);% 
            a = particles + particle_jump.*particle_directions;
            particles = a;
            particles = cyclical_box(particles,0,10,"particles");
    
            r = rand(1,7).*[1 zeros(1,6)];
            renewables = 0;%*mvnrnd(zeros(1,NN),eye(NN));
            
            attacks = zeros(no_of_attacks,NN);
            attacks(:,Local_Load) = particles;
            
            [Network_Failure_Test,...     % BINARY: test for under frequency load shedding during simulation
            extreme_frequencies,...
            Failure_sequence,...          % 20 X 16 array of protection system activations
            Nodal_Failures,...            % 4 X 6 array of protection system ativations at each node
            Dynamics_TFVGBR,...           % list of Frequency, Voltage Governor Battery and ROCOF dynamics 
            Time_to_AGC,...               % this is the time until the next agc signal
            Initial_Battery_State,...     % this is the initial battery output
            Cumulative_Attack,...         % total effect attack during simulation
            Incremental_Attack...         % each effective cyber attack during simulation
            ] = Updated_Two_Area_Model(...
                renewables,...          %this is the disturbance vector from the renewables
                Be,...                  %normal b matrix
                Be_Fault,...            %faulted B-matrix in case of line trip
                M,...                   %generator/motor inertia
                Local_Gens,...          %location of generators 
                Local_Load,...          %location of loads
                power(diurnal_cycle,:),...               %pre-contingency power vector
                loads(diurnal_cycle,:),...               %loads
                NN,...                  %number of nodes
                NR,...                  %number of loads
                initial_conditions(diurnal_cycle,:),...                  %initial conditions for dynamical system
                0,...%*proposal(8:9),...       %vector of battery maximum power outputs/inputs
                nom_freq_range,...      %nominal frequency range
                r,...                   %random values for agc time and intial battery output
                dist_location,...       %the location of the disturbances (loads or generators)
                attacks,...   	    %this is the initial load attack (scaled by initial load and 100MVA base
                proposal(1),  ...       %ratio of vulnerable load
                attack_interval,...         %interval of cyber attack
                maximum_load_ratios(diurnal_cycle,:) ... %max load ratio
                ) ;    
            
            if Network_Failure_Test> 0 
                Region = 1;
                break
            else
                Region = 0;
            end
            
            subcounter = subcounter+1;
            
        end
        
        if Region==1
            break
        end

        counter = counter+1;
    end
    
    prop_cumulative_attack = Cumulative_Attack;
    
      accept = Region*prod(lognpdf(abs(prop_cumulative_attack(Local_Load)),1,4))/...
      (prod(lognpdf(abs(current_cumulative_attack(Local_Load)),1,4)));
%     
   %  accept = Region*lognpdf(sum(prop_cumulative_attack(Local_Load)),1,4)/...
  %           (lognpdf(sum(current_cumulative_attack(Local_Load)),1,4));	

    
    if accept > rand()
        jj = jj+1;
        current_network_attack = proposal;
        current_particles = particles;
        time_of_day(end+1,:) = diurnal_cycle;
        current_cumulative_attack = prop_cumulative_attack;
        network_attacks{end+1} = proposal.*[1 1 60];
        ROU_Sequence{end+1} = Failure_sequence;
        Nodal_Activation(:,:,end+1) = Nodal_Failures;
        Dynamics{end+1} = Dynamics_TFVGBR;
        AGC_Time(end+1) = Time_to_AGC;
        Pre_Contin_Battery(:,end+1) = Initial_Battery_State;
        Renewables_Fluctuations(:,end+1) = renewables;
        Extreme_Frequencies(end+1,:) = extreme_frequencies;
        Total_Effective_Cyber_Attack(end+1,:) = prop_cumulative_attack;
        Cyber_Attack_Sequence{end+1} = Incremental_Attack;
        Commanded_Attack_Sequence{end+1} = attacks;
    end
    
    if mod(i,100 )==0 
        rate = jj/i
        i
        save_IEEE(q,...
        network_attacks,...
        ROU_Sequence,...
        Nodal_Activation,...
        Dynamics,...
        AGC_Time,...       %this is the time until the next agc signal
        Pre_Contin_Battery,...  
        Renewables_Fluctuations,...
        Extreme_Frequencies,...
        Total_Effective_Cyber_Attack,...
        Cyber_Attack_Sequence,...
        Commanded_Attack_Sequence,...
        time_of_day ...
        );
    end

end

end



% end

function x =  cyclical_box(x,minn,maxx,type)

    mintest = @(x,minn,maxx) maxx - mod(minn-x,(maxx-minn)+0.00001);
    maxtest = @(x,maxx,minn) minn+mod(x-maxx,(maxx-minn)+0.00001);

    maxpos = find(x>maxx);
    minpos =find(x<minn);
    if type=="network"
        x(maxpos) = maxtest(x(maxpos),minn(maxpos),maxx(maxpos));
        x(minpos) = mintest(x(minpos),minn(minpos),maxx(minpos));

%         x(maxpos) = minn(maxpos)+(x(maxpos)-maxx(maxpos));
%         x(minpos) = max(minpos) - (minn(minpos)-x(minpos));
    else
        x(maxpos) = maxtest(x(maxpos),minn,maxx);
        x(minpos) = mintest(x(minpos),minn,maxx);
%         x(maxpos) = minn+(x(maxpos)-maxx);
%         x(minpos) = maxx- (minn-x(minpos));
    end

end

