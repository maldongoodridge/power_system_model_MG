function save_IEEE(q,...
        network_parameters,...
        BAT_ROU_Sequence,...
        BAT_Nodal_Activation,...
        BAT_Dynamics,...
        Time_to_AGC,...       %this is the time until the next agc signal
        Initial_Battery,...  
        renewables,...
        extreme_freq,...
        effective_attack,...
        attack_sequence,...
	    commanded_attacks,...
        time_of_day ...
        )

n =   num2str((q-1)*100);

    DATA.(['D',n])  = network_parameters;
    DATA.(['BatSequence',n]) = BAT_ROU_Sequence;
    DATA.(['BatNodeAct',n]) = BAT_Nodal_Activation;
    DATA.(['AGC_Initial_Time',n]) = Time_to_AGC;
    DATA.(['Initial_Bat',n]) = Initial_Battery;
    DATA.(['Bat_Dynamics',n]) = BAT_Dynamics;
    DATA.(['renewables',n]) = renewables;
    DATA.(['extreme_frequencies',n]) = extreme_freq;
    DATA.(['effective_attack',n]) = effective_attack;
    DATA.(['attack_sequence',n]) = attack_sequence;
    DATA.(['Time_of_Day',n]) = time_of_day;
    DATA.(['commanded_attack_sequence',n]) = commanded_attacks;
    save(['Cyber_Attack',n,'MW_BatLogN4.mat'],'DATA','-v7.3')

end

