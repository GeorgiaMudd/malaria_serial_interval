function within_host_untreated(file_number)
    file_name = strcat('new_data_',int2str(file_number), '.mat');
    p = 0.88;
    mean_incubation = 10.1;

        final_models(1,file_number);
        load(file_name);
        for j = 1:801
            human_to_mosquito_infectivity(j) = c(1,j);
        end

    g = 1-p;                                                                    % Per-capita death rate of mosquitoes. Complement of the probability (p) of survival provided as an argument.
    incubation_max = 18;                                                        % Max incubation (EIP) in the mosquito.
    emergence_constant = 6;                                                     %IIP modeled as a constant six days.
    max_lifespan = ceil((log(10^-2)/(-1*g)));                                   %Calculates max lifespan of the mosquito based on the per-capita death rate.
    display(max_lifespan);
    normal_mosquitodeath = 1/normalizing(g);                                    %Calls normalizing function to calculate normalizing constant. 
    mosquito_death = @(x) normal_mosquitodeath * ((1 - exp(-g * (x))) - (1 - exp(-g * (x-1)))); %Anonymous function. Represents probability of transmission from mosquito to human. 
    normal_incubationprobability = 1/(normcdf(18,mean_incubation,2.472235)-normcdf(0,mean_incubation,2.472235)); %Normalizing constant for the probability of mosquito incubation on any given day.
    incubation_probability = @(x) normal_incubationprobability*((normcdf(x,mean_incubation,2.472235))-normcdf(x-1,mean_incubation,2.472235));

        secondary_probabilities = zeros(size(human_to_mosquito_infectivity,2)+ incubation_max + max_lifespan+emergence_constant, 1);
        for j = 1:size(human_to_mosquito_infectivity,2)
            for k = 1:max_lifespan
                transition_value = human_to_mosquito_infectivity(j)*(mosquito_death(k));
                for l = 0:incubation_max
                    secondary_probabilities(j+k+l+emergence_constant) = secondary_probabilities(j+k+l+emergence_constant) + transition_value * incubation_probability(l);
                end
            end
        end
        
    csvwrite(strcat('untreated_output',int2str(file_number),'.csv'),secondary_probabilities);
end


