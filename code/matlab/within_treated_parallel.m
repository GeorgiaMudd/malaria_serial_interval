
p = 0.88;
mean_incubation = 10.1;
total_probabilities = poisscdf(14,3.068807);
treatment_wait_probabilities = poisspdf(0:14, 3.068807);
treatment_wait_probabilities = treatment_wait_probabilities/total_probabilities;
human_to_mosquito_infectivity = zeros([1000 801]);
for i = 1:2
    treatment_wait =  randsample(0:14, 1, true, treatment_wait_probabilities);
    final_model(1,treatment_wait);
    load('new_data_1.mat');
    for j = 1:801
        human_to_mosquito_infectivity(i,j) = c(1,j);
    end
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
secondary_probabilities = zeros([2 (size(human_to_mosquito_infectivity,2)+ incubation_max + max_lifespan+emergence_constant)]);
for i = 1:2
    for j = 1:size(human_to_mosquito_infectivity,2)
        for k = 1:max_lifespan
            transition_value = human_to_mosquito_infectivity(j)*(mosquito_death(k));
            for l = 0:incubation_max
                secondary_probabilities(i,j+k+l+emergence_constant) = secondary_probabilities(i,j+k+l+emergence_constant) + transition_value * incubation_probability(l);
            end
        end
    end
end   
csvwrite('within_treated_simulation.csv',secondary_probabilities);
        