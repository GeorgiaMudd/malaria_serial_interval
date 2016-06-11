function secondary_infections(filename, sheetname, range)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%ADJUSTABLE PARAMETERS
g = 0.05;                                                                                      %Per-Capita Death Rate of Mosquitoes
incubation_max = 11;
max_lifespan = ceil((log(10^-2)/(-1*g)));                                                      %max_lifespan = 93
normal_mosquitodeath = 1/(0.9904);                                                             %Normalizing constant for function mosquito_death
normal_incubationprobability = 1/(normcdf(9,5,2)-normcdf(0,5,2));                                 %Normalizing constant for function incubation_probability
% Read in one column into array
human_mosquitoinf = xlsread(filename,sheetname, range);                                        %Reads average daily human to mosquito infectivities into a column vector

len_data = size(human_mosquitoinf,2);                                                          %Number of values in data. 
secondary_probabilities = zeros(len_data + incubation_max + max_lifespan, 1);                                   %Vector. Will contain mapped probabilities for secondary infection. 
display(len_data);
display(max_lifespan);
display(incubation_max);
%%FUNCTIONS
mosquito_death = @(x) normal_mosquitodeath * ((1 - exp(-g * (x))) - (1 - exp(-g * (x-1))));    %Anonymous function. Represents probability of transmission from mosquito to human. 
incubation_probability = @(x) normal_incubationprobability*((normcdf(x,5,2))-normcdf(x-1,5,2));%Anonymous function. Represents probability of incubation period on any day X. Taken from norm(5,2)
%Loop. Maps mosquito-to-human transmission probabilities onto
%human-to-mosquito transmission probabilities.
for i = 1:len_data
    for j = 1:max_lifespan
        transition_value = human_mosquitoinf(i)*(mosquito_death(j));
        for k = 0:incubation_max                                                                           %normcdf(11,5,2) = 0.9987;
            secondary_probabilities(i+j+k) = secondary_probabilities(i+j+k) + transition_value * incubation_probability(k);
        end
    end
end
filename_final = 'secondary_probabilities_test.xlsx';
xlswrite(filename_final, secondary_probabilities);
    
end


