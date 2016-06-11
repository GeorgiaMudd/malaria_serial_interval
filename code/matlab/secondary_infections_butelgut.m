function secondary_infections_butelgut(filename)
%UNTITLED Summary of this function goes here
%This function serves to map a series of probability distributions upon
%eachother. 
%ADJUSTABLE PARAMETERS

p = 0.86;                                                                                      %Probabilty of mosquito survival per day in Butelgut, Papa New Guinea
g = 1 - p;                                                                                     %Per-Capita Death Rate of Mosquitoes
incubation_max = 18;                                                                           %Max Incubation in the Mosquito
max_lifespan = ceil((log(10^-2)/(-1*g)));                                                      %max_lifespan = 33
emersion_constant = 6;                                                                         %Constant measuring the time until emersion of gametocytes from the liver.
normal_mosquitodeath = 1/(0.99015);                                                             %Normalizing constant for function mosquito_death
normal_incubationprobability = 1/(normcdf(18,8.95,2.472235)-normcdf(0,8.95,2.472235));                              %Normalizing constant for function incubation_probability
% Read in one column into array
human_mosquitoinf = xlsread(filename);                                        %Reads average daily human to mosquito infectivities into a column vector

len_data = size(human_mosquitoinf,2);                                                          %Number of values in data. 
secondary_probabilities = zeros(len_data + incubation_max + max_lifespan + emersion_constant, 1);%Vector. Will contain mapped probabilities for secondary infection. 
%%FUNCTIONS
mosquito_death = @(x) normal_mosquitodeath * ((1 - exp(-g * (x))) - (1 - exp(-g * (x-1))));    %Anonymous function. Represents probability of transmission from mosquito to human. 
incubation_probability = @(x) normal_incubationprobability*((normcdf(x,8.95,2.472235))-normcdf(x-1,8.95,2.472235));%Anonymous function. Represents probability of incubation period in mosquito on any day X. Taken from norm(5,2)
%Loop. Maps mosquito-to-human transmission probabilities onto
%human-to-mosquito transmission probabilities.
for i = 1:len_data
    for j = 1:max_lifespan
        transition_value = human_mosquitoinf(i)*(mosquito_death(j));
        for k = 0:incubation_max                                                               %normcdf(18,8.95,2.472235) = 0.9999;
            secondary_probabilities(i+j+k + emersion_constant) = secondary_probabilities(i+j+k + emersion_constant) + transition_value * incubation_probability(k);
        end
    end
end
%the time between bite and exiting liver stage(6 days)
%and the liver and fever (in geoff's model) and then the poisson (fever to
%clinical presentation)
maxtreatment_wait = 14;
feverday_probability = csvread('feverday_probabilities.csv');              %Read in file with probabilities of time between emergence and fever
display(size(feverday_probability,1));
bite_to_fever = zeros(size(feverday_probability,1) + emersion_constant,1);   %Create vector of zeros with size corresponding to the length of the imported file + emersion constant
for i = 1:size(feverday_probability)                                       %Adding feverday probabilities into the vector of zeros. Offset by emersion constant. 
    bite_to_fever(i+emersion_constant) = feverday_probability(i);
end

fever_to_clinic = @(x) poisscdf(x,3.068807)-poisscdf(x-1,3.068807);        %Anonymous function returning the probability of the time between fever and clinic. Range is 0:14
transition_vector = zeros(size(bite_to_fever,1) + maxtreatment_wait,1);
serial_interval = zeros(size(secondary_probabilities,1) + size(transition_vector,1),1);
for i = 1:size(secondary_probabilities)                                    %Filling new vector
    serial_interval(i) = secondary_probabilities(i);
end

for i = 1:size(bite_to_fever)                                              %Convolving to create vector of noise to be added/subtracted to the generation time to yield the serial interval
    for j = 0:maxtreatment_wait
        transition_vector(i+j) = transition_vector(i+j) + bite_to_fever(i)*fever_to_clinic(j);
    end
end
for i = 1:size(secondary_probabilities)
    for j=1:(size(transition_vector))
        serial_interval(i+j)=serial_interval(i+j)+transition_vector(j);               %Raises question of day 1. Ask Dr. Perkins
    end
    for j=1:(size(transition_vector))
        serial_interval(i+j) = serial_interval(i+j) - transition_vector(j);
    end
end

filename_final = 'new_drug_serialinterval_butelgut.csv';
csvwrite(filename_final, serial_interval);
   
end


