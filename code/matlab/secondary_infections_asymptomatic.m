function secondary_infections_asymptomatic(filename,numberofgenerations)
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ADJUSTABLE PARAMETERS
% p = 0.88;                                                                                      %Probabilty of mosquito survival per day in Namawala, Tanzania (0.83)
% g = 1 - p;                                                                                     %Per-Capita Death Rate of Mosquitoes
% incubation_max = 18;                                                                           %Max incubation in the mosquito (18)
emersion_constant = 6;                                                                         %Constant measuring the time until emersion of gametocytes from the liver.
% max_lifespan = ceil((log(10^-2)/(-1*g)));                                                      %max_lifespan = 28
% normal_mosquitodeath = 1/(0.9949);                                                             %Normalizing constant for function mosquito_death
% normal_incubationprobability = 1/(normcdf(18,10.1,2.472235)-normcdf(0,10.1,2.472235));                                 %Normalizing constant for function incubation_probability
% Read in one column into array
% human_mosquitoinf = xlsread(filename);                                                          %Reads average daily human to mosquito infectivities into a column vector
% 
% len_data = size(human_mosquitoinf,1);                                                     %Number of values in data. 
% secondary_probabilities = zeros(len_data + incubation_max + max_lifespan+emersion_constant, 1); %Vector. Will contain mapped probabilities for secondary infection. 
% display(len_data);
% display(max_lifespan);
% display(incubation_max);
% %%FUNCTIONS
% mosquito_death = @(x) normal_mosquitodeath * ((1 - exp(-g * (x))) - (1 - exp(-g * (x-1))));    %Anonymous function. Represents probability of transmission from mosquito to human. 
% incubation_probability = @(x) normal_incubationprobability*((normcdf(x,10.1, 2.472235))-normcdf(x-1,10.1,2.472235));%Anonymous function. Represents probability of incubation period on any day X. Taken from norm(5,2)
% %Loop. Maps mosquito-to-human transmission probabilities onto
% %human-to-mosquito transmission probabilities.
% for i = 1:len_data
%     for j = 1:max_lifespan
%         transition_value = human_mosquitoinf(i)*(mosquito_death(j));
%         for k = 0:incubation_max                                                                           %normcdf(18,11.15,2.472235) = 0.9972;
%             secondary_probabilities(i+j+k+emersion_constant) = secondary_probabilities(i+j+k+emersion_constant) + transition_value * incubation_probability(k);
%         end
%     end
% end
% %% Adding in Detection of Asymptomatic Individual %% 
% size_nb = 2.625162;                                                           %Estimated in R
% mu = 41.42700;                                                             %Estimated in R
% 
% p  = size_nb/(size_nb + mu);                                                     %Probability of Success in Negative Binomial
% detection_probability = @(x) nbinpdf(x, size_nb, p);                          %Probability of Detection Based on Modeling of Gametocyte Densities (>20) 
% max_detection = 150;                                                       %Max number of days to detect
% 
% gt_plus_noise = zeros(size(secondary_probabilities,1) + max_detection,1);
% 
% for i = 1:len_data
%     for j = 1:max_detection
%         gt_plus_noise(i+j) = gt_plus_noise(i+j) + secondary_probabilities(i)*detection_probability(j);
%     end
% end
% 
% asymptomatic_gt = zeros(size(gt_plus_noise,1)+abs(1-max_detection),1);
% for i = 1:size(gt_plus_noise,1)
%     for j = 1:max_detection
%         asymptomatic_gt(abs(1-max_detection)+1 + i - j) = asymptomatic_gt(abs(1-max_detection)+1+i-j) + gt_plus_noise(i)*detection_probability(j);
%     end
% end
% csvwrite('asymptomatic_generationtime.csv', asymptomatic_gt);  
%         
%% Serial Interval (Treated -> Asymptomatic)
generation = csvread(filename);
generation = generation/sum(generation);
generation_size = size(generation,1);
maxtreatment_wait = 14;
feverday_probability = csvread('feverday_probabilities.csv');
display(size(feverday_probability,1));
bite_to_fever = zeros(size(feverday_probability,1)+emersion_constant,1);
for i = 1:size(feverday_probability)
    bite_to_fever(i+emersion_constant) = feverday_probability(i);
end

fever_to_clinic = @(x) poisscdf(x,3.068807)-poisscdf(x-1,3.068807);        %Anonymous function returning the probability of the time between fever and clinic. Range is 0:14
bite_to_clinic = zeros(size(bite_to_fever,1) + maxtreatment_wait,1);
display(size(bite_to_clinic,1));
for i = 1:size(bite_to_fever,1)                                              %Convolving to create vector of noise to be added/subtracted to the generation time to yield the serial interval
    for j = 0:maxtreatment_wait
        bite_to_clinic(i+j) = bite_to_clinic(i+j) + bite_to_fever(i)*fever_to_clinic(j);
    end
end


%Adding Asymptomatic Component to Generation Time 
gt_plus_asymptomatic = zeros(generation_size+ max_detection,1);
for i  = 1:generation_size
    for j = 1:max_detection
        gt_plus_asymptomatic(i+j) = gt_plus_asymptomatic(i+j) + generation(i)*detection_probability(j);
    end
end

serialinterval_treatedtoasymptomatic = zeros(size(gt_plus_asymptomatic,1)+abs(1-size(bite_to_clinic,1)),1);

for i = 1:size(gt_plus_asymptomatic,1)
    for j = 1:size(bite_to_clinic)
        serialinterval_treatedtoasymptomatic(abs(1-size(bite_to_clinic,1))+1+i-j) = serialinterval_treatedtoasymptomatic(abs(1-size(bite_to_clinic,1))+1+i-j) + gt_plus_asymptomatic(i)*bite_to_clinic(j);
    end
end
filename_final = strcat('serialintervalgen', int2str(numberofgenerations), '_untreatedtoasymptomatic.csv');
csvwrite(filename_final, serialinterval_treatedtoasymptomatic);


SERIAL INTERVAL (Untreated -> Symptomatic)
gt_plus_symptomatic = zeros(generation_size+size(bite_to_clinic,1),1);
for i = 1:generation_size
    for j = 1:size(bite_to_clinic)
        gt_plus_symptomatic(i+j) = gt_plus_symptomatic(i+j) + generation(i)*bite_to_clinic(j);
    end
end

serialinterval_untreatedtosymptomatic = zeros(size(gt_plus_symptomatic,1)+abs(1-max_detection),1);
for i = 1:size(gt_plus_symptomatic,1)
    for j = 1:max_detection
        serialinterval_untreatedtosymptomatic(abs(1-max_detection)+1+i-j) = serialinterval_untreatedtosymptomatic(abs(1-max_detection)+1+i-j) + gt_plus_symptomatic(i)*detection_probability(j);
    end
end
filename_final = strcat('serialintervalgen', int2str(numberofgenerations), 'untreatedtosymptomatic.csv');
csvwrite(filename_final, serialinterval_untreatedtosymptomatic);

%% SERIAL INTERVAL (UNTREATED -> ASYMPTOMATIC)
gt_plus_asymptomatic = zeros(generation_size+max_detection,1);
for i = 1:generation_size
    for j = 1:max_detection
        gt_plus_asymptomatic(i+j) = gt_plus_asymptomatic(i+j) + generation(i)*detection_probability(j);
    end
end

serialinterval_untreatedtoasymptomatic = zeros(size(gt_plus_asymptomatic,1)+abs(1-max_detection),1);
for i = 1:size(gt_plus_asymptomatic,1)
    for j = 1:max_detection
        serialinterval_untreatedtoasymptomatic(abs(1-max_detection)+1+i-j) = serialinterval_untreatedtoasymptomatic(abs(1-max_detection)+1+i-j) + gt_plus_asymptomatic(i)*detection_probability(j);
    end
end
filename_final = strcat('serialintervalgen', int2str(numberofgenerations), '_untreatedtoasymptomatic.csv');
csvwrite(filename_final, serialinterval_untreatedtoasymptomatic);
    
    