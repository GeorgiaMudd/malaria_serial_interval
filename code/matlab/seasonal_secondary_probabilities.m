function seasonal_secondary_probabilities(file_number, amp, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%ADJUSTABLE PARAMETERS

%p = 0.88;                                                                                      %Mean Probabilty of mosquito survival per day
%g = 1 - p;                                                                                     %Per-Capita Death Rate of Mosquitoes
T = @(x,p) 5*sin(2*pi*(x+p-1)/365)+ 25;                                                           %Temperature as a function of time and start date. 
g = @(x,p) (-log(-0.000828*((T(x,p))^2)+0.0367*T(x,p) + 0.522));                                              %Mosquito Survivorship Function, function of temperature.
incubation_max = 18;                                                                           %Max incubation in the mosquito (18)
emersion_constant = 6;                                                                         %Constant measuring the time until emersion of gametocytes from the liver.
%max_lifespan = ceil((log(10^-2)/(-1*g)));                                                     %max_lifespan = 28
max_lifespan = ceil((log(10^-2)/(-1*0.1303)));
normal_mosquitodeath = 1/(0.9949);                                                             %Normalizing constant for function mosquito_death
normal_incubationprobability = 1/(normcdf(18,10.35,2.472235)-normcdf(0,10.35,2.472235));         %Normalizing constant for function incubation_probability
max_density = max(normpdf(1:365,180,sigma));
mosquito_density = @(x,p) (amp*normpdf(mod(x+p-1,365),180,sigma)/max_density + 1);                                                 %Mosquito Density Function; function of time and start date. 
% Read in one column into array
human_mosquitoinf = xlsread('Daily_Infectivities.xlsx','Sheet4', 'A1001:ADU1001');                                      %Reads average daily human to mosquito infectivities into a column vector
len_data = size(human_mosquitoinf,2);
human_mosquitoinfmatrix = zeros(1, len_data);

    for j = 1:len_data
        human_mosquitoinfmatrix(1,j) = human_mosquitoinf(j)*mosquito_density(j,file_number);
    end
    human_mosquitoinfmatrix(1,:) = human_mosquitoinfmatrix(1,:)/sum(human_mosquitoinfmatrix(1,:));
                                                 
secondary_probabilities = zeros(len_data + incubation_max + max_lifespan+emersion_constant,1); %Vector. Will contain mapped probabilities for secondary infection. 
display(size(secondary_probabilities,1));
%%FUNCTIONS
mosquito_death = @(x,p) normal_mosquitodeath * ((1 - exp(-g(x,p) * (x))) - (1 - exp(-g(x,p) * (x-1))));    %Anonymous function. Represents probability of transmission from mosquito to human. 
incubation_probability = @(x) normal_incubationprobability*((normcdf(x,10.35, 2.472235))-normcdf(x-1,10.35,2.472235));%Anonymous function. Represents probability of incubation period on any day X. Taken from norm(5,2)

%Loop. Maps mosquito-to-human transmission probabilities onto
%human-to-mosquito transmission probabilities.

    for i = 1:len_data
        for j = 0:incubation_max                                           %normcdf(18,11.15,2.472235) = 0.9972;
            transition_value = human_mosquitoinfmatrix(1,i)*incubation_probability(j);
            for k = 1:max_lifespan                                                                           
                secondary_probabilities(i+j+k+emersion_constant) = secondary_probabilities(i+j+k+emersion_constant) + transition_value * mosquito_death(k,j);
            end
        end
    end

csvwrite(strcat(int2str(amp),'_',int2str(sigma), 'seasonal',int2str(file_number),'.csv'),secondary_probabilities);

%the time between bite and exiting liver stage(6 days)
%and the liver and fever (in geoff's model) and then the poisson (fever to
%clinical presentation)
% maxtreatment_wait = 14;
% feverday_probability = csvread('feverday_probabilities.csv');              %Read in file with probabilities of time between emergence and fever
% display(size(feverday_probability,1));
% bite_to_fever = zeros(size(feverday_probability,1) + emersion_constant,1);   %Create vector of zeros with size corresponding to the length of the imported file + emersion constant
% for i = 1:size(feverday_probability)                                       %Adding feverday probabilities into the vector of zeros. Offset by emersion constant. 
%     bite_to_fever(i+emersion_constant) = feverday_probability(i);
% end
% 
% fever_to_clinic = @(x) poisscdf(x,3.068807)-poisscdf(x-1,3.068807);        %Anonymous function returning the probability of the time between fever and clinic. Range is 0:14
% bite_to_clinic = zeros(size(bite_to_fever,1) + maxtreatment_wait,1);
% generation_time = zeros(size(secondary_probabilities,1) + size(bite_to_clinic,1),1);
% display(size(bite_to_clinic,1));
% 
% for i = 1:size(bite_to_fever,1)                                              %Convolving to create vector of noise to be added/subtracted to the generation time to yield the serial interval
%     for j = 0:maxtreatment_wait
%         bite_to_clinic(i+j) = bite_to_clinic(i+j) + bite_to_fever(i)*fever_to_clinic(j);
%     end
% end
% for i = 1:size(secondary_probabilities)
%     for j=1:(size(bite_to_clinic,1))
%         generation_time(i+j)=generation_time(i+j) + secondary_probabilities(i)*bite_to_clinic(j);               %Raises question of day 1. Ask Dr. Perkins
%     end
% end
% %csvwrite('gt_plus_noise.csv', generation_time);
% serial_interval = zeros(size(generation_time,1)+abs(1 - size(bite_to_clinic,1)),1);
% for i = 1:size(generation_time,1)
%     for j = 1:size(bite_to_clinic)
%         serial_interval(abs(1-size(bite_to_clinic,1))+1 + i - j) = serial_interval(abs(1-size(bite_to_clinic,1))+1+i-j) + generation_time(i)*bite_to_clinic(j);
%     end
% end
% 
% filename_final = 'new_drug_serialinterval_general.csv';
% csvwrite(filename_final, serial_interval);   
end
