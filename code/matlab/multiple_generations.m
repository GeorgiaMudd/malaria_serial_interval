function multiple_generations(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

generationtime_individual = xlsread(filename);                             %Reads in file for the average serial interval distribution of an untreated case.
generationtime_individual = generationtime_individual/sum(generationtime_individual);   %Normalizing the distribution
len_data = size(generationtime_individual, 1);                             %Accesses the length of the serial interval distribution (853). 
display(len_data);
%multiple_generation2 = zeros(len_data + len_data,1);                       %Creates vector of zeros that will ultimately be the distribution of the generation time across 2 transmission cycles. 
%multiple_generation3 = zeros(len_data + len_data + len_data, 1);           %Creates vector of zeros that will ultimately be the distribution of the generation time across 3 transmission cycles. 
multiple_generation4 = zeros(len_data + len_data + len_data + len_data,1); %Creates vector of zeros that will ultimately be the distribution of the generation time across 4 transmission cycles. 
%multiple_generation5 = zeros(len_data + len_data + len_data + len_data + len_data,1); %Creates vector of zeros that will ultimately be the distribution of the generation time across 5 transmission cycles. 

%%CONVOLUTION%%
% for i = 1:len_data
%     for j = 1:len_data
%         multiple_generation2(i + j) = multiple_generation2(i+j) + generationtime_individual(i)*generationtime_individual(j);
%     end
% end
% for i = 1:len_data
%     for j = 1:len_data
%         for k = 1:len_data
%             multiple_generation3(i+j+k) = multiple_generation3(i+j+k) + generationtime_individual(i)*generationtime_individual(j)*generationtime_individual(k);
%         end
%     end
% end

for i = 1:len_data
    for j = 1:len_data
        for k = 1:len_data
            for l = 1:len_data
                multiple_generation4(i+j+k+l) = multiple_generation4(i+j+k+l) + generationtime_individual(i)*generationtime_individual(j)*generationtime_individual(k)*generationtime_individual(l);
            end
        end
    end
end


%%WRITING FILE%%
%filename_final = 'GenerationTime_Across_2Transmissions_Treated';
%filename_final_2 = 'GenerationTime_Across_3Transmissions_Treated';
filename_final_3 = 'GenerationTime_Across_4Transmissions_Treated';
%filename_final_4 = 'GenerationTime_Across_5Transmissions_Treated';
%csvwrite(filename_final, multiple_generation2);
%csvwrite(filename_final_2, multiple_generation3);
csvwrite(filename_final_3, multiple_generation4);
%csvwrite(filename_final_4, multiple_generation5);

end

