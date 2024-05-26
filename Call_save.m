%  Hybrid ensemble model for DNA splice junction prediction version 1.0                                               %
% Main paper: "A hybird approach of ensemble learning and gray wolf optimzation for DNA splice junction prediction "  %
%                                                                                                                     %
%  Eslam Hamouda and Mayada Tarek                                                                                     %
%_____________________________________________________________________________________________________________________%

function Call_save(file_name,Destination_Pos1,Destination_Pos2,Destination_Pos3,Destination_fitness_All)
save(file_name,'Destination_Pos1','Destination_Pos2','Destination_Pos3','Destination_fitness_All');
end 