%  Hybrid ensemble model for DNA splice junction prediction version 1.0                                               %
% Main paper: "A hybird approach of ensemble learning and gray wolf optimzation for DNA splice junction prediction "  %
%                                                                                                                     %
%  Eslam Hamouda and Mayada Tarek                                                                                     %
%_____________________________________________________________________________________________________________________%
   

function Run

clear
clc
n_runs=10;
Pop=20;
Iterations=80;
un_or_all_data=load_data();
idx = randperm(size(un_or_all_data,1)) ;
or_all_data=un_or_all_data;
or_all_data=or_all_data(idx,:);

 % the last colum is calss label, it is not ant attribute,  (so subtract 1 form feature size) 
Feature_size=size(or_all_data,2)-1;


Destination_fitness_All=[];
Convergence_All=[];
Destination_position_All=[];

Accurcy_all=[];

PrecisionN_all=[];
PrecisionIE_all=[];
PrecisionEI_all=[];


F_scoreN_all=[];
F_scoreIE_all=[];
F_scoreEI_all=[];


MCCN_All=[];
MCCIE_All=[];
MCCEI_All=[];


TPN_all=[];
TPIE_all=[];
TPEI_all=[];

AUCN_all=[];
AUCIE_all=[];
AUCEI_all=[];

modelAll=[];
tic
for t=1:n_runs   % for parallel exexcution  use  parfor 
 
% fold size
step=round(size(or_all_data,1)/10); 

for j=1:step:size(or_all_data,1)
Trainfolds=or_all_data;
if j+step-1<=size(or_all_data,1)
Test_set=or_all_data(j:j+step-1,:);
Trainfolds(j:j+step-1,:)=[];
else
   Test_set=or_all_data(j:end,:);
   Trainfolds(j:end,:)=[];
end



% 70% 30% expirment
% index=floor(size(or_all_data,1))*0.7;
% index=floor(index);
% Trainfolds=or_all_data(1:index,:);
% Test_set=or_all_data(index+1:end,:);


Destination_fitness_All=[];
Convergence_All=[];
Destination_position_All=[];
modelAll={};
k=1;

%foldcross valdition
step=round(size(Trainfolds,1)/9);

for i=1:step:size(Trainfolds,1)

Train_data=Trainfolds;

if i+step-1<=size(Trainfolds,1)
Train_data(i:i+step-1,:)=[];
else
    Train_data(i:end,:)=[];
end


[Destination_position,Destination_fitness,Convergence,model]=GWO(Pop,Iterations,Feature_size,Train_data);

Convergence_All=[Convergence_All;Convergence];
Destination_fitness_All=[Destination_fitness_All;Destination_fitness];
Destination_position_All=[Destination_position_All;Destination_position];
modelAll{k}=model;
k=k+1;
end

[Accurcy,TP_N,TP_IE,TP_EI,Precision_N,Precision_IE,Precision_EI,F_score_N,F_score_IE,F_score_EI, MCC_N,MCC_IE,MCC_EI, AUC_N,AUC_IE,AUC_EI] =Test_data(Destination_position_All,Test_set,modelAll);

Accurcy_all=[Accurcy_all;Accurcy];

if (TP_N~=0)
PrecisionN_all=[PrecisionN_all;Precision_N];
PrecisionIE_all=[PrecisionIE_all;Precision_IE];
PrecisionEI_all=[PrecisionEI_all;Precision_EI];


F_scoreN_all=[F_scoreN_all;F_score_N];
F_scoreIE_all=[F_scoreIE_all;F_score_IE];
F_scoreEI_all=[F_scoreEI_all;F_score_EI];


MCCN_All=[MCCN_All;MCC_N];
MCCIE_All=[MCCIE_All;MCC_IE];
MCCEI_All=[MCCEI_All;MCC_EI];


TPN_all=[TPN_all;TP_N];
TPIE_all=[TPIE_all;TP_IE];
TPEI_all=[TPEI_all;TP_EI];


AUCN_all=[AUCN_all;AUC_N];
AUCIE_all=[AUCIE_all;AUC_IE];
AUCEI_all=[AUCEI_all;AUC_EI];

end


end


mean_Time = sprintf('computional Time = %.2f  seconds\r\n',mean(toc));  
mean_Accury = sprintf('Accurecy = %.4f \r\n',mean(Accurcy_all));

mean_TP_N = sprintf('TP for N = %.4f \r\n',mean(TPN_all));
mean_TP_IE = sprintf('TP for IE = %.4f \r\n',mean(TPIE_all));
mean_TP_EI = sprintf('TP for EI = %.4f \r\n',mean(TPEI_all));

mean_Precision_N = sprintf('Precision for N = %.4f \r\n',mean(PrecisionN_all));
mean_Precision_IE = sprintf('Precision for IE = %.4f \r\n',mean(PrecisionIE_all));
mean_Precision_EI = sprintf('Precision for EI = %.4f \r\n',mean(PrecisionEI_all));

mean_F_score_N = sprintf('F-measure for N = %.4f \r\n',mean(F_scoreN_all));
mean_F_score_IE = sprintf('F-measure for IE = %.4f \r\n',mean(F_scoreIE_all));
mean_F_score_EI = sprintf('F-measure for EI = %.4f \r\n',mean(F_scoreEI_all));

mean_MCC_N = sprintf('MCC for N = %.4f \r\n',mean(MCCN_All));
mean_MCC_IE = sprintf('MCC for IE = %.4f \r\n',mean(MCCIE_All));
mean_MCC_EI = sprintf('MCC for EI = %.4f \r\n',mean(MCCEI_All));

mean_AUC_N= sprintf('AUC for N= %.4f \r\n',mean(AUCN_all));
mean_AUC_IE = sprintf('AUC for IE= %.4f \r\n',mean(AUCIE_all));
mean_AUC_EI= sprintf('AUC for EI= %.4f \r\n',mean(AUCEI_all));


mean_Fitness = sprintf('Fittness = %.4f \r\n',mean(Destination_fitness_All));

disp(mean_Accury);

disp(mean_TP_N);
disp(mean_TP_IE);
disp(mean_TP_EI);

disp(mean_Precision_N);
disp(mean_Precision_IE);
disp(mean_Precision_EI);

disp(mean_F_score_N);
disp(mean_F_score_IE);
disp(mean_F_score_EI);

disp(mean_MCC_N);
disp(mean_MCC_IE);
disp(mean_MCC_EI);

disp(mean_AUC_N);
disp(mean_AUC_IE);
disp(mean_AUC_EI);

disp(mean_Fitness);
disp(mean_Time);
disp('Done---------------------------');


file_name=strcat('GWO_','DNA','.txt');
fileID = fopen(file_name,'w');
fprintf(fileID,mean_Accury);

fprintf(fileID,mean_TP_N);
fprintf(fileID,mean_TP_IE);
fprintf(fileID,mean_TP_EI);

fprintf(fileID,mean_Precision_N);
fprintf(fileID,mean_Precision_IE);
fprintf(fileID,mean_Precision_EI);

fprintf(fileID,mean_F_score_N);
fprintf(fileID,mean_F_score_IE);
fprintf(fileID,mean_F_score_EI);

fprintf(fileID,mean_MCC_N);
fprintf(fileID,mean_MCC_IE);
fprintf(fileID,mean_MCC_EI);

fprintf(fileID,mean_AUC_N);
fprintf(fileID,mean_AUC_IE);
fprintf(fileID,mean_AUC_EI);

fprintf(fileID,mean_Fitness);
fprintf(fileID,mean_Time);
fprintf(fileID,'-------------------------\r\n');
fclose(fileID);

save(strcat('All_Convergence_','DNA'),'Convergence_All');
Destination_position_All=Destination_position_All>=0.5;
save(strcat('All_Destination_position_','DNA'),'Destination_position_All');


end