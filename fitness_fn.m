%  Hybrid ensemble model for DNA splice junction prediction version 1.0                                               %
% Main paper: "A hybird approach of ensemble learning and gray wolf optimzation for DNA splice junction prediction "  %
%                                                                                                                     %
%  Eslam Hamouda and Mayada Tarek                                                                                     %
%                                                                                                                     %
%_____________________________________________________________________________________________________________________%

function [Error, model]= fitness_fn(Agent,or_Train_Validation_data)


 NewAgent=Map_input(Agent);
 
selectedFeatures = find(NewAgent ==1 );
selectedFeatures=[selectedFeatures,181]; %the class label

or_Train_Validation_data=or_Train_Validation_data(:,selectedFeatures);

index=floor(size(or_Train_Validation_data,1))/4;
index=floor(index);
valdition=or_Train_Validation_data(1:index,:);
trainn=or_Train_Validation_data(index+1:end,:);

% SVM
t = templateSVM('SaveSupportVectors',true);
model = fitcecoc(trainn(:, 1:end-1), trainn(:, end) ,'Learners',t);

% DT
%model = fitctree(trainn(:, 1:end-1), trainn(:, end) );

% KNN 
%model = fitcknn(or_Train_Validation_data(:, 1:end-1), or_Train_Validation_data(:, end), 'NumNeighbors', 5);

%Naice base clasisfer
%model = fitcnb(or_Train_Validation_data(:, 1:end-1), or_Train_Validation_data(:, end),'ClassNames',{'1','2','3'},'DistributionNames','mn');



predictions = predict(model,valdition(:, 1:end-1)); 



%predictions=str2num(cell2mat(predictions));                               % used only for Naive classifer
 
Error=numel(find(predictions~=valdition(:,end)))/(size(valdition,1));

end