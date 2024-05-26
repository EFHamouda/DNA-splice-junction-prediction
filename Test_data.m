%  Hybrid ensemble model for DNA splice junction prediction version 1.0                                               %
% Main paper: "A hybird approach of ensemble learning and gray wolf optimzation for DNA splice junction prediction "  %
%                                                                                                                     %
%  Eslam Hamouda and Mayada Tarek                                                                                     %
%_____________________________________________________________________________________________________________________%

function  [Accurcy,TP_N,TP_IE,TP_EI,Precision_N,Precision_IE,Precision_EI,F_score_N,F_score_IE,F_score_EI, MCC_N,MCC_IE,MCC_EI, AUC_N,AUC_IE,AUC_EI] = Test_data(Destination_position_All,Test_set,modelAll)

numberofmodels=size(Destination_position_All,1);
predictionsAll=[];
scores=[];
for i= 1:numberofmodels

NewAgent=Map_input(Destination_position_All(i,:));
selectedFeatures = find(NewAgent ==1 );
selectedFeatures=[selectedFeatures,181]; %the class label

Current_Test_set=Test_set(:,selectedFeatures);
model = modelAll{i};


[predictions, scores] = predict(model, Current_Test_set(:, 1:end-1));





%predictions=str2num(cell2mat(predictions));  % use only in case of Naive Classifer
 
predictionsAll=[predictionsAll predictions];

end

Finalpredictions=mode(predictionsAll')';
Accurcy=numel(find(Finalpredictions==Test_set(:,end)))/(size(Test_set,1));
Accurcy=Accurcy*100;

% EI--> 1 , IE--> 2, N-->3


confMat = confusionmat(Test_set(:,end) ,Finalpredictions);

t = ["EI", "IE","N"];
v = [1,2,3]; 
New= string(categorical(Test_set(:,end),v,t));


if size(New,1)~=1

figure
rocObj = rocmetrics(New,scores,["EI";"IE";"N"]);
plot(rocObj,"ShowDiagonalLine",false,"ShowModelOperatingPoint",false);

pp=gcf;

saveas(pp,'ROC_Fig.fig');

figure
FinalpredictionsNewlables=string(categorical(Finalpredictions,v,t));
cm=confusionchart(New,FinalpredictionsNewlables,'Normalization','column-normalized');

saveas(cm,'Confusion_Matraix_fig.fig');

end

if size(confMat,1)==3
TP_EI=confMat(1,1);
TP_IE=confMat(2,2);
TP_N=confMat(3,3);

FP_EI=sum(confMat(1, :)) - confMat(1,1);
FP_IE=sum(confMat(2, :)) - confMat(2,2);
FP_N=sum(confMat(3, :)) - confMat(3,3);

FN_EI=sum(confMat(:, 1)) - confMat(1,1); 
FN_IE=sum(confMat(:, 2)) - confMat(2,2); 
FN_N=sum(confMat(:, 3)) - confMat(3,3); 

TN_EI=sum(confMat(:)) - sum(confMat(1,:)) - sum(confMat(:,1)) + confMat(1,1);
TN_IE=sum(confMat(:)) - sum(confMat(2,:)) - sum(confMat(:,2)) + confMat(2,2);
TN_N=sum(confMat(:)) - sum(confMat(3,:)) - sum(confMat(:,3)) + confMat(3,3);


TP_EI=TP_EI/sum(Test_set(:,end)==1);
TP_IE=TP_IE/sum(Test_set(:,end)==2);
TP_N=TP_N/sum(Test_set(:,end)==3);

FP_EI=FP_EI/sum(Test_set(:,end)==1);
FP_IE=FP_IE/sum(Test_set(:,end)==2);
FP_N=FP_N/sum(Test_set(:,end)==3);

FN_EI=FN_EI/sum(Test_set(:,end)==1); 
FN_IE=FN_IE/sum(Test_set(:,end)==2); 
FN_N=FN_N/sum(Test_set(:,end)==3); 

TN_EI=TN_EI/sum(Test_set(:,end)==1);
TN_IE=TN_IE/sum(Test_set(:,end)==2);
TN_N=TN_N/sum(Test_set(:,end)==3);


Precision_N=TP_N/(TP_N+FP_N);
Precision_IE=TP_IE/(TP_IE+FP_IE);
Precision_EI=TP_EI/(TP_EI+FP_EI);


Recal_N=TP_N/(TP_N+FN_N);
Recal_IE=TP_IE/(TP_IE+FN_IE);
Recal_EI=TP_EI/(TP_EI+FN_EI);

F_score_N=2*(Precision_N*Recal_N)/(Precision_N+Recal_N);
F_score_IE=2*(Precision_IE*Recal_IE)/(Precision_IE+Recal_IE);
F_score_EI=2*(Precision_EI*Recal_EI)/(Precision_EI+Recal_EI);

AUC_N=0.5*((TP_N/(TP_N+FN_N))+(TN_N/(TN_N+FP_N)));
AUC_IE=0.5*((TP_IE/(TP_IE+FN_IE))+(TN_IE/(TN_IE+FP_IE)));
AUC_EI=0.5*((TP_EI/(TP_EI+FN_EI))+(TN_EI/(TN_EI+FP_EI)));

MCC_N=(TP_N*TN_N-FP_N*FN_N)/sqrt((TP_N+FP_N)*(TP_N+FN_N)*(TN_N+FP_N)*(TN_N+FN_N));
MCC_IE=(TP_IE*TN_IE-FP_IE*FN_IE)/sqrt((TP_IE+FP_IE)*(TP_IE+FN_IE)*(TN_IE+FP_IE)*(TN_IE+FN_IE));
MCC_EI=(TP_EI*TN_EI-FP_EI*FN_EI)/sqrt((TP_EI+FP_EI)*(TP_EI+FN_EI)*(TN_EI+FP_EI)*(TN_EI+FN_EI));
else
   TP_N=0;
   TP_IE=0;
   TP_EI=0;
   Precision_N=0;
   Precision_IE=0;
   Precision_EI=0;
   F_score_N=0;
   F_score_IE=0;
   F_score_EI=0;
   MCC_N=0;
   MCC_IE=0;
   MCC_EI=0;
   AUC_N=0;
   AUC_IE=0;
   AUC_EI =0;
end

end 