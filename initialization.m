%  Hybrid ensemble model for DNA splice junction prediction version 1.0                                               %
% Main paper: "A hybird approach of ensemble learning and gray wolf optimzation for DNA splice junction prediction "  %
%                                                                                                                     %
%  Eslam Hamouda and Mayada Tarek                                                                                     %
%_____________________________________________________________________________________________________________________%


% This function creates the first random population of agents

function X=initialization(SearchAgents_no,dim)


 X=rand(SearchAgents_no,dim);
for i=1:SearchAgents_no

 for j=1:dim
    
 if(X(i,j)>0.5)
     X(i,j)=1;
 else
     X(i,j)=0;
 end
end

end