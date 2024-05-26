%  Hybrid ensemble model for DNA splice junction prediction version 1.0                                               %
% Main paper: "A hybird approach of ensemble learning and gray wolf optimzation for DNA splice junction prediction "  %
%                                                                                                                     %
%  Eslam Hamouda and Mayada Tarek                                                                                     %
%                                                                                                                     %
%   Based on code devloped by  S. Mirjalili, S. M. Mirjalili, A. Lewis, , Grey Wolf Optimizer (GWO)                   %
%   Advances in Engineering Software, Volume 69,2014                                                                  %
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2013.12.007.                                                          %
%_____________________________________________________________________________________________________________________%


% Grey Wolf Optimizer

function [Alpha_pos,Alpha_score,Convergence,BestSVMmodel]=GWO(SearchAgents_no,Max_iter,dim,train_data)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; 

Beta_pos=zeros(1,dim);
Beta_score=inf; 

Delta_pos=zeros(1,dim);
Delta_score=inf;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim);
Convergence=zeros(1,Max_iter);

l=0;% Loop counter
BestSVMmodel=[];
% Main loop
a=2;
% paramters used for choatic map functions
aa=0.5;
bb=0.2;

while l<Max_iter
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        
        % Calculate objective function for each search agent
        
        [fitness, model]=fitness_fn(Positions(i,:),train_data);
        
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score 
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
            BestSVMmodel=model;
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly from 2 to 0

% Value=20-l*((20-1e-10)/Max_iter);

%Chebyshev map
   % a=cos(i*acos(a));
    %G=(a+1)*Value)/2;

%Gauss/mouse map
 %   if a~=0
 %      
 %       a=mod(1/a,1);
%    end
%    a=a*Value;
  
%Sinusoidal map
  %   a = 2.3*a^2*sin(pi*a);
  %   a=a*Value;

%Circle map

 %   a=mod(a+bb-(aa/(2*pi))*sin(2*pi*a),1); 
 %   a=a*Value;

    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        

        for j=1:size(Positions,2)     
        
  
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand();
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)

            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1         

            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1

            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
         

            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand();
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
                 
 
        
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)

             
        end
    end
    l=l+1;    
        Convergence(l)=Alpha_score;

end
end


