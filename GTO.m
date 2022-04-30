
function [Silverback_Score,Silverback,convergence_curve]=GTO(pop_size,max_iter,lower_bound,upper_bound,variables_no,fobj)


% initialize Silverback
Silverback=[];
Silverback_Score=inf;

%Initialize the first random population of Gorilla
X=initialization(pop_size,variables_no,upper_bound,lower_bound);


convergence_curve=zeros(max_iter,1);
%% 挑选银背大猩猩
for i=1:pop_size   
    Pop_Fit(i)=fobj(X(i,:));%#ok
    if Pop_Fit(i)<Silverback_Score 
            Silverback_Score=Pop_Fit(i); 
            Silverback=X(i,:);
    end
end


GX=X(:,:);
lb=ones(1,variables_no).*lower_bound; 
ub=ones(1,variables_no).*upper_bound; 

%%  Controlling parameter

p=0.03;                       % 决定勘探阶段所使用的勘探方式
Beta=3;                       %               
w=0.8;                        % 决定开发阶段 是追随银背大猩猩还是 竞争成为新的银背大猩猩

%%Main loop
for It=1:max_iter 
    
    a=(cos(2*rand)+1)*(1-It/max_iter);    %eq 2
    C=a*(2*rand-1);                       % eq 4

%% Exploration:

    for i=1:pop_size
        if rand<p    
            GX(i,:) =(ub-lb)*rand+lb;              % 移动到未知的地方
        else  
            if rand>=0.5
                Z = unifrnd(-a,a,1,variables_no);  % eq 6  
                H=Z.*X(i,:);                       % eq 5
                GX(i,:)=(rand-a)*X(randi([1,pop_size]),:)+C.*H;   % 移动到其他大猩猩的位置
            else   
                GX(i,:)=X(i,:)-C.*(C*(X(i,:)- GX(randi([1,pop_size]),:))+rand*(X(i,:)-GX(randi([1,pop_size]),:)));  %在已知领域进行探索 
            end
        end
    end       
       
    GX = boundaryCheck(GX, lower_bound, upper_bound);
    
    % Group formation operation 位置好则更新，否则不更新
    for i=1:pop_size
         New_Fit= fobj(GX(i,:));          
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end
         if New_Fit<Silverback_Score 
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
         end
    end
    
%% Exploitation:  
    for i=1:pop_size
       if a>=w  
            g=2^C;
            delta= (abs(mean(GX)).^g).^(1/g);
            GX(i,:)=C*delta.*(X(i,:)-Silverback)+X(i,:);    %追随银背大猩猩
       else
           
           if rand>=0.5
              h=randn(1,variables_no);
           else
              h=randn(1,1);
           end
           r1=rand; 
           GX(i,:)= Silverback-(Silverback*(2*r1-1)-X(i,:)*(2*r1-1)).*(Beta*h); 
           
       end
    end
   
    GX = boundaryCheck(GX, lower_bound, upper_bound);
    
    % Group formation operation    
    for i=1:pop_size
         New_Fit= fobj(GX(i,:));
         if New_Fit<Pop_Fit(i)
            Pop_Fit(i)=New_Fit;
            X(i,:)=GX(i,:);
         end
         if New_Fit<Silverback_Score 
            Silverback_Score=New_Fit; 
            Silverback=GX(i,:);
         end
    end
Silverback_Score             
convergence_curve(It)=Silverback_Score;
fprintf("In Iteration %d, best estimation of the global optimum is %4.4f \n ", It,Silverback_Score );
         
end 