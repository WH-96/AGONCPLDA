%Artificial Gorilla Troops Optimizer

% Read the following publication first and cite if you use it

% @article{abdollahzadeh2021african,
%   title={Artificial Gorilla Troops Optimizer: A New Nature-Inspired Metaheuristic Algorithm for Global Optimization Problems},
%   author={Abdollahzadeh, Benyamin and Gharehchopogh, Farhad Soleimanian and Mirjalili, Seyedali},
%   journal={International Journal of Intelligent Systems},
%   pages={},
%   year={2021},
%   publisher={Willy},
%   url = {https://doi.org/10.1002/int.22535}
% }


clear all 
close all
clc

% Population size and stoppoing condition 



%%
SearchAgents_no=40; % Number of search agents

Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=500; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lower_bound,upper_bound,variables_no,fobj]=NCPHMDA5KCV(Function_name);      
[Silverback_Score,Silverback,convergence_curve]=GTO(SearchAgents_no,Max_iteration,lower_bound,upper_bound,variables_no,fobj);


figure 

%Draw objective space
subplot(1,2,2);
semilogy(convergence_curve,'Color','r')
title('Convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid off
box on
legend('GTO-D-L')
display(['The best solution obtained by GTO is : ', num2str(Silverback)]);
display(['The best optimal value of the objective funciton found by GTO is : ', num2str(Silverback_Score)]);
