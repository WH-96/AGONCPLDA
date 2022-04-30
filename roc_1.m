function [auc] = roc_1(deci,label_y,colour)
%deci 5折交叉验证得到的关联矩阵
%label_y 初始关联矩阵
[threshold,ind] = sort(deci,'descend');
roc_y = label_y(ind);

% T1=roc_y == 0;
% T_T1=find(T1~=1);
% T2=cumsum(roc_y == 0); %负样本被划分为正样本的个数  /把所有个体都预测为有关系
% T3=sum(roc_y == 0) ;   %真实负样本个数
stack_x = cumsum(roc_y == 0)/sum(roc_y == 0); %FPR

% T4=roc_y == 1;
% T5=cumsum(roc_y == 1); %真样本被划分为真的个数
% T6=sum(roc_y == 1);      %真样本个数
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);  %TPR


% x1=stack_x(2:length(roc_y));
% x2=stack_x(1:length(roc_y)-1);
% date_x=x1-x2;
% date_y=stack_y(2:length(roc_y));

auc=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));
% plot(stack_x,stack_y,colour);
% auc
% xlabel('False Positive Rate');
% ylabel('True Positive Rate');
end
