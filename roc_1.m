function [auc] = roc_1(deci,label_y,colour)
%deci 5�۽�����֤�õ��Ĺ�������
%label_y ��ʼ��������
[threshold,ind] = sort(deci,'descend');
roc_y = label_y(ind);

% T1=roc_y == 0;
% T_T1=find(T1~=1);
% T2=cumsum(roc_y == 0); %������������Ϊ�������ĸ���  /�����и��嶼Ԥ��Ϊ�й�ϵ
% T3=sum(roc_y == 0) ;   %��ʵ����������
stack_x = cumsum(roc_y == 0)/sum(roc_y == 0); %FPR

% T4=roc_y == 1;
% T5=cumsum(roc_y == 1); %������������Ϊ��ĸ���
% T6=sum(roc_y == 1);      %����������
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
