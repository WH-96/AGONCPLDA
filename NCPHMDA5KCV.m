
% disease 178
% lncRNA  115
%x(1): DAG�������ƾ���
%x(2):�������ƾ���
%x(3):�������ƾ���
%x(4):�������ƾ���
%x(5):�������ƾ���
%x(6):�������ƾ���
%x(7):WKNN�Ĳ�������r
%x(8):WKNN�ľ���k
%x(9):˫��wknn�ļ���Ȩ��
%x(10):˫��wknn�Ļ���Ȩ��
function [lb,ub,dim,fobj] = NCPHMDA5KCV(F)
switch F
    case 'F1'
        fobj = @F1;
        lb=[0,0,0,0,0,0,0,0,0,0];
        ub=[1,1,1,1,1,1,1,50,1,1];
        dim=10;  
end

end
function o = F1(x) % Tension/compression spring design



%% Load data set
% [~,disease]=xlsread(['Dataset5\disease_178.xlsx']);
% [~,lncRNA]=xlsread(['Dataset5\lncRNA_115.xlsx']);
load('F:\LDA�����¼\�ж���\LDA\NparameterLDA\Dataset5\disease.mat');
load('F:\LDA�����¼\�ж���\LDA\NparameterLDA\Dataset5\lncRNA.mat');

lncSim = load ('Dataset5\lncRNAsimilarity.txt');                            %lncrna�ı��������
interaction = importdata ('Dataset5\known_lncRNA_disease_interaction.txt'); %��֪�ļ���-rna��������
disSim = load('Dataset5\diseasesimilarity.txt');                            %����������������

dis_GS_Sim  = GSD( interaction );         % ������˹��������  �ɸ�˹�˵õ��ļ��������Ծ������ 1
lnc_GS_Sim  = GSM( interaction );         % �����˹��������   �ɸ�˹�˵õ��Ļ��������Ծ������ 1
%% Construct two weight matrices
WD=zeros(178);
index = find(0 ~= disSim);
WD(index) = 0.5;
WL=zeros(115);
index1 = find(0 ~= lncSim);
WL(index1) = 0.5;
[nl,nd] = size(interaction);
%% Retain the original correlation matrix and transpose
interaction = interaction';

%% Calculate the cosine similarity matrix and integrate
[id ,il] = cosSim( interaction );                    % Return the processed cosine similarity matrix


%%
[sd,sl] = integratedsimilarity2(lncSim,disSim,id,il,dis_GS_Sim,lnc_GS_Sim,x);  % Integrated similarity for diseases and miRNAs   
  KK=x(8);  % ���ڸ���
  r=x(7);  % ����Ȩ�ز���
  interaction=WKNKN( interaction',sl,sd,x(8),x(7),x(9),x(10)); % ������ȨKNN������ĳ�ʼ�����������0.0798
  interaction=interaction';
  save interaction interaction;
%% Calculated scoring matrix
[NCP]=NCPHLDA(interaction,sd,sl);


% %% result   ���������
% allresult(disease,lncRNA,interaction,NCP);

%% Five fold cross validation
index_1 = find(1 == interaction);
auc=zeros(1,1);
pp = length(index_1);
   for i = 1 : 1
%     i
    indices = crossvalind('Kfold', pp, 5); %Randomly divide the data sample into 5 parts
    for j = 1:5  %Cycle 5 times, take the i-th part as the test sample and the other two parts as the training samples
       
        index_2 = find(j == indices);
        load interaction;
        interaction(index_1(index_2)) = 0;
        [id ,il] = cosSim( interaction);     %����������������������� �������Ծ����ɳ�ʼ�����������õ� 
        dis_GS_Sim  = GSD(interaction');         % ������˹��������  
        lnc_GS_Sim  = GSM(interaction');         % �����˹��������    
        [sd,sl] = integratedsimilarity2(lncSim,disSim,id,il,dis_GS_Sim,lnc_GS_Sim,x); 
        interaction=WKNKN( interaction',sl,sd,x(8),x(7),x(9),x(10)); % ������ȨKNN������ĳ�ʼ�����������0.0798
        interaction=interaction';
        
        [result]=NCPHLDA(interaction,sd,sl);
        NCP(index_1(index_2)) = result(index_1(index_2)); 

    end
    pre_label_score = NCP(:);
    save pre_label_score_NCPHLDA_5kcv pre_label_score;
    load interaction;
    label_y = interaction(:);
    auc(i) = roc_1(pre_label_score,label_y,'red');
    auc=1/auc;
   end
  
cost=auc;


g(1)=x(1)+x(2)+x(3)-1;
g(2)=x(4)+x(5)+x(6)-1;
g(3)=x(9)+x(10)-1;
% Apply all inequality constraints as a penalty function
lam=10^15;
Z=0;   
for k=1:length(g)
    Z=Z+ lam*g(k)^2*getH(g(k));
end
o=cost+ Z;
end

function H=getH(g)
if g==0 
    H=0; 
else
    H=1; 
end
end
%    auc_avg = mean(auc);
%    std(auc);
 