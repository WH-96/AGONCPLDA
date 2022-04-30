clear
clc
load dataSet3.mat %初始关联矩阵非零0.0076 %DAG疾病相似矩阵非零0.6567

%%疾病226
%%基因285
%%已知条件：初始关联矩阵，DAG疾病相似矩阵
warning('off');

        lncSim=miRNASS( LD_adjmat, disSim );  %由DAG疾病相似矩阵得到的lncRNA功能表达相似性非零0.8141
        disSim01  = GSD( LD_adjmat );         %由高斯核得到的疾病相似性矩阵非零 1
        lncSim01  = GSM( LD_adjmat );         %由高斯核得到的基因相似性矩阵非零 1

        disSim02  = combineSim(disSim,disSim01); %整合得到的疾病相似性矩阵非零 1
        lncSim02  = combineSim(lncSim,lncSim01); %整合得到的基因相似性矩阵非零 1
    
        KK=10;  % 相邻个数
        r=0.4;  % 调节权重参数
        ld_adjmat_new=WKNKN( LD_adjmat, lncSim, disSim, KK, r ); % 经过加权KNN处理过的初始关联矩阵非零0.0798

        matPredict=NCPLDA(lncSim02, disSim02, ld_adjmat_new);


[NCP_rank,NCP_rank_known] =Rank_miRNAs( matPredict, LD_adjmat, lncRNA_Name, disease_Name);

Write_file( NCP_rank )
