clear
clc
load dataSet3.mat %��ʼ�����������0.0076 %DAG�������ƾ������0.6567

%%����226
%%����285
%%��֪��������ʼ��������DAG�������ƾ���
warning('off');

        lncSim=miRNASS( LD_adjmat, disSim );  %��DAG�������ƾ���õ���lncRNA���ܱ�������Է���0.8141
        disSim01  = GSD( LD_adjmat );         %�ɸ�˹�˵õ��ļ��������Ծ������ 1
        lncSim01  = GSM( LD_adjmat );         %�ɸ�˹�˵õ��Ļ��������Ծ������ 1

        disSim02  = combineSim(disSim,disSim01); %���ϵõ��ļ��������Ծ������ 1
        lncSim02  = combineSim(lncSim,lncSim01); %���ϵõ��Ļ��������Ծ������ 1
    
        KK=10;  % ���ڸ���
        r=0.4;  % ����Ȩ�ز���
        ld_adjmat_new=WKNKN( LD_adjmat, lncSim, disSim, KK, r ); % ������ȨKNN������ĳ�ʼ�����������0.0798

        matPredict=NCPLDA(lncSim02, disSim02, ld_adjmat_new);


[NCP_rank,NCP_rank_known] =Rank_miRNAs( matPredict, LD_adjmat, lncRNA_Name, disease_Name);

Write_file( NCP_rank )
