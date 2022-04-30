function allresult(disease,lncRNA,interaction,NCP)
a=(interaction-1)*(-1);                 % 无关联矩阵
NCP=NCP.*a;                             %去除有关联得分 已知关联全部为零
[m,n]=size(NCP);                        % m:疾病  n:基因
c=find(NCP~=0);
[~,k]=size(c);
p=0;
str=cell(k,3);
for i=1:m
    for j=1:n
        if NCP(i,j)~=0
            p=p+1;
            str(p,1)=disease(i,1);
            str(p,2)=lncRNA(j,1);           
            str{p,3}=NCP(i,j);
        end
    end
end
str=sortrows(str,[1 3],'descend');   %先对第一列进行降序-排序 在对第三列进行降序排序                   
xlswrite('.\allresult2.xlsx',str);
end
