function [MD_mat_new] = WKNKN( MD_mat, MM_mat, DD_mat, K, r )

[rows,cols]=size(MD_mat);
y_m=zeros(rows,cols);  
y_d=zeros(rows,cols);  

knn_network_m = KNN( MM_mat, K );  %由基因相似矩阵找出跟基因最亲近的K个基因
for i = 1 : rows   
         w=zeros(1,K);
        [sort_m,idx_m]=sort(knn_network_m(i,:),2,'descend'); %sort_m为与基因i最相近的10个基因
        sum_m=sum(sort_m(1,1:K));   
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_m(1,j); 
            y_m(i,:) =  y_m(i,:)+ w(1,j)* MD_mat(idx_m(1,j),:); 
        end                      
            y_m(i,:)=y_m(i,:)/sum_m;              
end

knn_network_d = KNN( DD_mat , K );  %由疾病相似矩阵找出跟疾病最亲近的K个疾病
for i = 1 : cols   
        w=zeros(1,K);
        [sort_d,idx_d]=sort(knn_network_d(i,:),2,'descend');
        sum_d=sum(sort_d(1,1:K));
        for j = 1 : K
            w(1,j)=r^(j-1)*sort_d(1,j);
            y_d(:,i) =  y_d(:,i)+ w(1,j)* MD_mat(:,idx_d(1,j)); 
        end                      
            y_d(:,i)=y_d(:,i)/sum_d;               
end

a1=0;
a2=1;
y_md=(y_m*a1+y_d*a2)/(a1+a2);  %由加权kNN得到的关联矩阵

 for i = 1 : rows
     for j = 1 : cols
         MD_mat_new(i,j)=max(MD_mat(i,j),y_md(i,j));
     end    
 end

end

function [ knn_network ] = KNN( network , k )
    [rows, cols] = size( network );
    network= network-diag(diag(network));  %diag(A) 返回A的主对角线元素的列向量  diag(V)返回包含主对角线上向量 v 的元素的对角矩阵。
    knn_network = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');%对矩阵的每一行按降序排序
    for i = 1 : rows
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
    end
end


