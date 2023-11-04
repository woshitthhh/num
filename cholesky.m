function [L,X,eps1]=cholesky(A,B)
%%%输入n*n的正定矩阵A和常数项列向量B，
Y=A\B;
[n,~]=size(A);
L=zeros(n);
for j=1:n
    for i=j:n
        if i==j
        L(j,j)=(A(j,j)-L(j,1:j-1)*(L(j,1:j-1))')^(1/2);
        elseif i>j
        L(i,j)=(A(i,j)-L(i,1:j-1)*(L(j,1:j-1))')/L(j,j);
        end
    end 
end
y= solution(L,B,2);  %%调用下三角回代法函数；
X= solution(L',y,1);  %%调用上三角回代法函数。
eps1=norm(X-Y,2);  %%计算误差。



