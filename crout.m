function [L,R,X,eps1]=crout(A,B)
%%%输入n*n的方阵A，和常数项列向量B
Y=A\B;
[n,~]=size(A);
L=zeros(n);
R=eye(n);
for k=1:n
    for i=k:n
        L(i,k)=A(i,k)-L(i,1:k-1)*R(1:k-1,k);
    end
    for j=k+1:n
        R(k,j)=(A(k,j)-L(k,1:k-1)*R(1:k-1,j))/L(k,k);
    end
end
y=solution(L,B,2);   %%调用下三角回代法；
X=solution(R,y,1);   %%调用上三角回代法；
eps1=norm(X-Y,2);  %%计算误差。
