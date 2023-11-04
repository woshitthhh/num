function [P,Q,X,eps1]=Tridiagonal(A,B)
%%%使用追赶法解三对角矩阵线性方程组：
Y=A\B;
[n,~]=size(A);
P=eye(n);
Q=zeros(n);
Q(1,1)=A(1,1);
for i=2:n
    Q(i-1,i)=A(i-1,i);
    P(i,i-1)=A(i,i-1)/Q(i-1,i-1);
    Q(i,i)=A(i,i)-P(i,i-1)*Q(i-1,i);
end
y= solution(P,B,2);  %%调用下三角回代法函数；
X= solution(Q,y,1);  %%调用上三角回代法函数。
eps1=norm(X-Y,2);  %%计算误差。
