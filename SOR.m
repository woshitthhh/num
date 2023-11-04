function [x,m,err]=SOR(A,b,x0,eps,a)%%x为解向量，m为迭代次数，a为松弛因子
[n,~]=size(A);
A_diag=diag(A);
m=0;m_max=50;%%计数器，若迭代次数大于50次则认为此方程组的迭代解法不收敛
x=x0;
while m<m_max
    for i=1:n
        x(i)=x0(i)+(b(i)-A(i,:)*x)*a/A_diag(i);
    end
    x
    err=norm(x-x0,inf);
    if(err<eps)%%判断若相邻的两次迭代结果间的误差精度是否小于eps
       break
    end
    x0=x;
    m=m+1;
end