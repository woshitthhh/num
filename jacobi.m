function [x,m,err]=jacobi(A,b,x0,eps)%%xΪ��������mΪ��������
[n,~]=size(A);
A_diag=diag(A);
A=A-diag(diag(A));
m=0;m_max=50;%%����������������������50������Ϊ�˷�����ĵ����ⷨ������
x=x0;
while m<m_max
    for i=1:n
        x(i)=(b(i)-A(i,:)*x0)/A_diag(i);
    end
    x
    err=norm(x-x0,inf);
    if(err<eps)%%�ж������ڵ����ε��������������Ƿ�С��eps
        break
    end
    x0=x;
    m=m+1;
end

