function [x,m,err]=Gauss_seidel2(A,b,x0,eps)%%xΪ��������mΪ��������
[n,~]=size(A);
A_diag=diag(A);
A=A-diag(diag(A));
m=0;m_max=50;%%����������������������50������Ϊ�˷�����ĵ����ⷨ������
x=x0;
while m<m_max
    for i=n:-1:1 %%ֻ��ÿ�ε������x����
        x(i)=(b(i)-A(i,:)*x)/A_diag(i);
    end
    x
    err=norm(x-x0,inf);
    if(err<eps)%%�ж������ڵ����ε��������������Ƿ�С��eps
        break
    end
    x0=x;
    m=m+1;
end