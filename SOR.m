function [x,m,err]=SOR(A,b,x0,eps,a)%%xΪ��������mΪ����������aΪ�ɳ�����
[n,~]=size(A);
A_diag=diag(A);
m=0;m_max=50;%%����������������������50������Ϊ�˷�����ĵ����ⷨ������
x=x0;
while m<m_max
    for i=1:n
        x(i)=x0(i)+(b(i)-A(i,:)*x)*a/A_diag(i);
    end
    x
    err=norm(x-x0,inf);
    if(err<eps)%%�ж������ڵ����ε��������������Ƿ�С��eps
       break
    end
    x0=x;
    m=m+1;
end