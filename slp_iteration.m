function [x,k]=slp_iteration(A,b,x0,type_iterate,eps,alpha)
%%  ����˵��
%A--���Է�����ϵ������
%b--�������Ҷ���������
%x0--��ʼ������
%type_iterate--������ʽ����ѡ'jacobi'���ſɱȵ���������'GS'��Gauss-Seidel����������'SOR'(��γ��ɳڷ�)��
%eps--���ȣ�����ѡ��1e-5��
%alpha--�ɳ����ӣ��ڲ���SOR������ʱ����Ҫ��һ�����ʵ��ɳ�����alpha
%x--������
%k--����������
%%  ���ø�ʽ��
     %first��[x,k]=slp_iteration(A,b,x0,type_iterate,eps)������jacobi��������GS������ʱ��
     %second�� [x,k]=slp_iteration(A,b,x0,type_iterate,eps,alpha)������SOR������ʱ��
%%
[~,m]=size(A);
A_diag=diag(A);  %��ȡ�Խ�Ԫ�أ�
A(1:m+1:end)=0;  %���Խ�Ԫ�ظ�ֵΪ0��
x=x0;
flag=1;  %ѭ����ǣ�
k=0;      %��������
K_max=100; %�����������������������������������Ϊ��������
%%
if min(abs(A_diag))==0
    disp('error������������ϵ�����󲻷��ϵ���Ҫ�󣺶Խ���Ԫ�ز���Ϊ�㣡');
end
switch type_iterate
    case 'jacobi'  
                          %%  Jacobi Methed
        while flag
              for i=1:m
                   x(i)=(b(i)-A(i,:)*x0)/A_diag(i);
              end
              if norm(x-x0,inf)<=eps   %�ж��ط��ﵽ����
                  break
              end
              x0=x;
              k=k+1;
              if k>=K_max
                  disp('�������������ڲ��������������Ǵ˵���������������')
                  break
              end
        end
    case 'GS'
                        %% Gauss-Seidel Methed
        while flag
              for i=1:m
                   x(i)=(b(i)-A(i,:)*x)/A_diag(i);
              end
              if norm(x-x0,inf)<=eps
                  break
              end
              x0=x;
              k=k+1;
              if k>=K_max
                  disp('�������������ڲ��������������Ǵ˵���������������')
                  break
              end
        end
    case 'SOR'
                  %% Successive Overrelaxation Methed
        while flag
              for i=1:m
                   x(i)=(1-alpha)*x0(i)+alpha*(b(i)-A(i,:)*x)/A_diag(i);
              end
              if norm(x-x0,inf)<=eps
                  break
              end
              x0=x;
              k=k+1;
              if k>=K_max
                  disp('�������������ڲ��������������Ǵ˵���������������')
                  break
              end
        end       
end
