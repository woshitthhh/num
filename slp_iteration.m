function [x,k]=slp_iteration(A,b,x0,type_iterate,eps,alpha)
%%  参数说明
%A--线性方程组系数矩阵；
%b--方程组右端列向量；
%x0--初始向量；
%type_iterate--迭代方式，可选'jacobi'（雅可比迭代法），'GS'（Gauss-Seidel迭代法）和'SOR'(逐次超松弛法)；
%eps--精度，例如选用1e-5；
%alpha--松弛因子；在采用SOR迭代法时，需要用一个合适的松弛因子alpha
%x--迭代解
%k--迭代次数；
%%  调用格式：
     %first：[x,k]=slp_iteration(A,b,x0,type_iterate,eps)：采用jacobi迭代法或GS迭代法时；
     %second： [x,k]=slp_iteration(A,b,x0,type_iterate,eps,alpha)：采用SOR迭代法时；
%%
[~,m]=size(A);
A_diag=diag(A);  %提取对角元素；
A(1:m+1:end)=0;  %将对角元素赋值为0；
x=x0;
flag=1;  %循环标记；
k=0;      %计数器；
K_max=100; %最大迭代次数，如果超过最大迭代次数，则认为不收敛；
%%
if min(abs(A_diag))==0
    disp('error！！！方程组系数矩阵不符合迭代要求：对角线元素不能为零！');
end
switch type_iterate
    case 'jacobi'  
                          %%  Jacobi Methed
        while flag
              for i=1:m
                   x(i)=(b(i)-A(i,:)*x0)/A_diag(i);
              end
              if norm(x-x0,inf)<=eps   %判断素法达到精度
                  break
              end
              x0=x;
              k=k+1;
              if k>=K_max
                  disp('在最大迭代次数内不收敛！，可能是此迭代方法不收敛。')
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
                  disp('在最大迭代次数内不收敛！，可能是此迭代方法不收敛。')
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
                  disp('在最大迭代次数内不收敛！，可能是此迭代方法不收敛。')
                  break
              end
        end       
end
