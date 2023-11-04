function z=two_point_boundary_value(left_node,right_node,M)
%一维两点边值问题
%Au=-a(x)*u''(x)+b(x)*u'(x)+c(x)*u(x)=f(x),   x\in  (x_left,x_right)
%u(left_node)=ul:u(right_node)=ur;

%有限差分离散格式1
%AhUj=-a(x_j)*(U_(j+1)-2U_j+U(j-1))/h^2+b(x_j)*U(_(j+1)-U(_(j-1))+c(xj)*U_j=f(xj),j=1,2,...,M
%U_0=u1;U_{M+1}=ur;
%离散格式2
%-(aj+1/2*h*bj)*(U_(j-1)-(2aj+h^2*cj)Uj-(aj-1/2*hbj)U_(j+1)=H^2*f(


%step 1;网格剖分,
%M:区间等距剖分
%h:网格尺寸，即剖分后每一段长度
%Point(1,:)
%Point(2,:)
Point=zeros(2,M+1);   %(第一行存节点坐标Point(1,:),Point(2,:)，第二行存判断是否是边界节点的符号, 若为正值（＞0），表示内部节点，若为负值（＜0），为边界节点(
h=(right_node-left_node)/M;

for i=1:M+1
    n=0;
    Point(1,i)=left_node+(i-1)*h;
    if i==1
        Point(2,i)=-1;
    elseif i==M+1
        Point(2,i)=-2;
    else 
        n=n+1;
        Point(2,i)=n;
    end
end
% Step2，建立参数函数，并且向量化
%构造精确解
uexact=@(x)x.^2.*(x-1).*sin(2.*pi.*x);%(a大于零，即离散后的矩阵是正定的）
a_func=@(x)exp(x);
b_func=@(x)cos(pi.*x);
c_func=@(x)x.^2;

a_vec=a_func(Point(1,:));
b_vec=b_func(Point(1,:));
c_vec=c_func(Point(1,:));

%计算右端项
syms x
ue=uexact(x);
a=a_func(x);b=b_func(x);c=c_func(x);
ux=diff(ue,x);
uxx=diff(ux,x);
rhs=-a*uxx+b*ux+c*ue;

f_func=matlabFunction(rhs);
%参数函数向量化
u_vec=uexact(Point(1,:));
a_vec=a_func(Point(1,:));
b_vec=b_func(Point(1,:));
c_vec=c_func(Point(1,:));
f_vec=f_func(Point(1,:));
%Step3:组装系数矩阵及右端向量
%方法1，内部节点是未知量，边界节点是对应已知量

%参数函数在内部节点的取值
a_vec_in=a_vec(2:M);
b_vec_in=b_vec(2:M);
c_vec_in=c_vec(2:M);
f_vec_in=f_vec(2:M);
u_vec_in=u_vec(2:M);
%-(aj+(1/2)*h*bj)*U_(j-1)+(2aj+h^2*cj)Uj-(aj-(1/2)*h*bj)U_(j+1)=h^2*fj
A=sparse(M-1,M-1);
F=zeros(M-1,1);
for k=1:M-1
    %k内部节点编号，表示当前节点是第k个内部节点
    %k+1是整体节点编号中是第k+1个
    %if Point(2,k)<0||Point(2<2+k)<0   
   %else    
   % end
   if k==1
        %内部节点，涉及左边界节点
        A(k,[k,k+1])=[(2*a_vec_in(k)+h^2*c_vec_in(k)),-(a_vec_in(k)-0.5*h*b_vec_in(k))];
            F(k,1)=(a_vec_in(k)+0.5*h*b_vec_in(k))*u_vec(1)+h^2*f_vec_in(k);
   elseif  k==M-1
        %内部节点,涉及右边界节点
        A(k,[k-1,k])=[-(a_vec_in(k)+0.5*h*b_vec_in(k)),(2*a_vec_in(k)+h^2*c_vec_in(k))];
        F(k,1)=(a_vec_in(k)-0.5*h*b_vec_in(k))*u_vec(M+1)+h^2*f_vec_in(k);
   else
        %纯内部节点
        A(k,[k-1,k,k+1])=[-(a_vec_in(k)+0.5*h*b_vec_in(k)),(2*a_vec_in(k)+h^2*c_vec_in(k)),-(a_vec_in(k)-0.5*h*b_vec_in(k))];
        F(k,1)=h^2*f_vec_in(k);
   end
end
%Step4 求解线性方程组
%AU=F  追赶法解矩阵
U=A\F;
%step 5 ：计算误差及收敛阶
        err_vec=u_vec_in'-U;
        L2eff=norm(err_vec,2);
        Inferr=norm(err_vec,inf);
        L1err=norm(err_vec,1);
        
%Step 6 画图
figure
plot(1:M-1,u_vec_in,'-ro',1:M-1,U','-b*')
legend('Exact Solution','Numerical Soultion')
        