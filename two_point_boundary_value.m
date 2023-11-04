function z=two_point_boundary_value(left_node,right_node,M)
%һά�����ֵ����
%Au=-a(x)*u''(x)+b(x)*u'(x)+c(x)*u(x)=f(x),   x\in  (x_left,x_right)
%u(left_node)=ul:u(right_node)=ur;

%���޲����ɢ��ʽ1
%AhUj=-a(x_j)*(U_(j+1)-2U_j+U(j-1))/h^2+b(x_j)*U(_(j+1)-U(_(j-1))+c(xj)*U_j=f(xj),j=1,2,...,M
%U_0=u1;U_{M+1}=ur;
%��ɢ��ʽ2
%-(aj+1/2*h*bj)*(U_(j-1)-(2aj+h^2*cj)Uj-(aj-1/2*hbj)U_(j+1)=H^2*f(


%step 1;�����ʷ�,
%M:����Ⱦ��ʷ�
%h:����ߴ磬���ʷֺ�ÿһ�γ���
%Point(1,:)
%Point(2,:)
Point=zeros(2,M+1);   %(��һ�д�ڵ�����Point(1,:),Point(2,:)���ڶ��д��ж��Ƿ��Ǳ߽�ڵ�ķ���, ��Ϊ��ֵ����0������ʾ�ڲ��ڵ㣬��Ϊ��ֵ����0����Ϊ�߽�ڵ�(
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
% Step2��������������������������
%���쾫ȷ��
uexact=@(x)x.^2.*(x-1).*sin(2.*pi.*x);%(a�����㣬����ɢ��ľ����������ģ�
a_func=@(x)exp(x);
b_func=@(x)cos(pi.*x);
c_func=@(x)x.^2;

a_vec=a_func(Point(1,:));
b_vec=b_func(Point(1,:));
c_vec=c_func(Point(1,:));

%�����Ҷ���
syms x
ue=uexact(x);
a=a_func(x);b=b_func(x);c=c_func(x);
ux=diff(ue,x);
uxx=diff(ux,x);
rhs=-a*uxx+b*ux+c*ue;

f_func=matlabFunction(rhs);
%��������������
u_vec=uexact(Point(1,:));
a_vec=a_func(Point(1,:));
b_vec=b_func(Point(1,:));
c_vec=c_func(Point(1,:));
f_vec=f_func(Point(1,:));
%Step3:��װϵ�������Ҷ�����
%����1���ڲ��ڵ���δ֪�����߽�ڵ��Ƕ�Ӧ��֪��

%�����������ڲ��ڵ��ȡֵ
a_vec_in=a_vec(2:M);
b_vec_in=b_vec(2:M);
c_vec_in=c_vec(2:M);
f_vec_in=f_vec(2:M);
u_vec_in=u_vec(2:M);
%-(aj+(1/2)*h*bj)*U_(j-1)+(2aj+h^2*cj)Uj-(aj-(1/2)*h*bj)U_(j+1)=h^2*fj
A=sparse(M-1,M-1);
F=zeros(M-1,1);
for k=1:M-1
    %k�ڲ��ڵ��ţ���ʾ��ǰ�ڵ��ǵ�k���ڲ��ڵ�
    %k+1������ڵ������ǵ�k+1��
    %if Point(2,k)<0||Point(2<2+k)<0   
   %else    
   % end
   if k==1
        %�ڲ��ڵ㣬�漰��߽�ڵ�
        A(k,[k,k+1])=[(2*a_vec_in(k)+h^2*c_vec_in(k)),-(a_vec_in(k)-0.5*h*b_vec_in(k))];
            F(k,1)=(a_vec_in(k)+0.5*h*b_vec_in(k))*u_vec(1)+h^2*f_vec_in(k);
   elseif  k==M-1
        %�ڲ��ڵ�,�漰�ұ߽�ڵ�
        A(k,[k-1,k])=[-(a_vec_in(k)+0.5*h*b_vec_in(k)),(2*a_vec_in(k)+h^2*c_vec_in(k))];
        F(k,1)=(a_vec_in(k)-0.5*h*b_vec_in(k))*u_vec(M+1)+h^2*f_vec_in(k);
   else
        %���ڲ��ڵ�
        A(k,[k-1,k,k+1])=[-(a_vec_in(k)+0.5*h*b_vec_in(k)),(2*a_vec_in(k)+h^2*c_vec_in(k)),-(a_vec_in(k)-0.5*h*b_vec_in(k))];
        F(k,1)=h^2*f_vec_in(k);
   end
end
%Step4 ������Է�����
%AU=F  ׷�Ϸ������
U=A\F;
%step 5 ��������������
        err_vec=u_vec_in'-U;
        L2eff=norm(err_vec,2);
        Inferr=norm(err_vec,inf);
        L1err=norm(err_vec,1);
        
%Step 6 ��ͼ
figure
plot(1:M-1,u_vec_in,'-ro',1:M-1,U','-b*')
legend('Exact Solution','Numerical Soultion')
        