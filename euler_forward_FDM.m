function [errall,ord]=euler_forward_FDM(interval,T)
%欧拉向前有限差分方法
%连续的一维对流扩散方程
%ut-a(x)*uxx+b(x)*ux+c(x)*u=f(x,t),x in (xl,xr),t in (0,T].
%u(xl,t)=ue(xl,t),u(xr,t)=ue(xr,t), t in [0,T]
%u(x,0)=v(x),x in (xl,xr)
 
%向前euler差分
%(U^{n}_j-U^{n-1}_j）/k=aj*(U^{n-1}_{j+i}-2*U^{n-1}_j+U^{n-1}_{j-i})/h^2
%           -bj*(U^{n-1}_{j+i}-U^{n-1}_{j-1})/2h-cj*U^{n-1}_j+f^{n-i}_j
%合并
%U^n_j=(aj*lamda-bj*k/2h)U^{n-1}_{j+1}+(-cj*k-2aj*lamda+1)U^{n-1}_{j}
%      +(aj*lamda+bj*k/2h)U^{n-1}_{j-1}+kf^{n-1}_j,j=2,...,M,n=1,...,N
%U^n_1=ue(xl,tn),U^n_{M+1}=ue(xr,tn),n=1,2,...,N
%U^0_j=v(xj),j=1,2,...,M+1
 
%interval=[0,1];T=1;
M=[10,20,40];N=10.*M.^2;
for i=1:length(M)
    [U,err,h]=euler_forward_FDM_solver1(interval,M(i),T,N(i));
    errall(i,:)=err;hall(i)=h;
    if i>1
         ord(i-1,1)=log(errall(i,1)/errall(i-1,1))/log(hall(i)/hall(i-1));
    end
end
end
 
function [U,err,h]=euler_forward_FDM_solver1(interval,M,T,N)
 
%网格剖分
%interval是区间
%P是空间剖分节点坐标
%h是空间步长
P=zeros(M+1,1);
xl=interval(1);xr=interval(2);
h=(xr-xl)/M;
for i=1:M+1
    P(i) = xl+(i-1)*h;
end
%k是时间步长，T是终止时刻
k=T/N;
lamda=k/(h^2);
 
%step2建立参数函数
ue=@(x,t)x.*(x.^2).*cos(x).*cos(2.*pi.*t);
a=@(x)exp(x);b=@(x)cos(x);c=@(x)x.^4;
syms x t
ues=ue(x,t);as=a(x);bs=b(x);cs=c(x);
ut=diff(ues,t);
ux=diff(ues,x);uxx=diff(ux,x);
fs=ut-as*uxx+bs*ux+cs*ues;
f=matlabFunction(fs,'vars',[x,t]);
 
vec_a=a(P(:));vec_b=b(P(:));vec_c=c(P(:));
V=ue(P(:),0);
 
%step3 组装系数矩阵和右端项
%U表示当前时刻tn数值解
%Uold表示上一时刻tp数值解
Uold=V;U=zeros(M+1,1);Uall=zeros(M+1,N);
for n=1:N
    tn=n*k;tp=(n-1)*k;
    for j=1:M+1
        if j==1||j==M+1
            %边界条件
            U(j,1)=ue(P(j),tn);
        else 
            %内部节点
            U(j,1)=(vec_a(j)*lamda-vec_b(j)*k/(2*h))*Uold(j+1,1)...
                +(-vec_c(j)*k-2*vec_a(j)*lamda+1)*Uold(j,1)...
                +(vec_a(j)*lamda+vec_b(j)*k/(2*h))*Uold(j-1,1)...
                +k*f(P(j),tp);
        end
    end
    Uold=U;
    Uall(:,n)=U;
end
 
%step4计算误差
Uexact_final=ue(P(:),T);
error_vector=U-Uexact_final;
errL1=norm(error_vector,1);
errL2=norm(error_vector,2);
errLinf=norm(error_vector,inf);
err=[errLinf,errL2,errL1];
 
%step5画图
 
figure
subplot(1,2,1)
plot(1:M+1,U,'--ro',1:M+1,Uexact_final,':b*')
legend(['numerical solution'],['exact  solution']);
 
subplot(1,2,2)
hold on
for n=1:N
    plot(1:M+1,Uall(:,n),'-r');
end
end
