function [errall,ord]=crank_nicolson(interval,T)
%crank-nicolson格式
M=[10,20,40];N=10.*M.^2;
for i=1:length(M)
    [U,err,h]=crank_nicolson_solver1(interval,M(i),T,N(i));
    errall(i,:)=err;hall(i)=h;
    if i>1
         ord(i,1)=log(errall(i,1)/errall(i-1,1))/log(hall(i)/hall(i-1));
    end
end
end
 
function [U,err,h]=crank_nicolson_solver1(interval,M,T,N)
 
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
ue=@(x,t)x.*(x-1).*sin(x).*cos(2.*pi.*t);
a=@(x)exp(x);b=@(x)cos(x);c=@(x)x.^4;
syms x t
ues=ue(x,t);as=a(x);bs=b(x);cs=c(x);
ut=diff(ues,t);
ux=diff(ues,x);uxx=diff(ux,x);
fs=ut-as*uxx+bs*ux+cs*ues;
f=matlabFunction(fs,'vars',[x,t]);
V=ue(P(:),0);
 
%step3 组装系数矩阵和右端项
%U表示当前时刻tn数值解
%Uold表示上一时刻tp数值解
Uold=V;
U=zeros(M+1,1);
Uall=zeros(M+1,N);
A=sparse(M+1,M+1);
B=sparse(M+1,M+1);
F=zeros(M+1,1);
G=zeros(M+1,1);
for n=1:N
    tn=n*k;tp=(n-1)*k;
    for i=1:M+1
        if i==1
            A(1,1)=1;
            G(1,1)=ue(P(1),tn);
        elseif i==M+1
            A(M+1,M+1)=1;
            G(M+1,1)=ue(P(M+1),tn);
        else
            A(i,[i-1,i,i+1])=[-lamda/2 1+lamda -lamda/2];
            B(i,[i-1,i,i+1])=[lamda/2 1-lamda lamda/2];
            F(i,1)=(k*0.5)*(f(P(i),tn)+f(P(i),tp));
        end
    end
    U=A\(B*Uold+F+G);
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
