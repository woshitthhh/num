function [errall,ord]=euler_backward_FDM(interval,T)
%ŷ��������޲�ַ���
%������һά������ɢ����
%ut-a(x)*uxx+b(x)*ux+c(x)*u=f(x,t),x in (xl,xr),t in (0,T].
%u(xl,t)=ue(xl,t),u(xr,t)=ue(xr,t), t in [0,T]
%u(x,0)=v(x),x in (xl,xr)
 
%���euler���
%(U^{n}_j-U^{n-1}_j��/k=aj*(U^{n}_{j+i}-2*U^{n}_j+U^{n}_{j-1})/h^2
%           +bj*(U^{n}_{j+i}-U^{n}_{j-1})/2h+cj*U^{n}_j+k*f^{n}_j
%(-aj*lamda+bj*k/2h)U^{n}_{j-1}+(2aj*lamda+cj*k+1)U^{n}_{j}+(-aj*lamda-bj*k/2h)U^{n}_{j+1}=U^{n-1}_j+kf^{n}_j
%U^n_1=ue(xl,tn),U^n_{M+1}=ue(xr,tn),n=1,2,...,N
%U^0_j=v(xj),j=1,2,...,M+1
%interval=[0,1];T=1;
M=[10,20,40];N=10.*M.^2;
for i=1:length(M)
    [U,err,h]=euler_backward_FDM_solver1(interval,M(i),T,N(i));
    errall(i,:)=err;hall(i)=h;
    if i>1
         ord(i,1)=log(errall(i,1)/errall(i-1,1))/log(hall(i)/hall(i-1));
    end
end
end
 
function [U,err,h]=euler_backward_FDM_solver1(interval,M,T,N)
 
%�����ʷ�
%interval������
%P�ǿռ��ʷֽڵ�����
%h�ǿռ䲽��
P=zeros(M+1,1);
xl=interval(1);xr=interval(2);
h=(xr-xl)/M;
for i=1:M+1
    P(i) = xl+(i-1)*h;
end
%k��ʱ�䲽����T����ֹʱ��
k=T/N;
lamda=k/(h^2);
 
%step2������������
ue=@(x,t)x.*(x-1).*sin(x).*cos(2.*pi.*t);
a=@(x)exp(x);b=@(x)cos(x);c=@(x)x.^4;
syms x t
ues=ue(x,t);as=a(x);bs=b(x);cs=c(x);
ut=diff(ues,t);
ux=diff(ues,x);uxx=diff(ux,x);
fs=ut-as*uxx+bs*ux+cs*ues;
f=matlabFunction(fs,'vars',[x,t]);
 
vec_a=a(P(:));vec_b=b(P(:));vec_c=c(P(:));
V=ue(P(:),0);
 
%step3 ��װϵ��������Ҷ���
%U��ʾ��ǰʱ��tn��ֵ��
%Uold��ʾ��һʱ��tp��ֵ��
Uold=V;U=zeros(M+1,1);Uall=zeros(M+1,N);
A=sparse(M+1,M+1);
F=zeros(M+1,1);
for n=1:N
    tn=n*k;tp=(n-1)*k;
    for j=1:M+1
        if j==1||j==M+1
            %�߽�����
            A(j,j)=1;
            F(j,1)=ue(P(j),tn);
        else
            A(j,[j-1 j j+1])=[-vec_a(j)*lamda-vec_b(j)*k/(2*h) vec_c(j)*k+2*vec_a(j)*lamda+1  -vec_a(j)*lamda+vec_b(j)*k/(2*h)];
            F(j,1)=Uold(j,1)+k*f(P(j),tn);
        end
    end
    U=A\F;
    Uold=U;
    Uall(:,n)=U;
end
 
%step4�������
Uexact_final=ue(P(:),T);
error_vector=U-Uexact_final;
errL1=norm(error_vector,1);
errL2=norm(error_vector,2);
errLinf=norm(error_vector,inf);
err=[errLinf,errL2,errL1];
 
%step5��ͼ
 
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
