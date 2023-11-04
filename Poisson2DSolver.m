function Poisson2DSolver(x_interval,y_interval,Mx,My)

%step1:�����ʷ�
%x_interval:x��������,x_interval=[xl,xr]
%y_interval:y��������,y_interval=[yl,yr]
%Mx:x����Ⱦ��ʷ�Mx��
%My:y����Ⱦ��ʷ�My��
%P:�洢2D�ڵ����꣬P��1������X�����ꡣP��2������Y�����꣬P��3�������������
hx=(x_interval(2)-x_interval(1))/Mx;
hy=(y_interval(2)-y_interval(1))/My;
M=Mx;
%�������ʷ�
h=hx;
%�����ʷ�
h=sqrt(hx^2+hy^2);

P=zeros(3,(Mx+1)*(My+1));
xl=x_interval(1);xr=x_interval(2);
yl=x_interval(1);yr=y_interval(2);

%�̶��ڵ��ţ��ö�ά�ڵ��Թ̶������з�ʽ�ŵ�
nn=0;
for i=1:Mx+1
    for j=1:My+1
        k=j+(i-1)*(My+1);
        P(1:2,k)=[xl+(i-1)*hx,yl+(j-1)*hy];
        if j==1&&i~=1&&i~=Mx+1
            %e1:������±߽絫�������߽�˵�
            P(3,k)=-1;
        elseif i==Mx+1&&j~=1&&j~=My+1
            %e2:������ұ߽絫�������߽�˵�
            P(3,k)=-2;
        elseif j==My+1&&i~=1&&i~=Mx+1
            %e3:������ϱ߽絫�������߽�˵�
            P(3,k)=-3;
        elseif i==1&&j~=1&&j~=My+1
            %e4;�������߽絫�������߽�˵�
            P(3,k)=-4;
        elseif i==1&&j==1
            P(3,k)=-11;
        elseif i==1&&j==My+1
            P(3,k)=-12;    
        elseif i==Mx+1&&j==My+1
            P(3,k)=-13;  
        elseif i==Mx+1&&j==1
            P(3,k)=-14;    
        else
            nn=nn+1;
            P(3,k)=nn;
        end
    end
end

figure
hold on
for k=1:size(P,2)
    if P(3,k)>0
        plot(P(1,k),P(2,k),'*r')
    elseif P(3,k)<0
        plot(P(1,k),P(2,k),'*b')
    end
    text(P(1,k)+0.02,P(2,k)+0.02,num2str(k))
end

%step2:��������������������
ue=@(x,y) x.*(x-1).*y.*(y-1).*exp(x+y);
syms x y
uex=diff(ue,x);
uexx=diff(uex,x);
uey=diff(ue,y);
ueyy=diff(uey,y);
f=-uexx-ueyy;
rhs_f=matlabFunction(f);
Ue=ue(P(1,:),P(2,:))';
rhsf=rhs_f(P(1,:),P(2,:))';

%step3:�ϳ�ϵ��������Ҷ�����(��ʽһ��ȫ���ڵ㶼��δ֪��)
A=sparse((Mx+1)*(My+1),(Mx+1)*(My+1));
F=zeros((Mx+1)*(My+1),1);
for k=1:(Mx+1)*(My+1)
    if P(3,k)>0 %�ڲ��ڵ�
        A(k,[k-(Mx+1),k-1,k,k+1,k+Mx+1])=[-1,-1,4,-1,-1];
        F(k,1)=h^2*rhsf(k,1);
    else
        A(k,k)=1;
        F(k,1)=Ue(k,1);
    end
end
%step4:��ⷽ����
U=A\F;

%step5:������������
err_vec=U-Ue;
errL2=norm(err_vec,2);
errL1=norm(err_vec,1);
errinf=norm(err_vec,inf);


%step6:��ͼ
for i=1:M+1
    X(:,i)=P(1,1+(i-1)*(M+1):i*(M+1));
    Y(:,i)=P(2,1+(i-1)*(M+1):i*(M+1));
    Z(:,i)=U(1+(i-1)*(M+1):i*(M+1),1);
    Ze(:,i)=Ue(1+(i-1)*(M+1):i*(M+1),1);
end

figure
surf(X,Y,Z)
title('Numerical')
figure
surf(X,Y,Ze)
title('Exact')
    




