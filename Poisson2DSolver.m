function Poisson2DSolver(x_interval,y_interval,Mx,My)

%step1:网格剖分
%x_interval:x所在区间,x_interval=[xl,xr]
%y_interval:y所在区间,y_interval=[yl,yr]
%Mx:x区间等距剖分Mx份
%My:y区间等距剖分My份
%P:存储2D节点坐标，P（1，：）X轴坐标。P（2，：）Y轴坐标，P（3，：）点的类型
hx=(x_interval(2)-x_interval(1))/Mx;
hy=(y_interval(2)-y_interval(1))/My;
M=Mx;
%正方形剖分
h=hx;
%矩形剖分
h=sqrt(hx^2+hy^2);

P=zeros(3,(Mx+1)*(My+1));
xl=x_interval(1);xr=x_interval(2);
yl=x_interval(1);yr=y_interval(2);

%固定节点编号：让二维节点以固定的排列方式排点
nn=0;
for i=1:Mx+1
    for j=1:My+1
        k=j+(i-1)*(My+1);
        P(1:2,k)=[xl+(i-1)*hx,yl+(j-1)*hy];
        if j==1&&i~=1&&i~=Mx+1
            %e1:区域的下边界但不包括边界端点
            P(3,k)=-1;
        elseif i==Mx+1&&j~=1&&j~=My+1
            %e2:区域的右边界但不包含边界端点
            P(3,k)=-2;
        elseif j==My+1&&i~=1&&i~=Mx+1
            %e3:区域的上边界但不包含边界端点
            P(3,k)=-3;
        elseif i==1&&j~=1&&j~=My+1
            %e4;区域的左边界但不包含边界端点
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

%step2:建立参数函数及向量化
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

%step3:合成系数矩阵和右端向量(方式一：全部节点都是未知量)
A=sparse((Mx+1)*(My+1),(Mx+1)*(My+1));
F=zeros((Mx+1)*(My+1),1);
for k=1:(Mx+1)*(My+1)
    if P(3,k)>0 %内部节点
        A(k,[k-(Mx+1),k-1,k,k+1,k+Mx+1])=[-1,-1,4,-1,-1];
        F(k,1)=h^2*rhsf(k,1);
    else
        A(k,k)=1;
        F(k,1)=Ue(k,1);
    end
end
%step4:求解方程组
U=A\F;

%step5:计算误差及收敛阶
err_vec=U-Ue;
errL2=norm(err_vec,2);
errL1=norm(err_vec,1);
errinf=norm(err_vec,inf);


%step6:画图
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
    




