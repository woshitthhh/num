function [X1, X2] = outpoint(f, h)
M = 1;         %罚值
eps =  0.001;    %误差
syms x1 x2;             %定义变量
[m,~]=size(h);
% 
% while true
    p_sum=0;
    for i=1:m
        p_sum=p_sum+h(i,1);
    end
    P = f + M*p_sum        %构造函数
    x1=-10:0.1:10;
    x2=x1;
    y=P(x1,x2);
    plot3(x1,x2,y)
    P_x1 = diff(P,x1);      %求x1偏导
    P_x2 = diff(P,x2);    %求x2偏导
%     X = solve(P_x1==0,P_x2==0,x1,x2); %解方程组
%     X1 = double(X.x1);
%     X2 = double(X.x2);
%     if M * abs(subs(p_sum,x1,X1)) < eps  %判断终止条件
%         break;
%     end
%     M = M * 2;              %加倍
end

