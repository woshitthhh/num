function [X1, X2] = outpoint(f, h)
M = 1;         %��ֵ
eps =  0.001;    %���
syms x1 x2;             %�������
[m,~]=size(h);
% 
% while true
    p_sum=0;
    for i=1:m
        p_sum=p_sum+h(i,1);
    end
    P = f + M*p_sum        %���캯��
    x1=-10:0.1:10;
    x2=x1;
    y=P(x1,x2);
    plot3(x1,x2,y)
    P_x1 = diff(P,x1);      %��x1ƫ��
    P_x2 = diff(P,x2);    %��x2ƫ��
%     X = solve(P_x1==0,P_x2==0,x1,x2); %�ⷽ����
%     X1 = double(X.x1);
%     X2 = double(X.x2);
%     if M * abs(subs(p_sum,x1,X1)) < eps  %�ж���ֹ����
%         break;
%     end
%     M = M * 2;              %�ӱ�
end

