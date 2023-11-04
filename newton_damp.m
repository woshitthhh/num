function [x,y] = newton_damp(a,f,X) %fΪ��ʼ������XΪ�Ա�����aΪ��ʼ�����㣬eps1Ϊ��ֹ���%
x = a';
syms b;
while norm(subs(gradient(f,X),X,x'),2) > 0.1
    d = subs(hessian(f,X),X,x')\subs(gradient(f,X),X,x');
    c = solve(diff(subs(f,X,x'+b*d),b)==0);
    x = x - c*d;
end
y=subs(f,X,x');