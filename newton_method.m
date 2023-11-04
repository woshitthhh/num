function [x,y] = newton_method(a,f,X) %fΪ��ʼ������XΪ�Ա�����aΪ��ʼ�����㣬eps1Ϊ��ֹ���%
x = a';
while norm(subs(gradient(f,X),X,x'),2) > 0.1
    x = x - subs(hessian(f,X),X,x')\subs(gradient(f,X),X,x');
end
y=subs(f,X,x');