function [x,y] = newton_method(a,f,X) %f为初始函数，X为自变量，a为初始迭代点，eps1为终止误差%
x = a';
while norm(subs(gradient(f,X),X,x'),2) > 0.1
    x = x - subs(hessian(f,X),X,x')\subs(gradient(f,X),X,x');
end
y=subs(f,X,x');