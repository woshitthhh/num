function [x,d]=gradient_decent(f,eps1,x0)

syms x1 x2
f1=diff(f,x1);
f2=diff(f,x2);
k=1;

x=x0;

grad_x1=subs(f1,x1,x(1));
grad_x=subs(grad_x1,x2,x(2));
grad_y1=subs(f2,x1,x(1));
grad_y=subs(grad_y1,x2,x(2));

d=-[grad_x,grad_y]';
d_norm=double(norm(d));

while d_norm>=eps1
    syms al
    x3 = x + al * d;
    al = solve(diff(subs(subs(f,x1,x3(1)),x2,x3(2)),al)==0,al);
    x = x + al * d;

    grad_x1=subs(f1,x1,x(1));
    grad_x=subs(grad_x1,x2,x(2));
    grad_y1=subs(f2,x1,x(1));
    grad_y=subs(grad_y1,x2,x(2));
    d=-[grad_x,grad_y]'
    k=k+1;
end