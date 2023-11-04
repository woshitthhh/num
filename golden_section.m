function  golden_section(f,eps,A)
B=A;
while A(2) - A(1) > eps
    a1 = A(1);
    b1 = A(2);
    a2 = a1 + 0.382*(b1-a1);
    b2 = a1 + 0.618*(b1-a1);
    if f(a2) > f(b2)
        A(1) = a2;
    else
        A(2) = b2;
    end
end
x = (A(1) + A(2))/2;
y = vpa(f(x),8);
fprintf("此函数在[%d,%d]内的最小值点为%d，最小值为%f\n",B(1),B(2),x,y);
[~,favl] = fminbnd('exp(-x)+x*x',0,1);
error = vpa(y-favl,16);
fprintf("误差为：%f\n",error);