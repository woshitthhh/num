function   s = picard(x0,y0,f,n)
syms A;
eps = 10^(-6);
A(1) = y0;
syms x;
for  i = 2:n
    A(i) = y0 + int(f(x,A(i-1)),x0,x);
end
h(x) = -log(1-x)+1;
for i = 1:n
    t = -1:0.1:1;
    g(x) = A(i);
    max_eps = max(abs(h(t) - g(t)));
    s = A(i);
    if(max_eps < eps)
        fprintf("迭代序列在第%d次收敛成功",i);
        break;
    end
end

