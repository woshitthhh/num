syms x y t
p_exact = @(x,y,t)x.*(x-1).*cos(2*pi.*t);
p_t = diff(p_exact,t);
p_x = diff(p_exact,x);
p_xx = diff(p_x,x);
p_y = diff(p_exact,y);
p_yy = diff(p_y,y);
f = p_t + 0.01*(p_xx+p_yy)