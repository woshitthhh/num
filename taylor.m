syms x;
x = -15*pi:0.5*pi:15*pi;
y1 = sin(x);
plot(x,y1);
hold on;
for k = 1:100 
    for j = 1:61
        x = -15*pi + (j-1)*0.5*pi;
        i = 1:k;
        s = x^(2*i-1)/(2*i-1)*(-1)^(i+1);
        y2(j) = sum(s);
    end
        plot(x,y2);
        hold on;
end
