x=-1:0.25:1.25;
y=[0.2209 0.3295 0.8826 1.4392 2.0003 2.5645 3.1334 3.7061 4.2836 3.1334];
sum=0;
for i=1:1:10
    yy(i)=1.786*x(i)+1.946;
    sum=sum+(y(i)-yy(i))^2;
end
sum
