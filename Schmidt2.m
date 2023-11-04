  function c=Schmidt2(a)
  [~,m] = size(a);
  syms t;
 %Õý½»»¯
  c(1)=a(1)/(L2_dot(a(1),a(1),exp(-t^2),-inf,inf))^0.5;
  d = 0;
 for i=2:m
     for j=1:i-1 
        d = d + L2_dot(a(i),c(i-1),exp(-t^2),-inf,inf)*c(i-1);
     end
     c(i)= (a(i) - d)/L2_dot(a(i),a(i),exp(-t^2),-inf,inf)^0.5;
 end
