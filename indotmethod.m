clear all
syms d l rho M t phi G ;  
%t=30;
%phi=pi/3;
%G=80;
%lmin=2;
%lmax=1;
%dmin=0.05;
%dmax=0.01;
syms f(d,l) h1(d,l) h2(d,l);
f(d,l)=rho*pi*d*d*l/4;
h1(d,l)=pi*d^3*t-16*M;
h2(d,l)=pi*G*d^4*phi-32*M*l;
%h3(d,l)=d-dmin;
%h4(d,l)=dmax-d;
%h5(d,l)=l-lmin;
%h6(d,l)=lmax-l;
syms lambda p(d,l,lambda);
p(d,l,lambda)=f(d,l)+lambda*h1(d,l)^2+lambda*h2(d,l)^2%+lambda*h3(d,l)^2+lambda*h4(d,l)^2+lambda*h5(d,l)^2+lambda*h6(d,l)^2
syms pd(d,l) pl(d,l)
pd(d,l)=diff(p(d,l,lambda),d)
pl(d,l)=diff(p(d,l,lambda),l)
[d1,l1,lambda1]=solve(pd,pl)
syms he(d,l,lambda)
he(d,l,lambda)=[diff(pd(d,l),d),diff(pd(d,l),l);diff(pl(d,l),d),diff(pl(d,l),l)]
he(d1,l1,lambda1)