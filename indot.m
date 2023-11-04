clear all
syms d1 d2 l rho g Q W k E;  
%t=30;
%phi=pi/3;
%G=80;
%lmin=2;
%lmax=1;
%dmin=0.05;
%dmax=0.01;
syms f(d1,d2) h1(d1,d2) h2(d1,d2);
f(d1,d2)=pi*l*rho*(2*d1^2+d2^2)/4;
h1(d1,d2)=pi*E*g/(10.67*Q*l^3*W^2*k^2)-(1/d1^4+2.38/d2^4);
h2(d1,d2)=-pi*E*g/(10.67*Q*l^3*W^2*k^2)+(1/d1^4+2.38/d2^4);
%h3(d,l)=d-dmin;
%h4(d,l)=dmax-d;
%h5(d,l)=l-lmin;
%h6(d,l)=lmax-l;
syms lambda p(d1,d2,lambda);
p(d1,d2,lambda)=f(d1,d2)+lambda*h1(d1,d2)^2+lambda*h2(d1,d2)^2%+lambda*h3(d,l)^2+lambda*h4(d,l)^2+lambda*h5(d,l)^2+lambda*h6(d,l)^2
syms pd(d1,d2) pl(d1,d2)
pd(d1,d2)=diff(p(d1,d2,lambda),d1)
pl(d1,d2)=diff(p(d1,d2,lambda),d2)
[d,m,lambda1]=solve(pd,pl)