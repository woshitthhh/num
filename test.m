n1 = 1.5;
L1 = 15.5;
syms alpha bete rho varphi %n1 L1
f1 = sin(alpha);
f2 = n1*sin(bete);
f3 = rho*sin(varphi)-L1*tan(alpha);
f4 = (rho*cos(varphi)-L1)*tan(bete);
f5 = f1-f2 == 0;
f6 = f3-f4 == 0;
eq = [f5,f6];
[x11,x22] = solve(eq,[alpha,bete])