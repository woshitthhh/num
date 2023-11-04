function z = L2_dot(x,y,z,a,b) %zÎªÈ¨º¯Êý
syms t;
z = int(x*y*z,t,a,b);
