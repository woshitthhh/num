function z = L2_dot(x,y,z,a,b) %zΪȨ����
syms t;
z = int(x*y*z,t,a,b);
