function [D,L] = cholesky_plus(A,b)
[~,n] = size(A);
d = eig(A);
D = zeros(n,n);
L = eye(n);

D(1,1) = d(1);

