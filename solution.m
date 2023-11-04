function [x,eps] = solution(A,b,type)  %若为上三角阵type=1，下三角阵type=2%
x_right = (A\b)';
if type == 2
    A = rot90(A,2);
    b = rot90(b,2);
end

[~ ,n] = size(A);
y = b;
if A(n,n) == 0 && b(n) == 0                         %系数矩阵与增广矩阵秩相同且均不满秩
    fprintf("此方程组的解不唯一")    
elseif A(n,n) == 0 && b(n) ~= 0                   %系数矩阵与增广矩阵秩不相同
    fprintf("此方程组无解")
else
    x(n) = b(n)/A(n,n); 
    for i = n-1:-1:1
        for j = n:-1:i+1
            y(i) = y(i)  - A(i,j) * x(j);
        end
        x(i) = y(i)/A(i,i);
    end
end
x=x';

if type == 2
    x = rot90(x,2);
end

%检验误差%

eps = vector_norm(x-x_right,"L1");