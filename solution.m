function [x,eps] = solution(A,b,type)  %��Ϊ��������type=1����������type=2%
x_right = (A\b)';
if type == 2
    A = rot90(A,2);
    b = rot90(b,2);
end

[~ ,n] = size(A);
y = b;
if A(n,n) == 0 && b(n) == 0                         %ϵ�������������������ͬ�Ҿ�������
    fprintf("�˷�����ĽⲻΨһ")    
elseif A(n,n) == 0 && b(n) ~= 0                   %ϵ����������������Ȳ���ͬ
    fprintf("�˷������޽�")
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

%�������%

eps = vector_norm(x-x_right,"L1");