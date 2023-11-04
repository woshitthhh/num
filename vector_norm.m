function [n] = vector_norm(a,normtype)
a_length = length(a);
n = 0;

%一范数%
if strcmp(normtype,"L1")
    for i = 1:a_length
        n = n + abs(a(i));
    end

    
 %二范数%   
elseif strcmp(normtype,"L2")
    for i = 1:a_length
        n = n +a(i)^2;
    end
    n = n^(0.5);
    if n == norm(a,2)
        fprintf("计算正确")
    else
         fprintf("计算错误")
    end
    
%无穷范数%    
elseif strcmp(normtype,"inf")
    n = max(abs(a));
    if n == norm(a,inf)
        fprintf("计算正确")
    else
         fprintf("计算错误")
    end
end





