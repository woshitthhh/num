function [s] = matrix_norm(A,normtype)
[m,n] = size(A);
s = 0;

%一范数%
if strcmp(normtype,"L1")
    for i = 1:m
        b(i) = vector_norm(A(:,i),"L1");
    end
    s = max(b);
    if s == norm(A,1)
        fprintf("计算正确")
    else
         fprintf("计算错误")
    end
    
%二范数%     
elseif strcmp(normtype,"L2")    
    s = (max(eig(A'*A)))^(0.5);
    if s == norm(A,2)
        fprintf("计算正确")
    else
         fprintf("计算错误")
    end
    
%无穷范数%     
elseif strcmp(normtype,"inf")
    for i = 1:n
        b(i) = 
    end
    s = max(a);
    if s == norm(A,inf)
        fprintf("计算正确")
    else
         fprintf("计算错误")
    end
end

    