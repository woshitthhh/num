function [s] = matrix_norm(A,normtype)
[m,n] = size(A);
s = 0;

%һ����%
if strcmp(normtype,"L1")
    for i = 1:m
        b(i) = vector_norm(A(:,i),"L1");
    end
    s = max(b);
    if s == norm(A,1)
        fprintf("������ȷ")
    else
         fprintf("�������")
    end
    
%������%     
elseif strcmp(normtype,"L2")    
    s = (max(eig(A'*A)))^(0.5);
    if s == norm(A,2)
        fprintf("������ȷ")
    else
         fprintf("�������")
    end
    
%�����%     
elseif strcmp(normtype,"inf")
    for i = 1:n
        b(i) = 
    end
    s = max(a);
    if s == norm(A,inf)
        fprintf("������ȷ")
    else
         fprintf("�������")
    end
end

    