function [n] = vector_norm(a,normtype)
a_length = length(a);
n = 0;

%һ����%
if strcmp(normtype,"L1")
    for i = 1:a_length
        n = n + abs(a(i));
    end

    
 %������%   
elseif strcmp(normtype,"L2")
    for i = 1:a_length
        n = n +a(i)^2;
    end
    n = n^(0.5);
    if n == norm(a,2)
        fprintf("������ȷ")
    else
         fprintf("�������")
    end
    
%�����%    
elseif strcmp(normtype,"inf")
    n = max(abs(a));
    if n == norm(a,inf)
        fprintf("������ȷ")
    else
         fprintf("�������")
    end
end





