function [x] = R_norm(X,a)
if a ~= inf
    x = sum(abs(X).^a)^(1/a);
else
    x = max(abs(X));
end
