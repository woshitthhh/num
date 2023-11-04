function [y] = sort_norm(n)
X = rand(1,n)

syms y;
for i = 1:9
    y(i) = i;
end
y(10) = "inf";

y

a = R_norm(X,y(1))
for i = 1:10
    for j = i+1:10
        if R_norm(X,y(j)) < R_norm(X,y(j+1))
            y([j,j+1]) = y([j+1,j]);
        end
    end
end
