function [x,eps1] = Gauss_method(A,b,type)
x_right = (A\b)';
[~,n] = size(A);                                                %方阵A的维度
L = zeros(n,n,n-1);     %对系数矩阵进行行列变换的三维矩阵
%普通方法%
 if type == 1      
         for k = 1:n-1
                    for i = 2:n
                         L(i,k,k) = -A(i,k)/A(k,k);
                    end
                    for m = 1:n
                        L(m,m,k) = 1;
                    end
                    A = L(:,:,k)*A;
                    b = L(:,:,k)*b;
         end
 %列主元素法%
 elseif type == 2
              for k = 1:n-1
                    [~,index] = max(abs(A(k:n,k)));  %寻找每列的最大元素的位置
                    A([k,index+k],:) = A([index+k,k],:) ;      %交换   
                    b([k,index+k],:) = b([index+k,k],:); 
                    for i = 2:n
                         L(i,k,k) = -A(i,k)/A(k,k);
                    end
                    for m = 1:n
                        L(m,m,k) = 1;
                    end
                    A = L(:,:,k)*A;
                    b = L(:,:,k)*b;
              end
%总体主元素法%
 elseif type == 3
     right_order = 1:n;
    for k = 1:n-1
         value=max(max(abs(A(k:n,k:n))));    
        [c,d]=find(value==abs(A(k:n,k:n)));
        c1 = c(1); d1=d(1);
        A([k,c1+k-1],:) = A([c1+k-1,k],:);
        A(:,[k,d1+k-1]) = A(:,[d1+k-1,k]);
        b([1,c1+k-1],:) = b([c1+k-1,1],:);
        right_order([1,d1+k-1]) =  right_order([d1+k-1,1]);
                    for i = 2:n
                         L(i,k,k) = -A(i,k)/A(k,k);
                    end
                    for m = 1:n
                        L(m,m,k) = 1;
                    end
                    A = L(:,:,k)*A;
                    b = L(:,:,k)*b;
    end
 end

y = b;
if A(n,n) == 0 && b(n) == 0                         %系数矩阵与增广矩阵秩相同且均不满秩
    fprintf("此方程组的解不唯一")    
elseif A(n,n) == 0 && b(n) ~= 0                   %系数矩阵与增广矩阵秩不相同
    fprintf("此方程组无解")
elseif A(n,n) ~= 0 && all(b == zeros(n,1))%齐次线性方程组系数矩阵满秩，只有零解
    x = zeros(1,n);
else
    x(n) = b(n)/A(n,n); 
    for i = n-1:-1:1
        for j = n:-1:i+1
            y(i) = y(i)  - A(i,j) * x(j);
        end
        x(i) = y(i)/A(i,i);
    end
end

if type == 3
for i=1:n
        for j=1:n-i
            if(right_order(j)>right_order(j+1))
                [right_order(j),right_order(j+1)]=swap(right_order(j),right_order(j+1));
                [x(j),x(j+1)]=swap(x(j),x(j+1));
            end
        end
end
end
eps1 = vector_norm(x-x_right,"L1");
end



