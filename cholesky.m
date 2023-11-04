function [L,X,eps1]=cholesky(A,B)
%%%����n*n����������A�ͳ�����������B��
Y=A\B;
[n,~]=size(A);
L=zeros(n);
for j=1:n
    for i=j:n
        if i==j
        L(j,j)=(A(j,j)-L(j,1:j-1)*(L(j,1:j-1))')^(1/2);
        elseif i>j
        L(i,j)=(A(i,j)-L(i,1:j-1)*(L(j,1:j-1))')/L(j,j);
        end
    end 
end
y= solution(L,B,2);  %%���������ǻش���������
X= solution(L',y,1);  %%���������ǻش���������
eps1=norm(X-Y,2);  %%������



