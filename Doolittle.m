function [L,R,X,eps1]=Doolittle(A,B)
%%%����n*n�ķ���A��������������B.
Y=A\B;
[n,~]=size(A);
L=eye(n);
R=zeros(n);
for k=1:n
    for j=k:n
        R(k,j)=A(k,j)-L(k,1:k-1)*R(1:k-1,j);
    end
    for i=k+1:n
        L(i,k)=(A(i,k)-L(i,1:k-1)*R(1:k-1,k))/R(k,k);
    end
end
y=solution(L,B,2);  %%���������ǻش�����
X=solution(R,y,1);   %%���������ǻش�����
eps1=norm(X-Y,2);  %%������