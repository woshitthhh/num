function [L,R,X,eps1]=crout(A,B)
%%%����n*n�ķ���A���ͳ�����������B
Y=A\B;
[n,~]=size(A);
L=zeros(n);
R=eye(n);
for k=1:n
    for i=k:n
        L(i,k)=A(i,k)-L(i,1:k-1)*R(1:k-1,k);
    end
    for j=k+1:n
        R(k,j)=(A(k,j)-L(k,1:k-1)*R(1:k-1,j))/L(k,k);
    end
end
y=solution(L,B,2);   %%���������ǻش�����
X=solution(R,y,1);   %%���������ǻش�����
eps1=norm(X-Y,2);  %%������
