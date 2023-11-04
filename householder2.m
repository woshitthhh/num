function [R,b,X]=householder2(A,B)
%%%输入方程组系数矩阵A和右端列向量B;
%%%输出上三角矩阵R，右端列向量b(前n个元素组成的列向量)，方程的解X。
x=A\B;
[m,n]=size(A);
for k=1:n-1            %对右下角矩阵及进行householder变换
    temp_A=A(k:m,k:n);  %取出A的右下角的块矩阵；
    temp_b=B(k:m,1);    %取出B下方对应的;
    a=temp_A(:,1);      %取出矩阵temp_A的第一列；
    g=zeros(m-k+1,1);
    g(1)=1;             %构造单位向量（1，0，...，0）'
    a_2=norm(a,2);  
    u=a-sign(a(1))*a_2*g;     
    w=u/norm(u,2);  
    Ak=temp_A-2*w*(w'*temp_A);  %（I-2*w*w')*temp_A，
    A(k:m,k:n)=Ak;               
    bk=temp_b-2*w*(w'*temp_b);       %（I-2*w*w')*temp_b;
    B(k:m,1)=bk;
end
R=A(1:n,1:n);
b=B(1:n,1);
X=solution(R,b,1);  
ess=1e-6;
if(norm(X-x)<ess)
    fprintf("误差允许");
end
