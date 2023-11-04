function [R,b,X]=householder2(A,B)
%%%���뷽����ϵ������A���Ҷ�������B;
%%%��������Ǿ���R���Ҷ�������b(ǰn��Ԫ����ɵ�������)�����̵Ľ�X��
x=A\B;
[m,n]=size(A);
for k=1:n-1            %�����½Ǿ��󼰽���householder�任
    temp_A=A(k:m,k:n);  %ȡ��A�����½ǵĿ����
    temp_b=B(k:m,1);    %ȡ��B�·���Ӧ��;
    a=temp_A(:,1);      %ȡ������temp_A�ĵ�һ�У�
    g=zeros(m-k+1,1);
    g(1)=1;             %���쵥λ������1��0��...��0��'
    a_2=norm(a,2);  
    u=a-sign(a(1))*a_2*g;     
    w=u/norm(u,2);  
    Ak=temp_A-2*w*(w'*temp_A);  %��I-2*w*w')*temp_A��
    A(k:m,k:n)=Ak;               
    bk=temp_b-2*w*(w'*temp_b);       %��I-2*w*w')*temp_b;
    B(k:m,1)=bk;
end
R=A(1:n,1:n);
b=B(1:n,1);
X=solution(R,b,1);  
ess=1e-6;
if(norm(X-x)<ess)
    fprintf("�������");
end
