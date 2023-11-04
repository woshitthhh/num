function [weight,bias,k] = PLR_2(input,target,weight0,bias0) %行向量的形式输入input，target.weight0，bias0分别为初始的权矩阵与偏差
[m,~] = size(input);
[~,n] = size(target);
output = zeros(m,n);
k = 0; %计数循环迭代次数
x = weight0;



for i = 1:m
    output(i,:) = hardlim(x*input(i,:)'+bias0)'; %weight0和bias0未迭代时的output
end

while ~isequal(output,target) %每迭代一步 判断所得output与目标分类target是否一致 一致则跳出while 否则继续迭代
    x0=input(:,1);
    y0=input(:,2);
    plot(x0,y0,".");
    hold on
    
    x_1=-5:0.1:5;
    if x(2)~=0
        y_1=(-x(1)*x_1-bias0)/x(2);
        plot(x_1,y_1);
    else
        y_1=-5:0.1:5;
        x_1=-bias0/x(1)*ones(1,101);
        plot(x_1,y_1);
    end
    
    
    if mod(k+1,m) ~= 0
        j = mod(k+1,m);
    else
        j = m;
    end

    r = target(j,:)'-hardlim(x*input(j,:)'+bias0); 
    x = x + r*input(j,:);
    bias0 = bias0 + r;
    k = k+1;

    for i = 1:m
        output(i,:) = hardlim(x*input(i,:)'+bias0); 
    end
end

weight = x;
bias = bias0;

if weight(2)~=0
    y_1=(-weight(1)*x_1-bias)/weight(2);
    plot(x_1,y_1);
else
    y_1=-5:0.1:5;
    x_1=-bias/weight(1)*ones(1,101);
    plot(x_1,y_1);
end
    







