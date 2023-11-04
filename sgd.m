% 随机梯度下降三维算例matlab代码

% 定义目标函数: f(x1, x2) = x1^2 + 2*x2^2
f = @(X) X(1)^2 + 2 * X(2)^2;

% 初始化参数
x = [3; 5]; % 初始自变量向量
alpha = 0.02; % 学习率
epoch = 300; % 迭代次数
tolerance = 1e-6; % 收敛容忍度
iter = 1; % 迭代计数器

% 计算梯度下降
while iter <= epoch
    % 获取当前位置
    current_x = x;
    % 获取当前损失值 f(current_x)
    current_f = f(current_x);
    
    % 计算当前位置的各方向梯度
    dx1 = 2 * current_x(1);
    dx2 = 4 * current_x(2);
    
    % 随机选取一个方向进行更新（两个方向同等概率）
    if rand() < 0.5
        x(1) = x(1) - alpha * dx1;
    else
        x(2) = x(2) - alpha * dx2;
    end
    
    % 计算新位置处的损失值
    new_f = f(x);
    
    % 更新计数器和误差容忍度
    iter = iter + 1;
    delta_f = abs(new_f - current_f);
    
    % 可视化更新结果
    fprintf('迭代次数: %d, 更新前位置: (%0.2f, %0.2f), 更新后位置: (%0.2f, %0.2f), 损失: %0.5f\n')
end
