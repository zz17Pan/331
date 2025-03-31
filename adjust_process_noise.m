function Q = adjust_process_noise(velocities, base_params)
%ADJUST_PROCESS_NOISE 自适应调整过程噪声以适应不同的目标运动状态
%   velocities: 目标速度状态 [径向速度, 方位角速度, 俯仰角速度]
%   base_params: 基础参数
%   返回调整后的过程噪声协方差矩阵Q

% 提取速度分量
vr = velocities(1);   % 径向速度
vaz = velocities(2);  % 方位角速度
vel = velocities(3);  % 俯仰角速度

% 保存历史速度用于加速度评估
persistent prev_velocities;
persistent prev_timestamp;

if isempty(prev_velocities)
    prev_velocities = velocities;
    prev_timestamp = now;
else
    % 估计时间差 (秒)
    dt = max(0.05, (now - prev_timestamp) * 86400);
    prev_timestamp = now;
end

% 计算总速度大小
v_total = sqrt(vr^2 + vaz^2 + vel^2);

% 估计加速度大小（如果有历史数据）
acc_estimate = 0;
if ~isempty(prev_velocities)
    velocity_diff = velocities - prev_velocities;
    acc_estimate = norm(velocity_diff);
    prev_velocities = velocities;
end

% 提取基础噪声参数
q_pos = base_params.pos;       % 位置
q_vel = base_params.vel;       % 速度
q_acc = base_params.acc;       % 加速度
q_angle = base_params.angle;   % 角度
q_ang_vel = base_params.ang_vel;  % 角速度
q_ang_acc = base_params.ang_acc;  % 角加速度
nominal_speed = base_params.nominal_speed;  % 标称速度

% 1. 速度适应因子 - 更平缓的速度与噪声映射
speed_ratio = v_total / nominal_speed;
vel_factor = min(1.5, max(0.9, 0.9 + 0.4 * speed_ratio));  

% 2. 加速度适应因子 - 更平稳的加速度映射
if acc_estimate > 0
    acc_adaptation = min(1.8, 1.0 + acc_estimate / 2.0);
else
    % 无加速度历史数据，仍使用速度作为启发式指标
    acc_adaptation = max(1.0, speed_ratio);
end
% 总加速度因子
acc_factor = min(1.8, acc_adaptation); 

% 3. 角速度适应因子 - 更温和的角度影响
angular_rate_magnitude = sqrt(vaz^2 + vel^2);
if angular_rate_magnitude > 10.0
    % 高角速度情况
    ang_factor = min(1.5, 1.0 + angular_rate_magnitude / 25.0);
elseif angular_rate_magnitude > 2.0
    % 中等角速度情况
    ang_factor = 1.0 + (angular_rate_magnitude - 2.0) / 20.0;
else
    % 低角速度情况
    ang_factor = 0.9;
end

% 创建适应后的过程噪声矩阵
Q = zeros(9, 9);

% 距离方向过程噪声 - 提高距离预测的适应性，避免陷入错误估计
if abs(vr) < 5.0
    % 低速状态，提高距离噪声，让系统更快适应变化
    range_precision_factor = 1.3;
    fprintf('距离高适应性模式: 低速状态\n');
else
    % 高速状态，应用适应因子但同样保持较高噪声
    range_precision_factor = max(1.2, vel_factor);
    fprintf('距离高适应性模式: 高速追踪\n');
end

% 提高距离相关噪声，使系统能够更快纠正错误估计
Q(1,1) = q_pos * range_precision_factor * 1.2;    % 提高距离位置噪声
Q(2,2) = q_vel * vel_factor * 1.1;                % 提高距离速度噪声
Q(3,3) = q_acc * acc_factor;                      % 保持正常加速度噪声

% 方位角方向过程噪声
if abs(vaz) < 3.0
    az_precision_factor = 0.9;
    fprintf('方位角高精度模式\n');
else
    az_precision_factor = ang_factor;
    fprintf('方位角高适应性模式\n');
end

Q(4,4) = q_angle * az_precision_factor;          % 方位角位置噪声
Q(5,5) = q_ang_vel * ang_factor;                 % 方位角速度噪声
Q(6,6) = q_ang_acc * ang_factor;                 % 方位角加速度噪声

% 俯仰角方向过程噪声
if abs(vel) < 3.0
    el_precision_factor = 0.9;
    fprintf('俯仰角高精度模式\n');
else
    el_precision_factor = ang_factor;
    fprintf('俯仰角高适应性模式\n');
end

Q(7,7) = q_angle * el_precision_factor;          % 俯仰角位置噪声
Q(8,8) = q_ang_vel * ang_factor;                 % 俯仰角速度噪声
Q(9,9) = q_ang_acc * ang_factor;                 % 俯仰角加速度噪声

% 减小状态间耦合，防止错误传播
if speed_ratio > 1.2 || acc_factor > 1.5
    coupling_factor = 0.03;  % 减小高动态状态耦合
    fprintf('增强状态间关联度: 高动态状态\n');
else
    coupling_factor = 0.02;  % 进一步减小低动态状态耦合
end

% 距离和角度彻底解耦 - 关键修改：减少距离与角度的关联
for i = 1:3
    for j = 1:3
        if i ~= j
            block_i = (i-1)*3 + (1:3);
            block_j = (j-1)*3 + (1:3);
            
            if (i == 1 && (j == 2 || j == 3)) || ((i == 2 || i == 3) && j == 1)
                % 距离与角度相关块使用更小的耦合因子
                super_small_coupling = 0.01;
                cross_coupling = zeros(3,3);
                for r = 1:3
                    for c = 1:3
                        cross_coupling(r,c) = super_small_coupling * sqrt(Q(block_i(r),block_i(r)) * Q(block_j(c),block_j(c)));
                    end
                end
            else
                % 角度与角度相关块使用正常耦合因子
                cross_coupling = zeros(3,3);
                for r = 1:3
                    for c = 1:3
                        cross_coupling(r,c) = coupling_factor * sqrt(Q(block_i(r),block_i(r)) * Q(block_j(c),block_j(c)));
                    end
                end
            end
            
            % 应用耦合
            Q(block_i, block_j) = cross_coupling;
        end
    end
end

% 确保Q是对称正定的
Q = (Q + Q') / 2;  % 强制对称
min_eig = min(eig(Q));
if min_eig < 1e-6
    Q = Q + 1e-6 * eye(size(Q));  % 确保正定
end

fprintf('过程噪声自适应调整: 速度因子=%.2f, 加速度因子=%.2f, 角度因子=%.2f\n', ...
    vel_factor, acc_factor, ang_factor);

end 