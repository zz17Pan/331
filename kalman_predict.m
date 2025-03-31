function kf = kalman_predict(kf, dt, params)
%KALMAN_PREDICT 卡尔曼滤波器预测步骤
%   kf: 卡尔曼滤波器状态结构体
%   dt: 时间步长
%   params: 系统参数

% 提取参数
max_acceleration = params.sim.max_acceleration;

% 创建状态转移矩阵
F = eye(9);

% 提取当前状态
r = kf.x(1);        % 距离
vr = kf.x(2);       % 径向速度
ar = kf.x(3);       % 径向加速度
az = kf.x(4);       % 方位角
vaz = kf.x(5);      % 方位角速度
aaz = kf.x(6);      % 方位角加速度
el = kf.x(7);       % 俯仰角
vel = kf.x(8);      % 俯仰角速度
ael = kf.x(9);      % 俯仰角加速度

% 静态变量保存历史状态用于机动检测
persistent prev_states;
persistent frame_counter;
persistent maneuver_status;

% 初始化静态变量
if isempty(prev_states)
    prev_states = zeros(9, 5);  % 存储最近5帧状态
    frame_counter = 1;
    maneuver_status = 0;  % 0-无机动, 1-轻微机动, 2-大幅机动
else
    frame_counter = frame_counter + 1;
end

% 更新历史状态队列
prev_states = [prev_states(:, 2:end), kf.x];

% 机动检测 - 增强检测逻辑，避免误判
maneuver_detected = false;
maneuver_level = 0;

% 只有有足够历史数据才进行检测
if frame_counter > 3
    % 1. 加速度幅值检测
    current_acc = sqrt(ar^2 + (r*aaz*cosd(el))^2 + (r*ael)^2);
    if current_acc > 2.0  % 适当提高阈值，减少误报
        maneuver_level = max(maneuver_level, 1);
        if current_acc > 4.0
            maneuver_level = 2;
        end
    end
    
    % 2. 轨迹弯曲度检测(如果历史足够长)
    if frame_counter > 4
        % 取最近的几帧位置估计计算轨迹弯曲度
        positions = zeros(3, min(5, frame_counter));
        for i = 1:size(positions, 2)
            r_i = prev_states(1, end-size(positions, 2)+i);
            az_i = prev_states(4, end-size(positions, 2)+i);
            el_i = prev_states(7, end-size(positions, 2)+i);
            
            % 转换为笛卡尔坐标
            x_i = r_i * cosd(el_i) * cosd(az_i);
            y_i = r_i * cosd(el_i) * sind(az_i);
            z_i = r_i * sind(el_i);
            
            positions(:, i) = [x_i; y_i; z_i];
        end
        
        % 计算轨迹曲率度量
        if size(positions, 2) >= 3
            vec1 = positions(:, 2) - positions(:, 1);
            vec2 = positions(:, 3) - positions(:, 2);
            
            % 归一化向量
            vec1 = vec1 / norm(vec1);
            vec2 = vec2 / norm(vec2);
            
            % 通过点积计算角度变化
            cos_angle = dot(vec1, vec2);
            cos_angle = min(1.0, max(-1.0, cos_angle));  % 确保在[-1,1]范围内
            angle_change = acosd(cos_angle);
            
            % 显著曲率表示机动
            if angle_change > 12.0  % 适当提高阈值
                maneuver_level = max(maneuver_level, 1);
                if angle_change > 20.0
                    maneuver_level = 2;
                end
            end
        end
    end
    
    % 3. 速度变化检测
    if frame_counter > 3
        prev_v = sqrt(prev_states(2, end-1)^2 + ...
                     (prev_states(1, end-1) * prev_states(5, end-1) * cosd(prev_states(7, end-1)))^2 + ...
                     (prev_states(1, end-1) * prev_states(8, end-1))^2);
        current_v = sqrt(vr^2 + (r*vaz*cosd(el))^2 + (r*vel)^2);
        
        v_change = abs(current_v - prev_v) / dt;
        if v_change > 2.5  % 适当提高阈值，减少误报
            maneuver_level = max(maneuver_level, 1);
            if v_change > 5.0
                maneuver_level = 2;
            end
        end
    end
    
    % 平滑机动状态，避免状态频繁跳变
    maneuver_status = round(0.7 * maneuver_status + 0.3 * maneuver_level);
    
    % 记录机动状态
    if maneuver_status > 0
        maneuver_detected = true;
        if maneuver_status == 1
            fprintf('轻微机动状态检测\n');
        else
            fprintf('大幅机动状态检测\n');
        end
    end
end

% 状态转移矩阵 F 构建
% 距离相关状态转移
F(1, 1) = 1;
F(1, 2) = dt;
F(1, 3) = 0.5 * dt^2;

F(2, 2) = 1;
F(2, 3) = dt;

F(3, 3) = 1;

% 方位角相关状态转移
F(4, 4) = 1;
F(4, 5) = dt;
F(4, 6) = 0.5 * dt^2;

F(5, 5) = 1;
F(5, 6) = dt;

F(6, 6) = 1;

% 俯仰角相关状态转移
F(7, 7) = 1;
F(7, 8) = dt;
F(7, 9) = 0.5 * dt^2;

F(8, 8) = 1;
F(8, 9) = dt;

F(9, 9) = 1;

% 基于机动状态调整过程噪声
if maneuver_detected
    % 增强的机动处理
    if maneuver_status == 2
        % 大幅机动状态下，增大过程噪声以适应快速变化
        alpha = 1.8;  % 均衡因子，不要过大
        beta = 0.6;   % 均衡因子，不要过小
    else
        % 轻微机动状态下，轻微增大过程噪声
        alpha = 1.4;
        beta = 0.8;
    end
else
    % 正常模式，保持默认状态转移
    alpha = 1.0;
    beta = 1.0;
end

% 应用状态转移矩阵
kf.x = F * kf.x;

% 确保更新后的距离和速度方向一致 - 关键修改
if frame_counter <= 3  % 只在初始几帧强制保持一致性
    if sign(kf.x(2)) ~= sign(vr) && abs(kf.x(2)) > 0.5
        fprintf('前3帧强制保持速度方向一致性: %.2f -> %.2f\n', kf.x(2), vr);
        kf.x(2) = vr;  % 保持速度方向一致性
    end
else
    if sign(kf.x(2)) ~= sign(vr) && abs(kf.x(2)) > 8.0 && abs(vr) > 8.0
        % 只有当速度反向且幅值都很大时才输出警告，但不强制修正
        fprintf('警告: 预测速度方向(%.2f)与先验值(%.2f)不一致\n', kf.x(2), vr);
    end
end

% 更新状态协方差矩阵 P
kf.P = F * kf.P * F' + kf.Q;

% 确保 P 是对称的
kf.P = (kf.P + kf.P') / 2;

% 应用状态约束（如果需要）
if isfield(params, 'constrain_state') && params.constrain_state
    kf = apply_state_constraints(kf, params);
end

end 