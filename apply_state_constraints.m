function x_constrained = apply_state_constraints(x, params)
%APPLY_STATE_CONSTRAINTS 应用物理约束到状态向量
%   x: 状态向量
%   params: 系统参数
%   x_constrained: 应用约束后的状态向量

% 复制状态向量
x_constrained = x;

% 距离约束 - 不能为负或过小
if x_constrained(1) < 0.1
    fprintf('应用约束: 距离 %.2f -> 0.1\n', x_constrained(1));
    x_constrained(1) = 0.1;
end

% 最大距离约束 (可选，取决于应用场景)
max_distance = 1000.0;  % 假设最大距离为1000米
if x_constrained(1) > max_distance
    fprintf('应用约束: 距离 %.2f -> %.2f\n', x_constrained(1), max_distance);
    x_constrained(1) = max_distance;
end

% 速度约束
max_speed = 40.0;  % 最大速度 m/s
if abs(x_constrained(2)) > max_speed
    fprintf('应用约束: 径向速度 %.2f -> %.2f\n', x_constrained(2), sign(x_constrained(2))*max_speed);
    x_constrained(2) = sign(x_constrained(2)) * max_speed;
end

% 加速度约束
max_accel = 20.0;  % 最大加速度 m/s^2
if abs(x_constrained(3)) > max_accel
    fprintf('应用约束: 径向加速度 %.2f -> %.2f\n', x_constrained(3), sign(x_constrained(3))*max_accel);
    x_constrained(3) = sign(x_constrained(3)) * max_accel;
end

% 角度约束 - 确保在[-180, 180]范围内
x_constrained(4) = wrapTo180(x_constrained(4));  % 方位角
x_constrained(7) = wrapTo180(x_constrained(7));  % 俯仰角

% 角速度约束
max_angular_vel = 30.0;  % 最大角速度 度/秒
if abs(x_constrained(5)) > max_angular_vel
    fprintf('应用约束: 方位角速度 %.2f -> %.2f\n', x_constrained(5), sign(x_constrained(5))*max_angular_vel);
    x_constrained(5) = sign(x_constrained(5)) * max_angular_vel;
end
if abs(x_constrained(8)) > max_angular_vel
    fprintf('应用约束: 俯仰角速度 %.2f -> %.2f\n', x_constrained(8), sign(x_constrained(8))*max_angular_vel);
    x_constrained(8) = sign(x_constrained(8)) * max_angular_vel;
end

% 角加速度约束
max_angular_accel = 10.0;  % 最大角加速度 度/秒^2
if abs(x_constrained(6)) > max_angular_accel
    fprintf('应用约束: 方位角加速度 %.2f -> %.2f\n', x_constrained(6), sign(x_constrained(6))*max_angular_accel);
    x_constrained(6) = sign(x_constrained(6)) * max_angular_accel;
end
if abs(x_constrained(9)) > max_angular_accel
    fprintf('应用约束: 俯仰角加速度 %.2f -> %.2f\n', x_constrained(9), sign(x_constrained(9))*max_angular_accel);
    x_constrained(9) = sign(x_constrained(9)) * max_angular_accel;
end

% 运动学一致性检查 (可选)
% 如果有明确的物理约束，可以在这里添加

return;
end
