function maneuver_detected = detect_maneuver(ukf)
%DETECT_MANEUVER 检测目标是否发生机动
%   ukf: UKF滤波器状态
%   maneuver_detected: 是否检测到机动

% 获取最近状态历史
states = ukf.history.states;
pointer = ukf.history.pointer;

% 如果历史不足，无法判断
if pointer < 4
    maneuver_detected = false;
    return;
end

% 获取最近几帧的状态
recent_indices = mod((pointer-4:pointer-1)-1, 10) + 1;
recent_states = states(:, recent_indices);

% 计算加速度变化
accel_changes = [];
for i = 2:size(recent_states, 2)
    % 径向加速度变化
    accel_changes(end+1) = abs(recent_states(3,i) - recent_states(3,i-1));
    % 方位角加速度变化
    accel_changes(end+1) = abs(recent_states(6,i) - recent_states(6,i-1));
    % 俯仰角加速度变化
    accel_changes(end+1) = abs(recent_states(9,i) - recent_states(9,i-1));
end

% 判断机动
accel_threshold = 2.0;
if any(accel_changes > accel_threshold)
    maneuver_detected = true;
else
    % 检查速度变化趋势
    vel_changes = [];
    for i = 2:size(recent_states, 2)
        % 径向速度变化
        vel_changes(end+1) = abs(recent_states(2,i) - recent_states(2,i-1));
        % 方位角速度变化
        vel_changes(end+1) = abs(recent_states(5,i) - recent_states(5,i-1));
        % 俯仰角速度变化
        vel_changes(end+1) = abs(recent_states(8,i) - recent_states(8,i-1));
    end
    
    vel_threshold = 1.5;
    maneuver_detected = any(vel_changes > vel_threshold);
end
end 