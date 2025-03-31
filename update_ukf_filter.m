function [kf, measurement_innovation] = update_ukf_filter(kf, z, dt, frame_idx)
%UPDATE_UKF_FILTER 更新无迹卡尔曼滤波器
%   kf: UKF滤波器结构
%   z: 测量向量 [r, az, el]
%   dt: 时间步长
%   frame_idx: 当前帧索引

% 更新时间步长
kf.params.dt = dt;

% 更新帧计数
kf.params.frame_count = kf.params.frame_count + 1;

% 检查z是否包含NaN
if any(isnan(z))
    warning('测量值包含NaN，使用仅预测模式');
    % 仅预测，不进行测量更新
    [kf, ~] = predict_only(kf);
    measurement_innovation = NaN(3, 1);
    return;
end

% 存储预更新状态和协方差用于创新计算
x_pre = kf.x;
P_pre = kf.P;

% =================== 状态预测 ===================
% 提取状态变量
r = kf.x(1);     % 距离
vr = kf.x(2);    % 径向速度
ar = kf.x(3);    % 径向加速度
az = kf.x(4);    % 方位角 (度)
vaz = kf.x(5);   % 方位角速度 (度/秒)
aaz = kf.x(6);   % 方位角加速度 (度/秒²)
el = kf.x(7);    % 俯仰角 (度)
vel = kf.x(8);   % 俯仰角速度 (度/秒)
ael = kf.x(9);   % 俯仰角加速度 (度/秒²)

% 执行预测步骤
[kf, ~] = ukf_predict(kf);

% =================== 测量更新 ===================
% 设置测量协方差矩阵 (R) - 考虑距离的动态影响
base_R = kf.params.R;
% 根据目标距离调整测量噪声 - 更远的目标噪声更大
range_factor = min(1 + max(0, (kf.x(1) - 10) / 20), 3.0);
R = base_R .* diag([range_factor, range_factor*1.5, range_factor*1.5]);

% 执行UKF更新步骤
[kf, innovation, ~, nis] = ukf_update(kf, z, R);

% 计算测量创新并记录
measurement_innovation = innovation;

% 在创新序列计算部分确保正确处理角度差异
innovation(2) = wrapTo180(innovation(2));
innovation(3) = wrapTo180(innovation(3));

% 记录NIS值
if ~isfield(kf.params, 'NIS_history')
    kf.params.NIS_history = [];
end
kf.params.NIS_history(end+1) = nis;
if length(kf.params.NIS_history) > 5
    kf.params.NIS_history = kf.params.NIS_history(end-4:end);
end

% 优化自适应机制
if length(kf.params.NIS_history) >= 3
    avg_nis = mean(kf.params.NIS_history);
    
    % NIS持续过高表明滤波器可能发散
    if avg_nis > 15
        fprintf('NIS持续过高(%.2f)，增加过程噪声并调整状态\n', avg_nis);
        
        % 增大过程噪声
        Q_scale = min(4.0, max(2.0, avg_nis/10));
        kf.params.Q = kf.params.initial_Q * Q_scale;
        
        % 对于严重发散情况，直接调整状态
        if avg_nis > 30
            % 记录当前测量
            if ~isfield(kf.params, 'z_history')
                kf.params.z_history = {};
            end
            kf.params.z_history{end+1} = z;
            
            % 如果有足够历史，使用平滑状态校正
            if length(kf.params.z_history) >= 3
                z_recent = kf.params.z_history{end};
                
                % 增加测量对状态的影响
                meas_weight = 0.7;
                kf.x(1) = kf.x(1) * (1-meas_weight) + z_recent(1) * meas_weight;
                
                % 小心处理角度更新
                az_diff = wrapTo180(z_recent(2) - kf.x(4));
                kf.x(4) = wrapTo180(kf.x(4) + az_diff * meas_weight);
                
                el_diff = wrapTo180(z_recent(3) - kf.x(7));
                kf.x(7) = wrapTo180(kf.x(7) + el_diff * meas_weight);
                
                % 减小速度和加速度估计
                vel_scale = 0.3;
                kf.x(2) = kf.x(2) * vel_scale;  % 径向速度
                kf.x(3) = kf.x(3) * vel_scale;  % 径向加速度
                kf.x(5) = kf.x(5) * vel_scale;  % 方位角速度
                kf.x(6) = kf.x(6) * vel_scale;  % 方位角加速度
                kf.x(8) = kf.x(8) * vel_scale;  % 俯仰角速度
                kf.x(9) = kf.x(9) * vel_scale;  % 俯仰角加速度
                
                fprintf('状态校正完成: 新状态 [r=%.2f, az=%.2f, el=%.2f]\n', kf.x(1), kf.x(4), kf.x(7));
            end
            
            % 保持历史记录合理大小
            if length(kf.params.z_history) > 5
                kf.params.z_history = kf.params.z_history(end-4:end);
            end
        end
    elseif avg_nis < 3 && length(kf.params.NIS_history) >= 5
        % NIS持续较低，可以逐渐恢复正常噪声水平
        kf.params.Q = 0.9 * kf.params.Q + 0.1 * kf.params.initial_Q;
    end
end

% 确保物理约束
kf = apply_physical_constraints(kf);

% =================== 自适应机制 ===================
% 检测目标机动和测量异常
if kf.params.adaptive_enabled && kf.params.frame_count > 3
    % 计算创新序列的归一化平方大小
    [kf, modified] = adaptive_filter_update(kf, innovation, nis, z, x_pre);
    
    % 如果自适应更新修改了状态，重新应用约束
    if modified
        kf = apply_physical_constraints(kf);
    end
end

% =================== 稳定性保证 ===================
% 检查协方差矩阵是否为半正定
if kf.params.stabilization_enabled
    kf = ensure_covariance_stability(kf);
end

% 更新状态和测量历史记录
kf.params.state_history{end+1} = kf.x;
kf.params.measurement_history{end+1} = z;
kf.params.innovation_history{end+1} = innovation;

% 修剪历史记录限制长度
if length(kf.params.state_history) > 20
    kf.params.state_history = kf.params.state_history(end-19:end);
    kf.params.measurement_history = kf.params.measurement_history(end-19:end);
    kf.params.innovation_history = kf.params.innovation_history(end-19:end);
end

% 每隔几帧输出跟踪状态
log_interval = 5;
if mod(frame_idx, log_interval) == 0 || frame_idx < 5
    fprintf('帧 %d: 跟踪状态 [r=%.2fm, vr=%.2fm/s, az=%.2f°, el=%.2f°], NIS=%.2f\n', ...
            frame_idx, kf.x(1), kf.x(2), kf.x(4), kf.x(7), nis);
end

% 检查跟踪质量 - 输出警告如果NIS值异常高
last_nis_values = kf.params.NIS_history(max(1, end-2):end);
if mean(last_nis_values) > 10
    fprintf('警告: 连续高NIS值 (%.2f)，可能存在测量-模型不匹配\n', mean(last_nis_values));
    
    % 如果跟踪质量极差，尝试救援措施
    if mean(last_nis_values) > 20 && length(kf.params.measurement_history) > 3
        disp('执行跟踪救援...');
        % 获取最近3次有效测量的平均值
        recent_meas = cat(2, kf.params.measurement_history{max(1,end-2):end});
        avg_meas = mean(recent_meas, 2);
        
        % 偏向最新测量的修正状态
        kf.x(1) = 0.7 * z(1) + 0.3 * kf.x(1);  % 距离偏向测量值
        kf.x(4) = 0.8 * z(2) + 0.2 * kf.x(4);  % 方位角强偏向测量值
        kf.x(7) = 0.8 * z(3) + 0.2 * kf.x(7);  % 俯仰角强偏向测量值
        
        % 降低速度和加速度估计
        kf.x(5) = 0.5 * kf.x(5);  % 减半方位角速度
        kf.x(8) = 0.5 * kf.x(8);  % 减半俯仰角速度
        kf.x(6) = 0.3 * kf.x(6);  % 大幅减小方位角加速度
        kf.x(9) = 0.3 * kf.x(9);  % 大幅减小俯仰角加速度
        
        % 增大协方差表示不确定性增加
        kf.P = kf.P * 1.5;
        fprintf('跟踪重置: 新状态 [r=%.2fm, az=%.2f°, el=%.2f°]\n', kf.x(1), kf.x(4), kf.x(7));
    end
end

end

%% ========== 辅助函数 ==========

function [kf, z_pred] = ukf_predict(kf)
% 无迹卡尔曼滤波器预测步骤

% 提取参数
x = kf.x;
P = kf.P;
n = kf.params.n;
alpha = kf.params.alpha;
beta = kf.params.beta;
kappa = kf.params.kappa;
lambda = kf.params.lambda;
dt = kf.params.dt;
Q = kf.params.Q;

% 计算UKF权重
wm = zeros(2*n+1, 1);  % Sigma点均值权重
wc = zeros(2*n+1, 1);  % Sigma点协方差权重

wm(1) = lambda / (n + lambda);
wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);

for i = 2:2*n+1
    wm(i) = 1 / (2*(n + lambda));
    wc(i) = 1 / (2*(n + lambda));
end

% 构建UKF Sigma点
sigma_points = zeros(n, 2*n+1);
sigma_points(:,1) = x;

% 计算矩阵平方根
% 尝试Cholesky分解，如果失败则使用特征值分解
try
    S = chol(P, 'lower');
catch
    % 如果P不是正定的，使用特征值分解
    [V, D] = eig(P);
    D(D<0) = 0;  % 移除任何负特征值
    S = V * sqrt(D);
end

for i = 1:n
    sigma_points(:,i+1)     = x + sqrt(n + lambda) * S(:,i);
    sigma_points(:,i+n+1)   = x - sqrt(n + lambda) * S(:,i);
end

% 使用动态模型传播Sigma点
propagated_sigma_points = zeros(size(sigma_points));
for i = 1:2*n+1
    propagated_sigma_points(:,i) = state_transition_model(sigma_points(:,i), dt);
end

% 计算预测状态和协方差
x_pred = zeros(n, 1);
for i = 1:2*n+1
    x_pred = x_pred + wm(i) * propagated_sigma_points(:,i);
end

P_pred = Q;  % 初始化为过程噪声协方差
for i = 1:2*n+1
    state_diff = propagated_sigma_points(:,i) - x_pred;
    % 调整角度差异，确保它们在-180到180度之间
    state_diff(4) = wrapTo180(state_diff(4));  % 方位角
    state_diff(7) = wrapTo180(state_diff(7));  % 俯仰角
    
    P_pred = P_pred + wc(i) * (state_diff * state_diff');
end

% 确保协方差矩阵的对称性
P_pred = (P_pred + P_pred')/2;

% 更新滤波器状态
kf.x = x_pred;
kf.P = P_pred;

% 计算预测的测量值
z_pred = kf.x([1,4,7]);  % [r, az, el]
end

function [kf, innovation, Pxz, nis] = ukf_update(kf, z, R)
% 无迹卡尔曼滤波器更新步骤

% 提取参数
x = kf.x;
P = kf.P;
n = kf.params.n;
m = kf.params.m;
alpha = kf.params.alpha;
beta = kf.params.beta;
kappa = kf.params.kappa;
lambda = kf.params.lambda;

% 计算UKF权重
wm = zeros(2*n+1, 1);  % Sigma点均值权重
wc = zeros(2*n+1, 1);  % Sigma点协方差权重

wm(1) = lambda / (n + lambda);
wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);

for i = 2:2*n+1
    wm(i) = 1 / (2*(n + lambda));
    wc(i) = 1 / (2*(n + lambda));
end

% 构建UKF Sigma点
sigma_points = zeros(n, 2*n+1);
sigma_points(:,1) = x;

% 计算矩阵平方根
try
    S = chol(P, 'lower');
catch
    [V, D] = eig(P);
    D(D<0) = 0;  % 移除任何负特征值
    S = V * sqrt(D);
end

for i = 1:n
    sigma_points(:,i+1)     = x + sqrt(n + lambda) * S(:,i);
    sigma_points(:,i+n+1)   = x - sqrt(n + lambda) * S(:,i);
end

% 将Sigma点转换到测量空间 [r, az, el]
meas_sigma_points = zeros(m, 2*n+1);
for i = 1:2*n+1
    meas_sigma_points(:,i) = measurement_model(sigma_points(:,i));
end

% 计算预测测量和测量协方差
z_pred = zeros(m, 1);
for i = 1:2*n+1
    z_pred = z_pred + wm(i) * meas_sigma_points(:,i);
end

% 调整角度测量预测，保证在有效范围内
z_pred(2) = wrapTo180(z_pred(2));  % 方位角
z_pred(3) = wrapTo180(z_pred(3));  % 俯仰角

% 测量协方差
S = R;  % 初始化为测量噪声协方差
for i = 1:2*n+1
    meas_diff = meas_sigma_points(:,i) - z_pred;
    % 调整角度差异
    meas_diff(2) = wrapTo180(meas_diff(2));  % 方位角
    meas_diff(3) = wrapTo180(meas_diff(3));  % 俯仰角
    
    S = S + wc(i) * (meas_diff * meas_diff');
end

% 交叉协方差
Pxz = zeros(n, m);
for i = 1:2*n+1
    state_diff = sigma_points(:,i) - x;
    meas_diff = meas_sigma_points(:,i) - z_pred;
    
    % 调整角度差异
    state_diff(4) = wrapTo180(state_diff(4));  % 状态中的方位角
    state_diff(7) = wrapTo180(state_diff(7));  % 状态中的俯仰角
    meas_diff(2) = wrapTo180(meas_diff(2));    % 测量中的方位角
    meas_diff(3) = wrapTo180(meas_diff(3));    % 测量中的俯仰角
    
    Pxz = Pxz + wc(i) * (state_diff * meas_diff');
end

% 确保测量协方差矩阵为正定
S = (S + S')/2;  % 确保对称性
min_eig = min(eig(S));
if min_eig < 1e-6
    S = S + (1e-6 - min_eig + 1e-10) * eye(size(S));
end

% 计算卡尔曼增益
K = Pxz / S;

% 计算测量创新
innovation = z - z_pred;

% 调整角度创新确保在-180到180度之间
innovation(2) = wrapTo180(innovation(2));  % 方位角
innovation(3) = wrapTo180(innovation(3));  % 俯仰角

% 计算归一化创新平方 (NIS)
nis = innovation' / S * innovation;

% 更新状态和协方差
x_update = x + K * innovation;
P_update = P - K * S * K';

% 确保协方差矩阵对称
P_update = (P_update + P_update')/2;

% 防止角度溢出
x_update(4) = wrapTo180(x_update(4));  % 方位角
x_update(7) = wrapTo180(x_update(7));  % 俯仰角

% 更新滤波器状态
kf.x = x_update;
kf.P = P_update;
end

function [kf, modified] = adaptive_filter_update(kf, innovation, nis, z, x_prev)
% 基于创新的自适应滤波器更新

modified = false;
detection_threshold = kf.params.maneuver_detection_threshold;
innovation_threshold = kf.params.innovation_threshold;

% 检查NIS是否超过阈值
nis_moving_avg = mean(kf.params.NIS_history(max(1, end-3):end));
high_nis = nis_moving_avg > detection_threshold * 3;

% 计算归一化创新
norm_innov = abs(innovation) ./ sqrt(diag(kf.params.R));
high_innov = norm_innov > innovation_threshold;

% 检测角度变化剧烈（可能的快速机动）
angle_changes = [abs(innovation(2)), abs(innovation(3))];
rapid_angle_change = any(angle_changes > 15);  % 超过15度的角度变化

% 如果创新或NIS超过阈值，调整过程噪声
if high_nis || any(high_innov) || rapid_angle_change
    % 增加Q矩阵（过程噪声）以适应可能的机动
    % 根据创新序列的不同元素选择性地增加不同状态的噪声
    
    % 获取原始Q矩阵
    Q_base = kf.params.Q;
    Q_new = Q_base;
    
    % 根据创新大小动态调整放大因子
    innov_factor = max(1.0, min(5.0, norm(innovation) / 3));
    
    if high_innov(1) || rapid_angle_change  % 距离创新大
        % 主要增加距离相关状态的噪声
        Q_new(1:3, 1:3) = Q_base(1:3, 1:3) * (1.0 + innov_factor);
        modified = true;
    end
    
    if high_innov(2) || angle_changes(1) > 10  % 方位角创新大
        % 主要增加方位角相关状态的噪声
        Q_new(4:6, 4:6) = Q_base(4:6, 4:6) * (1.0 + innov_factor * 1.5);
        modified = true;
    end
    
    if high_innov(3) || angle_changes(2) > 10  % 俯仰角创新大
        % 主要增加俯仰角相关状态的噪声
        Q_new(7:9, 7:9) = Q_base(7:9, 7:9) * (1.0 + innov_factor * 1.5);
        modified = true;
    end
    
    if modified
        % 更新Q矩阵
        kf.params.Q = Q_new;
        
        % 记录检测到可能的机动
        fprintf('检测到可能的机动! NIS=%.2f, 归一化创新=%.2f/%.2f/%.2f\n', ...
                nis, norm_innov(1), norm_innov(2), norm_innov(3));
                
        % 如果角度变化特别剧烈，可以考虑调整状态
        if any(angle_changes > 25)  % 超过25度的巨大角度变化
            fprintf('检测到剧烈角度变化: %.2f/%.2f 度，进行状态调整\n', angle_changes(1), angle_changes(2));
            
            % 调整角速度和角加速度估计以适应剧烈变化
            % 方位角变化剧烈
            if angle_changes(1) > 25 
                % 计算角速度，但添加约束
                estimated_az_speed = wrapTo180(z(2) - x_prev(4)) / kf.params.dt;
                kf.x(5) = 0.7 * estimated_az_speed + 0.3 * kf.x(5);  % 70%新估计，30%当前值
                % 角加速度向角速度方向调整
                sign_factor = sign(estimated_az_speed) * sign(kf.x(5));
                if sign_factor >= 0  % 符号一致，保持加速度
                    kf.x(6) = 0.5 * kf.x(6);  % 减小加速度
                else  % 符号不一致，反转加速度
                    kf.x(6) = -0.5 * kf.x(6);  % 反转并减小加速度
                end
            end
            
            % 俯仰角变化剧烈
            if angle_changes(2) > 25
                % 计算角速度，但添加约束
                estimated_el_speed = wrapTo180(z(3) - x_prev(7)) / kf.params.dt;
                kf.x(8) = 0.7 * estimated_el_speed + 0.3 * kf.x(8);  % 70%新估计，30%当前值
                % 角加速度向角速度方向调整
                sign_factor = sign(estimated_el_speed) * sign(kf.x(8));
                if sign_factor >= 0  % 符号一致，保持加速度
                    kf.x(9) = 0.5 * kf.x(9);  % 减小加速度
                else  % 符号不一致，反转加速度
                    kf.x(9) = -0.5 * kf.x(9);  % 反转并减小加速度
                end
            end
        end
    end
else
    % 如果跟踪正常，逐渐恢复正常的过程噪声
    if kf.params.frame_count > 10 && ~isempty(kf.params.NIS_history) && mean(kf.params.NIS_history(end-min(3,length(kf.params.NIS_history))+1:end)) < 1.5
        % 逐渐恢复到基础噪声水平
        kf.params.Q = 0.95 * kf.params.Q + 0.05 * kf.params.initial_Q;
        modified = true;
    end
end

% 保持跟踪协方差矩阵的大小
current_P_trace = trace(kf.P);
trace_ratio = current_P_trace / kf.params.last_P_trace;

% 如果协方差迅速增长，可能是不稳定的迹象
if trace_ratio > 3.0
    fprintf('警告: 协方差迅速增长 (比率=%.2f)，进行调整\n', trace_ratio);
    % 限制协方差增长
    scale_factor = sqrt(3.0 / trace_ratio);
    kf.P = kf.P * scale_factor;
    modified = true;
end

% 更新协方差跟踪记录
kf.params.last_P_trace = trace(kf.P);
end

function x_next = state_transition_model(x, dt)
% 运动模型 - 恒加速度+角度约束模型

% 提取状态
r = x(1);      % 距离
vr = x(2);     % 径向速度
ar = x(3);     % 径向加速度
az = x(4);     % 方位角
vaz = x(5);    % 方位角速度
aaz = x(6);    % 方位角加速度
el = x(7);     % 俯仰角
vel = x(8);    % 俯仰角速度
ael = x(9);    % 俯仰角加速度

% 初始化下一个状态
x_next = zeros(9, 1);

% 使用恒加速度模型 - 距离动力学
x_next(1) = r + vr * dt + 0.5 * ar * dt^2;  % 距离更新
x_next(2) = vr + ar * dt;                  % 径向速度更新
x_next(3) = ar;                            % 径向加速度（常数）

% 角度动力学 - 方位角
x_next(4) = wrapTo180(az + vaz * dt + 0.5 * aaz * dt^2);  % 方位角更新（确保在-180到180度范围内）
x_next(5) = vaz + aaz * dt;                             % 方位角速度更新
x_next(6) = aaz;                                        % 方位角加速度（常数）

% 角度动力学 - 俯仰角
x_next(7) = wrapTo180(el + vel * dt + 0.5 * ael * dt^2);  % 俯仰角更新（确保在-180到180度范围内）
x_next(8) = vel + ael * dt;                             % 俯仰角速度更新
x_next(9) = ael;                                        % 俯仰角加速度（常数）

end

function z = measurement_model(x)
% 测量模型 - 从状态到测量的转换

% 测量向量 [r, az, el]
z = zeros(3, 1);
z(1) = x(1);             % 距离
z(2) = wrapTo180(x(4));  % 方位角，确保在-180到180范围内
z(3) = wrapTo180(x(7));  % 俯仰角，确保在-180到180范围内
end

function [kf, z_pred] = predict_only(kf)
    % 直接实现预测步骤，避免递归调用ukf_predict
    dt = kf.params.dt;
    x = kf.x;
    P = kf.P;
    n = kf.params.n;
    alpha = kf.params.alpha;
    beta = kf.params.beta;
    kappa = kf.params.kappa;
    Q = kf.params.Q;
    
    % 计算lambda和权重
    lambda = alpha^2 * (n + kappa) - n;
    
    % 计算权重
    wm = zeros(2*n+1, 1); % 均值权重
    wc = zeros(2*n+1, 1); % 协方差权重
    
    wm(1) = lambda / (n + lambda);
    wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
    for i = 2:2*n+1
        wm(i) = 1 / (2*(n + lambda));
        wc(i) = 1 / (2*(n + lambda));
    end
    
    % 确保P是对称的
    P = (P + P')/2;
    
    % 计算矩阵平方根
    try
        % 尝试Cholesky分解
        S = chol(P, 'lower');
    catch
        % 如果失败，使用SVD
        [U, D, ~] = svd(P);
        S = U * sqrt(D);
    end
    
    % 生成sigma点
    X = zeros(n, 2*n+1);
    X(:,1) = x;
    
    for i = 1:n
        X(:,i+1) = x + sqrt(n + lambda) * S(:,i);
        X(:,i+n+1) = x - sqrt(n + lambda) * S(:,i);
    end
    
    % 传播sigma点
    Y = zeros(n, 2*n+1);
    for i = 1:2*n+1
        % 提取当前sigma点
        xi = X(:,i);
        
        % 提取状态变量
        r = xi(1);      % 距离
        vr = xi(2);     % 径向速度
        ar = xi(3);     % 径向加速度
        az = xi(4);     % 方位角
        vaz = xi(5);    % 方位角速度
        aaz = xi(6);    % 方位角加速度
        el = xi(7);     % 俯仰角
        vel = xi(8);    % 俯仰角速度
        ael = xi(9);    % 俯仰角加速度
        
        % 应用物理约束
        max_speed = 40;
        max_accel = 20;
        max_ang_speed = 30;
        max_ang_accel = 15;
        
        % 限制速度和加速度
        vr = sign(vr) * min(abs(vr), max_speed);
        ar = sign(ar) * min(abs(ar), max_accel);
        vaz = sign(vaz) * min(abs(vaz), max_ang_speed);
        aaz = sign(aaz) * min(abs(aaz), max_ang_accel);
        vel = sign(vel) * min(abs(vel), max_ang_speed);
        ael = sign(ael) * min(abs(ael), max_ang_accel);
        
        % 恒加速度模型更新
        r_new = max(0.1, r + vr*dt + 0.5*ar*dt^2);
        vr_new = vr + ar*dt;
        ar_new = ar;
        
        az_new = az + vaz*dt + 0.5*aaz*dt^2;
        vaz_new = vaz + aaz*dt;
        aaz_new = aaz;
        
        el_new = el + vel*dt + 0.5*ael*dt^2;
        vel_new = vel + ael*dt;
        ael_new = ael;
        
        % 角度归一化
        az_new = wrapTo180(az_new);
        el_new = wrapTo180(el_new);
        
        % 更新传播后的sigma点
        Y(:,i) = [r_new; vr_new; ar_new; az_new; vaz_new; aaz_new; el_new; vel_new; ael_new];
    end
    
    % 计算预测状态均值
    x_pred = zeros(n, 1);
    for i = 1:2*n+1
        x_pred = x_pred + wm(i) * Y(:,i);
    end
    
    % 角度归一化
    x_pred(4) = wrapTo180(x_pred(4));
    x_pred(7) = wrapTo180(x_pred(7));
    
    % 计算预测协方差
    P_pred = zeros(n, n);
    for i = 1:2*n+1
        diff = Y(:,i) - x_pred;
        
        % 角度差异特殊处理
        diff(4) = wrapTo180(diff(4));
        diff(7) = wrapTo180(diff(7));
        
        P_pred = P_pred + wc(i) * (diff * diff');
    end
    
    % 添加过程噪声
    P_pred = P_pred + Q;
    
    % 确保协方差矩阵对称
    P_pred = (P_pred + P_pred')/2;
    
    % 更新滤波器状态
    kf.x = x_pred;
    kf.P = P_pred;
    
    % 预测测量值
    z_pred = [x_pred(1); x_pred(4); x_pred(7)];
end

function kf = apply_physical_constraints(kf)
    % 提取状态
    x = kf.x;
    
    % 定义约束
    min_range = 0.1;        % 最小距离(m)
    max_range = 500.0;      % 最大距离(m)
    max_speed = 30.0;       % 最大速度(m/s)
    max_accel = 20.0;       % 最大加速度(m/s²)
    max_ang_speed = 20.0;   % 最大角速度(度/s)
    max_ang_accel = 10.0;   % 最大角加速度(度/s²)
    
    % 应用距离约束
    x(1) = max(min_range, min(max_range, x(1)));
    
    % 应用速度约束
    x(2) = sign(x(2)) * min(abs(x(2)), max_speed);
    x(5) = sign(x(5)) * min(abs(x(5)), max_ang_speed);
    x(8) = sign(x(8)) * min(abs(x(8)), max_ang_speed);
    
    % 应用加速度约束
    x(3) = sign(x(3)) * min(abs(x(3)), max_accel);
    x(6) = sign(x(6)) * min(abs(x(6)), max_ang_accel);
    x(9) = sign(x(9)) * min(abs(x(9)), max_ang_accel);
    
    % 确保角度一致性
    x(4) = wrapTo180(x(4));
    x(7) = wrapTo180(x(7));
    
    % 更新状态
    kf.x = x;
    
    % 确保协方差矩阵是对称正定的
    P = kf.P;
    P = (P + P')/2;
    
    % 检查特征值
    [V, D] = eig(P);
    D = diag(D);
    if any(D < 1e-6)
        D(D < 1e-6) = 1e-6;
        P = V * diag(D) * V';
        P = (P + P')/2;
    end
    
    kf.P = P;
end

function kf = ensure_covariance_stability(kf)
% 确保协方差矩阵的数值稳定性

P = kf.P;
min_eig_threshold = kf.params.min_eig_threshold;
diag_loading_factor = kf.params.diag_loading_factor;

% 确保协方差矩阵对称
P = (P + P')/2;

% 检查最小特征值
min_eig = min(eig(P));
if min_eig < min_eig_threshold
    % 应用对角加载以确保半正定性
    P = P + (min_eig_threshold - min_eig + diag_loading_factor) * eye(size(P));
end

% 更新协方差矩阵
kf.P = P;
end 