function kf = kalman_update(kf, z, params)
%KALMAN_UPDATE 执行卡尔曼滤波器的更新步骤，增强对机动目标的跟踪能力
%   kf: 卡尔曼滤波器状态结构体
%   z: 测量向量 [距离; 方位角; 俯仰角]
%   params: 系统参数结构体
%   更新后的卡尔曼滤波器状态

% 保存更新前的状态用于后续分析
prev_x = kf.x;

% 创建残差保存变量
kf.residual = zeros(3, 1);

% 计算预测的测量值
h_x = [kf.x(1); kf.x(4); kf.x(7)];  % 预测的测量 [r; az; el]

% 检查测量值与预测值是否存在明显差异
range_diff = z(1) - h_x(1);
az_diff = wrapTo180(z(2) - h_x(2));
el_diff = wrapTo180(z(3) - h_x(3));

% 保存历史测量数据用于平滑和异常值检测
persistent z_history;
if isempty(z_history)
    z_history = repmat(z, 1, 5);
else
    z_history = [z_history(:,2:end), z];
end

% 检测并修正异常测量值 - 使用中值滤波器思想
if size(z_history, 2) >= 5
    % 检查距离测量
    sorted_ranges = sort(z_history(1,:));
    median_range = sorted_ranges(3);  % 中位数
    
    % 如果当前测量偏离中位数过大，用中值加权替换
    if abs(z(1) - median_range) > 3.0
        weight = 0.3;  % 仅给异常值很小的权重
        z(1) = weight * z(1) + (1-weight) * median_range;
        fprintf('检测到距离测量异常，应用中值滤波修正: %.2f -> %.2f\n', z_history(1,end), z(1));
    end
    
    % 检查角度测量
    sorted_az = sort(z_history(2,:));
    median_az = sorted_az(3);
    if abs(wrapTo180(z(2) - median_az)) > 2.0
        weight = 0.4;
        z(2) = wrapTo180(weight * z(2) + (1-weight) * median_az);
        fprintf('检测到方位角测量异常，应用中值滤波修正: %.2f -> %.2f\n', z_history(2,end), z(2));
    end
    
    sorted_el = sort(z_history(3,:));
    median_el = sorted_el(3);
    if abs(wrapTo180(z(3) - median_el)) > 2.0
        weight = 0.4;
        z(3) = wrapTo180(weight * z(3) + (1-weight) * median_el);
        fprintf('检测到俯仰角测量异常，应用中值滤波修正: %.2f -> %.2f\n', z_history(3,end), z(3));
    end
end

% 特别处理测量噪声异常大的情况
max_acceptable_range_diff = 10.0; % 降低最大可接受的距离差异
if abs(range_diff) > max_acceptable_range_diff
    fprintf('距离测量异常大: 测量=%.2f, 预测=%.2f, 差异=%.2f\n', z(1), h_x(1), range_diff);
    
    % 应用渐进式测量矫正，防止单个大测量值完全破坏估计
    correction_ratio = max_acceptable_range_diff / abs(range_diff);
    z(1) = h_x(1) + range_diff * correction_ratio;
    fprintf('应用渐进测量矫正: 修正测量距离为 %.2f\n', z(1));
    
    % 同时暂时增大测量噪声
    kf.R(1,1) = kf.R(1,1) * 2.0;
end

% 特别处理角度测量异常的情况
max_acceptable_angle_diff = 5.0; % 角度允许最大差异(度)
if abs(az_diff) > max_acceptable_angle_diff
    correction_ratio = max_acceptable_angle_diff / abs(az_diff);
    z(2) = h_x(2) + az_diff * correction_ratio;
    fprintf('方位角测量异常大，应用渐进矫正: %.2f -> %.2f\n', h_x(2) + az_diff, z(2));
    kf.R(2,2) = kf.R(2,2) * 1.5;
end

if abs(el_diff) > max_acceptable_angle_diff
    correction_ratio = max_acceptable_angle_diff / abs(el_diff);
    z(3) = h_x(3) + el_diff * correction_ratio;
    fprintf('俯仰角测量异常大，应用渐进矫正: %.2f -> %.2f\n', h_x(3) + el_diff, z(3));
    kf.R(3,3) = kf.R(3,3) * 1.5;
end

% 计算创新序列（测量残差）
y = z - h_x;

% 保存残差供其他函数使用
kf.residual = y;

% 检查距离测量的可靠性
% 如果多次连续观测到距离增加，但估计的速度为负（或反之），可能存在问题
persistent range_trend_counter;
if isempty(range_trend_counter)
    range_trend_counter = 0;
end

% 如果测量距离与预测距离变化方向相反，且速度与变化不一致
if range_diff * kf.x(2) < 0 && abs(range_diff) > 1.5  % 降低阈值，提高敏感度
    range_trend_counter = range_trend_counter + 1;
    if range_trend_counter >= 2  % 降低触发计数器阈值
        fprintf('警告: 连续%d帧检测到距离变化趋势与速度不一致!\n', range_trend_counter);
        
        % 如果多次连续不一致，调整速度符号及大小
        if abs(kf.x(2)) > 1.5
            old_vel = kf.x(2);
            kf.x(2) = sign(range_diff) * abs(kf.x(2)) * 0.9;
            fprintf('修正速度方向: %.2f -> %.2f\n', old_vel, kf.x(2));
            range_trend_counter = 0; % 重置计数器
        end
    end
else
    range_trend_counter = max(0, range_trend_counter - 1); %逐渐减少计数器
end

% 处理角度的周期性 (限制在±180度范围内)
y(2) = wrapTo180(y(2));  % 方位角
y(3) = wrapTo180(y(3));  % 俯仰角

% 动态调整测量噪声协方差
kf.R = adjust_measurement_noise(y, kf.params.base_r);

% 附加精度增强：针对目前状态应用特定优化
% 1. 针对不同距离区间采用不同的测量噪声模型
if kf.x(1) < 30.0
    % 近距离，角度测量更精确，降低角度噪声
    kf.R(2,2) = kf.R(2,2) * 0.7;
    kf.R(3,3) = kf.R(3,3) * 0.7;
    fprintf('近距离目标: 减小角度测量噪声\n');
elseif kf.x(1) > 100.0
    % 远距离，距离测量更精确，降低距离噪声
    kf.R(1,1) = kf.R(1,1) * 0.8;
    fprintf('远距离目标: 减小距离测量噪声\n');
end

% 2. 根据速度和加速度调整噪声模型
if abs(kf.x(2)) > 10.0 || abs(kf.x(3)) > 3.0
    % 高速或高加速度情况下，增大距离测量噪声，减小角度测量噪声
    % 因为高速运动下角度测量通常更可靠
    kf.R(1,1) = kf.R(1,1) * 1.2;
    kf.R(2,2) = kf.R(2,2) * 0.9;
    kf.R(3,3) = kf.R(3,3) * 0.9;
    fprintf('高速/加速状态: 调整测量噪声模型\n');
end

% 计算创新协方差
S = kf.H * kf.P * kf.H' + kf.R;

% 确保S矩阵正定性
S = 0.5 * (S + S');  % 确保对称
[V, D] = eig(S);    % 特征值分解
min_eigenvalue = 1e-6;
for i = 1:size(D, 1)
    if D(i, i) < min_eigenvalue
        D(i, i) = min_eigenvalue;
    end
end
S = V * D * V';  % 重建协方差矩阵

% 计算归一化创新平方 (NIS) - 用于检测机动
nis = y' / S * y;

% 更新机动检测器
maneuver_detector = kf.params.maneuver_detector;
window_size = maneuver_detector.window_size;

% 更新创新历史
maneuver_detector.history = [maneuver_detector.history(2:end), nis];

% 计算平均创新序列
mean_nis = mean(maneuver_detector.history);

% 检测机动
prev_maneuver_state = maneuver_detector.is_maneuvering;
if mean_nis > maneuver_detector.threshold
    maneuver_detector.is_maneuvering = true;
    if ~prev_maneuver_state
        fprintf('检测到目标开始机动！增强跟踪参数\n');
    end
else
    % 如果连续多帧NIS较低，退出机动状态
    if mean_nis < maneuver_detector.threshold / 2 && maneuver_detector.is_maneuvering
        maneuver_detector.is_maneuvering = false;
        fprintf('目标机动结束，恢复标准跟踪参数\n');
    end
end

% 更新检测器
kf.params.maneuver_detector = maneuver_detector;

% 自适应增益因子 - 根据机动状态调整
beta = 1.0;  % 默认增益因子

% 根据机动状态和创新大小调整增益因子
if maneuver_detector.is_maneuvering
    % 机动状态下，使用互变模型滤波(IMM)思想，提高对状态变化的敏感度
    beta = sqrt(max(1.0, min(3.0, nis / kf.params.chi2_threshold)));
    fprintf('目标正在机动，跟踪增益因子 beta=%.2f\n', beta);
else
    % 非机动状态，保守跟踪
    if nis > kf.params.chi2_threshold
        % 异常数据点的有限响应
        beta = sqrt(min(1.5, nis / kf.params.chi2_threshold));
    end
end

% 对非常大的创新应用渐进跟踪，避免滤波器被异常值带偏
innovation_magnitude = norm(y);
max_innovation = 5.0;  % 进一步降低最大可接受的创新大小
if innovation_magnitude > max_innovation
    progression_factor = max_innovation / innovation_magnitude;
    y = y * progression_factor;
    fprintf('应用渐进跟踪: 创新大小 %.2f 减小到 %.2f\n', innovation_magnitude, max_innovation);
end

% 精确的卡尔曼增益计算
% 标准卡尔曼增益
K_standard = kf.P * kf.H' / S;

% 针对机动/非机动状态的差异化处理
K = K_standard;

% 差异化增益策略 - 根据残差分量进行优化
% 1. 如果距离残差小但角度残差大，增强角度分量的更新
if abs(y(1)) < 1.0 && (abs(y(2)) > 1.0 || abs(y(3)) > 1.0)
    angle_enhance = 1.3;
    K(:,2) = K(:,2) * angle_enhance;  % 增强方位角影响
    K(:,3) = K(:,3) * angle_enhance;  % 增强俯仰角影响
    fprintf('增强角度更新: 距离残差小但角度残差大\n');
end

% 2. 如果角度残差小但距离残差大，增强距离分量的更新
if (abs(y(2)) < 0.5 && abs(y(3)) < 0.5) && abs(y(1)) > 2.0
    range_enhance = 1.3;
    K(:,1) = K(:,1) * range_enhance;  % 增强距离影响
    fprintf('增强距离更新: 角度残差小但距离残差大\n');
end

% 根据当前估计状态针对性调整增益
% 对加速度分量增强更新
if maneuver_detector.is_maneuvering
    % 增强加速度分量的更新敏感度
    acc_enhance_factor = 1.5; 
    K(3,:) = K(3,:) * acc_enhance_factor;  % 距离加速度
    K(6,:) = K(6,:) * acc_enhance_factor;  % 方位角加速度
    K(9,:) = K(9,:) * acc_enhance_factor;  % 俯仰角加速度
end

% 增强距离分量的更新，特别是当距离误差明显存在时
if range_trend_counter >= 2
    % 如果连续多帧检测到距离变化与速度不一致，增强速度的更新
    K(2,1) = K(2,1) * 1.6; 
    fprintf('增强速度对距离测量的响应\n');
end

% 针对高精度距离和角度要求的特殊处理
% 对于距离误差，确保有足够的校正能力
if abs(y(1)) > 1.0  % 如果距离误差超过1米
    dist_correction_enhance = min(1.5, 1.0 + abs(y(1))/10);  
    K(1,1) = K(1,1) * dist_correction_enhance;
    fprintf('增强距离校正: 系数%.2f\n', dist_correction_enhance);
end

% 对于角度误差，确保有足够的校正能力
if abs(y(2)) > 0.8  % 方位角误差接近1度
    az_correction_enhance = min(1.8, 1.0 + abs(y(2))/5);
    K(4,2) = K(4,2) * az_correction_enhance;
    fprintf('增强方位角校正: 系数%.2f\n', az_correction_enhance);
end

if abs(y(3)) > 0.8  % 俯仰角误差接近1度
    el_correction_enhance = min(1.8, 1.0 + abs(y(3))/5);
    K(7,3) = K(7,3) * el_correction_enhance;
    fprintf('增强俯仰角校正: 系数%.2f\n', el_correction_enhance);
end

% 状态更新
kf.x = kf.x + K * y;

% 协方差更新，使用Joseph公式以增强数值稳定性
I = eye(size(K, 1));
% 创建9x9的完整测量噪声矩阵
full_R = zeros(9, 9);
% 将3x3测量噪声设置到对应位置
indices = [1, 4, 7];  % 对应的状态索引(距离、方位角、俯仰角)
for i = 1:3
    for j = 1:3
        full_R(indices(i), indices(j)) = kf.R(i, j);
    end
end
% 使用修正后的协方差矩阵
kf.P = (I - K * kf.H) * kf.P * (I - K * kf.H)' + K * (kf.H * full_R * kf.H') * K';

% 应用蒙特卡洛变换修正协方差估计 - 防止协方差过度收敛
if maneuver_detector.is_maneuvering
    % 机动状态下增加滤波的不确定性
    cov_inflation_factor = 1.2; 
    kf.P = kf.P * cov_inflation_factor;
    
    % 特别增加加速度的协方差
    acc_indices = [3, 6, 9];  % 加速度状态的索引
    acc_inflation = 1.4; 
    for i = acc_indices
        kf.P(i,i) = kf.P(i,i) * acc_inflation;
    end
else
    % 非机动状态下逐渐减小不确定性，提高稳定精度
    stable_deflation = 0.98;  
    kf.P = kf.P * stable_deflation;
end

% 保证协方差矩阵的正定性
% 对角线元素如果太小，设为最小阈值
min_variance = 1e-6;
for i = 1:size(kf.P, 1)
    if kf.P(i,i) < min_variance
        kf.P(i,i) = min_variance;
    end
end

% 检查更新后的状态是否合理
vel_change = kf.x(2) - prev_x(2);
if abs(vel_change) > 3.0  % 降低速度变化阈值
    fprintf('速度变化较大: %.2f -> %.2f (变化: %.2f)\n', prev_x(2), kf.x(2), vel_change);
    
    % 限制单步速度变化幅度
    max_vel_change = 3.0;
    if abs(vel_change) > max_vel_change
        kf.x(2) = prev_x(2) + sign(vel_change) * max_vel_change;
        fprintf('限制速度变化幅度: %.2f\n', kf.x(2));
    end
end

% 平滑位置估计，减少抖动
% 对距离的平滑仅在残差小的情况下应用，防止滞后
if abs(y(1)) < 1.0
    persistent smooth_r;
    if isempty(smooth_r)
        smooth_r = kf.x(1);
    end
    
    smooth_factor = 0.7;  % 使用更高的平滑因子
    smooth_r = smooth_factor * smooth_r + (1-smooth_factor) * kf.x(1);
    kf.x(1) = smooth_r;
end

% 角度平滑，只在残差小的情况下应用
if abs(y(2)) < 1.0
    persistent smooth_az;
    if isempty(smooth_az)
        smooth_az = kf.x(4);
    end
    
    az_smooth_factor = 0.7;
    smooth_az = wrapTo180(az_smooth_factor * smooth_az + (1-az_smooth_factor) * kf.x(4));
    kf.x(4) = smooth_az;
end

if abs(y(3)) < 1.0
    persistent smooth_el;
    if isempty(smooth_el)
        smooth_el = kf.x(7);
    end
    
    el_smooth_factor = 0.7;
    smooth_el = wrapTo180(el_smooth_factor * smooth_el + (1-el_smooth_factor) * kf.x(7));
    kf.x(7) = smooth_el;
end

% 应用状态约束，确保状态变量在物理有效范围内
% 保存Kalman增益
kf.K = K;

% 应用约束前记录一下状态
pre_constrained_x = kf.x;

% 应用约束
kf = apply_state_constraints(kf, kf.params);

% 输出创新信息
fprintf('创新序列: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°, NIS=%.2f\n', ...
    y(1), y(2), y(3), nis);

% 输出当前估计状态
fprintf('状态估计: r=%.2f m, vr=%.2f m/s, ar=%.2f m/s²\n', ...
    kf.x(1), kf.x(2), kf.x(3));
fprintf('角度估计: az=%.2f°(v=%.2f°/s, a=%.2f°/s²), el=%.2f°(v=%.2f°/s, a=%.2f°/s²)\n', ...
    kf.x(4), kf.x(5), kf.x(6), kf.x(7), kf.x(8), kf.x(9));

end