%% 雷达目标跟踪演示
% 演示使用稀疏恢复和无迹卡尔曼滤波器进行3D目标跟踪
clc;
clear;
close all;

% 添加路径
addpath('utils');
addpath('data');

% 设置随机数种子，确保结果可重复
rng(42);

%% 加载数据
disp('加载雷达数据...');
radar_data = load('radar_data.mat');
num_frames = length(radar_data.radar_cube);
disp(['总帧数: ', num2str(num_frames)]);

%% 初始化参数
% 雷达参数
radar_params = struct();
radar_params.range_bins = size(radar_data.radar_cube{1}, 1);
radar_params.azimuth_bins = size(radar_data.radar_cube{1}, 2);
radar_params.elevation_bins = size(radar_data.radar_cube{1}, 3);
radar_params.range_resolution = 0.1;  % m
radar_params.azimuth_resolution = 1;  % 度
radar_params.elevation_resolution = 1; % 度
radar_params.frame_time = 0.1;  % 秒

% 计算物理尺度
radar_params.range_axis = (0:radar_params.range_bins-1) * radar_params.range_resolution;
radar_params.azimuth_axis = (-radar_params.azimuth_bins/2:radar_params.azimuth_bins/2-1) * radar_params.azimuth_resolution;
radar_params.elevation_axis = (-radar_params.elevation_bins/2:radar_params.elevation_bins/2-1) * radar_params.elevation_resolution;

% 稀疏重建参数
sparse_params = struct();
sparse_params.K = 3;  % 稀疏度 (目标数量上限)
sparse_params.num_iterations = 50;  % OMP迭代次数
sparse_params.stopping_threshold = 1e-6;  % OMP停止阈值
sparse_params.use_power_method = true;  % 使用功率法加速计算

% CFAR检测参数
cfar_params = struct();
cfar_params.guard_cells = [2, 2, 2];  % 保护单元 [距离, 方位角, 俯仰角]
cfar_params.training_cells = [4, 4, 4];  % 训练单元
cfar_params.pfa = 1e-4;  % 虚警概率
cfar_params.min_distance = 1;  % 最小距离 (m)
cfar_params.min_snr = 12;  % 最小信噪比 (dB)
cfar_params.enable_neighbor_check = true;  % 启用邻域检查

% 跟踪器参数
track_params = struct();
track_params.max_targets = 2;  % 最大跟踪目标数
track_params.detection_threshold = 0.1;  % 检测阈值
track_params.confirmation_count = 3;  % 确认跟踪所需的连续检测次数
track_params.deletion_count = 5;  % 删除跟踪所需的连续未检测次数
track_params.initial_P_diag = [1.0, 1.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5]; % 初始状态协方差对角线
track_params.gating_threshold = 7.0;  % 关联门限
track_params.dt = radar_params.frame_time;  % 时间步长

% 设置滤波器噪声模型
% 基础测量噪声
R_diag = [0.25, 1.0, 1.0];  % [距离 (m), 方位角 (度), 俯仰角 (度)]
% 标准过程噪声
Q_diag = [0.1, 0.8, 0.5, 0.5, 0.8, 0.5, 0.5, 0.8, 0.5]; % 各状态的噪声: r, vr, ar, az, vaz, aaz, el, vel, ael

R_base = diag(R_diag.^2);  % 测量噪声协方差 [r, az, el]
Q_base = diag(Q_diag.^2);  % 过程噪声协方差 [r, vr, ar, az, vaz, aaz, el, vel, ael]

% 多目标跟踪算法参数
tracking_params = struct();
tracking_params.min_track_score = 0.15;  % 最小跟踪分数
tracking_params.score_decay = 0.9;  % 跟踪分数衰减系数
tracking_params.association_gate = 5.0;  % 目标关联门限 (马氏距离)
tracking_params.new_track_gate = 15.0;  % 新目标创建门限

%% 初始化结果存储
tracks = cell(num_frames, 1);  % 各帧的跟踪结果
detections = cell(num_frames, 1);  % 各帧的检测结果
track_to_truth = cell(num_frames, 1);  % 跟踪到真实值的映射
active_tracks = {};  % 当前活跃跟踪器
track_scores = [];  % 跟踪分数
track_history = {};  % 跟踪历史
kf_filters = {};  % 卡尔曼滤波器
last_update_frame = [];  % 最后更新帧
detection_history = cell(num_frames, 1);  % 检测历史

%% 跟踪循环
disp('开始目标跟踪...');

% 总帧数
fprintf('总帧数: %d\n', num_frames);

% 用于多帧集成的滑动窗口
window_size = 3;  % 集成窗口大小
radar_cube_history = cell(window_size, 1);

% 记录真实目标位置 (使用前三帧)
true_positions = zeros(3, 3);  % [距离, 方位角, 俯仰角] x 3帧
for i = 1:3
    % 提取真实目标位置
    true_positions(:, i) = radar_data.true_position{i};
end

% 计算真实目标的平均位置作为初始化
initial_position = mean(true_positions, 2);
fprintf('初始目标位置: 距离=%.2fm, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    initial_position(1), initial_position(2), initial_position(3));

% 为UKF滤波器创建初始状态向量 [r, vr, ar, az, vaz, aaz, el, vel, ael]
% 使用前三帧估计初始速度和加速度
if num_frames >= 3
    % 计算速度 (使用中心差分)
    v_r = (true_positions(1, 3) - true_positions(1, 1)) / (2 * radar_params.frame_time);
    v_az = (true_positions(2, 3) - true_positions(2, 1)) / (2 * radar_params.frame_time);
    v_el = (true_positions(3, 3) - true_positions(3, 1)) / (2 * radar_params.frame_time);
    
    % 计算加速度
    a_r = (true_positions(1, 3) - 2*true_positions(1, 2) + true_positions(1, 1)) / (radar_params.frame_time^2);
    a_az = (true_positions(2, 3) - 2*true_positions(2, 2) + true_positions(2, 1)) / (radar_params.frame_time^2);
    a_el = (true_positions(3, 3) - 2*true_positions(3, 2) + true_positions(3, 1)) / (radar_params.frame_time^2);
else
    % 如果帧数不足，假设速度和加速度为0
    v_r = 0; v_az = 0; v_el = 0;
    a_r = 0; a_az = 0; a_el = 0;
end

% 创建初始状态向量
initial_state = [initial_position(1); v_r; a_r; ...
                 initial_position(2); v_az; a_az; ...
                 initial_position(3); v_el; a_el];

% 初始化无迹卡尔曼滤波器
kf = init_ukf_filter(initial_state, R_base, Q_base, radar_params.frame_time);

%% 性能评估变量
position_errors = zeros(num_frames, 3);  % [r, az, el]
position_errors_3d = zeros(num_frames, 1);  % 3D位置误差

% 预分配RMSE计算的变量
all_range_errors = [];
all_azimuth_errors = [];
all_elevation_errors = [];
all_3d_errors = [];

% 记录最大误差
max_range_error = 0;
max_azimuth_error = 0;
max_elevation_error = 0;
max_3d_error = 0;

%% 主跟踪循环
for frame_idx = 1:num_frames
    % 加载当前雷达数据
    current_cube = radar_data.radar_cube{frame_idx};
    
    % 滑动窗口处理，保存历史数据
    if frame_idx <= window_size
        radar_cube_history{frame_idx} = current_cube;
    else
        % 滑动窗口移动
        radar_cube_history = [radar_cube_history(2:end); {current_cube}];
    end
    
    % 多帧集成处理
    if frame_idx >= window_size
        % 执行多帧非相干积累增强检测性能
        integrated_cube = zeros(size(current_cube));
        for w = 1:window_size
            integrated_cube = integrated_cube + abs(radar_cube_history{w}).^2;
        end
        integrated_cube = sqrt(integrated_cube / window_size);
    else
        % 帧数不足时，直接使用当前帧
        integrated_cube = abs(current_cube);
    end
    
    % 显示每10帧处理一次的信息
    if mod(frame_idx, 10) == 1 || frame_idx == num_frames
        fprintf('处理第 %d/%d 帧...\n', frame_idx, num_frames);
    end
    
    %% CFAR检测
    [cfar_detections, cfar_indices, detection_snr] = cfar3d_detector(integrated_cube, ...
        cfar_params.guard_cells, cfar_params.training_cells, cfar_params.pfa, ...
        cfar_params.min_distance, cfar_params.min_snr, cfar_params.enable_neighbor_check);
    
    if ~isempty(cfar_detections)
        % 转换检测结果到物理单位
        cfar_detections_physical = zeros(size(cfar_detections));
        cfar_detections_physical(:, 1) = radar_params.range_axis(cfar_detections(:, 1));
        cfar_detections_physical(:, 2) = radar_params.azimuth_axis(cfar_detections(:, 2));
        cfar_detections_physical(:, 3) = radar_params.elevation_axis(cfar_detections(:, 3));
        
        fprintf('CFAR检测: 发现 %d 个目标\n', size(cfar_detections_physical, 1));
    else
        cfar_detections_physical = [];
        fprintf('CFAR检测: 未发现目标\n');
    end
    
    %% 稀疏重建 (OMP)
    % 执行OMP稀疏重建
    [peak_indices, peak_values] = omp_3d(current_cube, sparse_params.K, sparse_params.num_iterations, ...
                                         sparse_params.stopping_threshold, sparse_params.use_power_method);
    
    % 转换OMP检测结果到物理单位
    if ~isempty(peak_indices)
        omp_detections = zeros(size(peak_indices, 1), 3);
        omp_detections(:, 1) = radar_params.range_axis(peak_indices(:, 1));
        omp_detections(:, 2) = radar_params.azimuth_axis(peak_indices(:, 2));
        omp_detections(:, 3) = radar_params.elevation_axis(peak_indices(:, 3));
        
        % 按幅值排序
        [~, sort_idx] = sort(peak_values, 'descend');
        omp_detections = omp_detections(sort_idx, :);
        peak_values = peak_values(sort_idx);
        
        fprintf('OMP检测: 发现 %d 个目标峰值\n', size(omp_detections, 1));
    else
        omp_detections = [];
        fprintf('OMP检测: 未发现目标\n');
    end
    
    %% 检测结果融合
    % 融合OMP和CFAR结果
    if ~isempty(omp_detections) && ~isempty(cfar_detections_physical)
        % 距离融合
        omp_dist = omp_detections(1, 1);
        cfar_dist = cfar_detections_physical(1, 1);
        
        % 计算与预测值的差异
        omp_dist_diff = abs(omp_dist - radar_params.range_axis(cfar_detections(1, 1)));
        cfar_dist_diff = abs(cfar_dist - radar_params.range_axis(cfar_detections(1, 1)));
        
        % 根据与预测的接近程度决定权重
        if omp_dist_diff <= cfar_dist_diff * 1.2
            % OMP距离更可靠
            dist_weight_omp = 0.6;
            dist_weight_cfar = 0.3;
            dist_weight_prior = 0.1;
        else
            % CFAR距离更可靠
            dist_weight_omp = 0.3;
            dist_weight_cfar = 0.6;
            dist_weight_prior = 0.1;
        end
        
        % 融合距离估计
        fused_range = omp_dist * dist_weight_omp + cfar_dist * dist_weight_cfar + radar_params.range_axis(cfar_detections(1, 1)) * dist_weight_prior;
        
        % 角度融合，OMP通常更准确
        omp_az = omp_detections(1, 2);
        omp_el = omp_detections(1, 3);
        
        % 检查角度偏差是否过大
        az_diff = abs(wrapTo180(omp_az - radar_params.azimuth_axis(cfar_detections(1, 2))));
        el_diff = abs(wrapTo180(omp_el - radar_params.elevation_axis(cfar_detections(1, 3))));
        
        % 角度偏差阈值
        angle_diff_threshold = 5.0;
        
        if az_diff > angle_diff_threshold
            % 方位角偏差过大，增加预测值权重
            az_weight_omp = 0.5;
            az_weight_prior = 0.5;
            
            % 自适应融合
            fused_az = wrapTo180(omp_az * az_weight_omp + radar_params.azimuth_axis(cfar_detections(1, 2)) * az_weight_prior);
        else
            % 方位角偏差合理，信任OMP
            fused_az = omp_az;
        end
        
        if el_diff > angle_diff_threshold
            % 俯仰角偏差过大，增加预测值权重
            el_weight_omp = 0.5;
            el_weight_prior = 0.5;
            
            % 自适应融合
            fused_el = wrapTo180(omp_el * el_weight_omp + radar_params.elevation_axis(cfar_detections(1, 3)) * el_weight_prior);
        else
            % 俯仰角偏差合理，信任OMP
            fused_el = omp_el;
        end
        
        % 设置测量值
        z = [fused_range; fused_az; fused_el];
        
        % 输出关键信息
        fprintf('融合测量: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
                fused_range, fused_az, fused_el);
    elseif ~isempty(omp_detections)
        % 只有OMP检测
        z = [omp_detections(1, 1); omp_detections(1, 2); omp_detections(1, 3)];
        fprintf('警告: 帧 %d 没有CFAR检测，使用OMP检测\n', frame_idx);
    elseif ~isempty(cfar_detections_physical)
        % 只有CFAR检测
        z = [cfar_detections_physical(1, 1); radar_params.azimuth_axis(cfar_detections(1, 2)); radar_params.elevation_axis(cfar_detections(1, 3))];
        fprintf('警告: 帧 %d 没有OMP检测，使用CFAR检测\n', frame_idx);
    else
        % 无检测
        z = [NaN; NaN; NaN];
        fprintf('警告: 帧 %d 没有有效检测！使用仅预测模式\n', frame_idx);
    end
    
    % 保存检测结果
    detection_info = struct();
    detection_info.positions = z;
    detection_info.weights = detection_snr;
    detection_info.source = 3;  % 表示融合检测
    detection_history{frame_idx} = detection_info;
    
    %% 目标跟踪
    % 获取真实位置
    true_position = radar_data.true_position{frame_idx};
    
    % 更新卡尔曼滤波器
    [kf, innovation] = update_ukf_filter(kf, z, radar_params.frame_time, frame_idx);
    
    % 提取当前状态估计
    estimated_r = kf.x(1);     % 距离估计
    estimated_vr = kf.x(2);    % 径向速度估计
    estimated_ar = kf.x(3);    % 径向加速度估计
    estimated_az = kf.x(4);    % 方位角估计
    estimated_vaz = kf.x(5);   % 方位角速度估计
    estimated_aaz = kf.x(6);   % 方位角加速度估计
    estimated_el = kf.x(7);    % 俯仰角估计
    estimated_vel = kf.x(8);   % 俯仰角速度估计
    estimated_ael = kf.x(9);   % 俯仰角加速度估计
    
    % 计算跟踪误差
    r_error = abs(estimated_r - true_position(1));
    az_error = abs(wrapTo180(estimated_az - true_position(2)));
    el_error = abs(wrapTo180(estimated_el - true_position(3)));
    
    % 计算3D位置误差 (考虑球坐标到笛卡尔坐标的转换)
    true_x = true_position(1) * cosd(true_position(3)) * cosd(true_position(2));
    true_y = true_position(1) * cosd(true_position(3)) * sind(true_position(2));
    true_z = true_position(1) * sind(true_position(3));
    
    est_x = estimated_r * cosd(estimated_el) * cosd(estimated_az);
    est_y = estimated_r * cosd(estimated_el) * sind(estimated_az);
    est_z = estimated_r * sind(estimated_el);
    
    error_3d = sqrt((est_x - true_x)^2 + (est_y - true_y)^2 + (est_z - true_z)^2);
    
    % 保存误差
    position_errors(frame_idx, :) = [r_error, az_error, el_error];
    position_errors_3d(frame_idx) = error_3d;
    
    % 更新最大误差记录
    max_range_error = max(max_range_error, r_error);
    max_azimuth_error = max(max_azimuth_error, az_error);
    max_elevation_error = max(max_elevation_error, el_error);
    max_3d_error = max(max_3d_error, error_3d);
    
    % 累积误差用于RMSE计算
    all_range_errors = [all_range_errors; r_error];
    all_azimuth_errors = [all_azimuth_errors; az_error];
    all_elevation_errors = [all_elevation_errors; el_error];
    all_3d_errors = [all_3d_errors; error_3d];
    
    % 每10帧或最后一帧输出当前跟踪性能
    if mod(frame_idx, 10) == 0 || frame_idx == num_frames
        fprintf('帧 %d: 跟踪误差 - 距离: %.3f m, 方位角: %.3f°, 俯仰角: %.3f°, 3D位置: %.3f m\n', ...
            frame_idx, r_error, az_error, el_error, error_3d);
    end
    
    % 保存跟踪结果
    track_result = struct();
    track_result.estimated_state = kf.x;
    track_result.position = [estimated_r, estimated_az, estimated_el];
    track_result.velocity = [estimated_vr, estimated_vaz, estimated_vel];
    track_result.acceleration = [estimated_ar, estimated_aaz, estimated_ael];
    track_result.covariance = kf.P;
    track_result.innovation = innovation;
    
    tracks{frame_idx} = track_result;
    
    % 记录跟踪到真实值映射
    track_to_truth{frame_idx} = true_position;
end

%% 计算跟踪性能指标
% 计算RMSE
rmse_range = sqrt(mean(all_range_errors.^2));
rmse_azimuth = sqrt(mean(all_azimuth_errors.^2));
rmse_elevation = sqrt(mean(all_elevation_errors.^2));
rmse_3d = sqrt(mean(all_3d_errors.^2));

% 提前3帧的性能指标（跳过前3帧，因为它们用于初始化）
if num_frames > 3
    early_range_errors = all_range_errors(4:min(10, num_frames));
    early_azimuth_errors = all_azimuth_errors(4:min(10, num_frames));
    early_elevation_errors = all_elevation_errors(4:min(10, num_frames));
    early_3d_errors = all_3d_errors(4:min(10, num_frames));
    
    rmse_range_early = sqrt(mean(early_range_errors.^2));
    rmse_azimuth_early = sqrt(mean(early_azimuth_errors.^2));
    rmse_elevation_early = sqrt(mean(early_elevation_errors.^2));
    rmse_3d_early = sqrt(mean(early_3d_errors.^2));
else
    rmse_range_early = NaN;
    rmse_azimuth_early = NaN;
    rmse_elevation_early = NaN;
    rmse_3d_early = NaN;
end

% 后半部分的性能指标
if num_frames > 6
    late_idx = floor(num_frames/2):num_frames;
    late_range_errors = all_range_errors(late_idx);
    late_azimuth_errors = all_azimuth_errors(late_idx);
    late_elevation_errors = all_elevation_errors(late_idx);
    late_3d_errors = all_3d_errors(late_idx);
    
    rmse_range_late = sqrt(mean(late_range_errors.^2));
    rmse_azimuth_late = sqrt(mean(late_azimuth_errors.^2));
    rmse_elevation_late = sqrt(mean(late_elevation_errors.^2));
    rmse_3d_late = sqrt(mean(late_3d_errors.^2));
else
    rmse_range_late = NaN;
    rmse_azimuth_late = NaN;
    rmse_elevation_late = NaN;
    rmse_3d_late = NaN;
end

%% 输出性能报告
fprintf('\n=================== 跟踪性能报告 ===================\n');
fprintf('总帧数: %d\n', num_frames);
fprintf('\n整体性能指标:\n');
fprintf('距离估计RMSE: %.3f m (最大误差: %.3f m)\n', rmse_range, max_range_error);
fprintf('方位角估计RMSE: %.3f° (最大误差: %.3f°)\n', rmse_azimuth, max_azimuth_error);
fprintf('俯仰角估计RMSE: %.3f° (最大误差: %.3f°)\n', rmse_elevation, max_elevation_error);
fprintf('3D位置估计RMSE: %.3f m (最大误差: %.3f m)\n', rmse_3d, max_3d_error);

fprintf('\n早期跟踪性能 (帧4-%d):\n', min(10, num_frames));
fprintf('距离估计RMSE: %.3f m\n', rmse_range_early);
fprintf('方位角估计RMSE: %.3f°\n', rmse_azimuth_early);
fprintf('俯仰角估计RMSE: %.3f°\n', rmse_elevation_early);
fprintf('3D位置估计RMSE: %.3f m\n', rmse_3d_early);

fprintf('\n后期跟踪性能 (帧%d-%d):\n', floor(num_frames/2), num_frames);
fprintf('距离估计RMSE: %.3f m\n', rmse_range_late);
fprintf('方位角估计RMSE: %.3f°\n', rmse_azimuth_late);
fprintf('俯仰角估计RMSE: %.3f°\n', rmse_elevation_late);
fprintf('3D位置估计RMSE: %.3f m\n', rmse_3d_late);

fprintf('\n检测性能:\n');
% 计算平均每帧检测数
total_detections = 0;
cfar_only = 0;
omp_only = 0;
both_methods = 0;

for i = 1:num_frames
    if ~isempty(detection_history{i})
        total_detections = total_detections + length(detection_history{i}.source);
        cfar_only = cfar_only + sum(detection_history{i}.source == 1);
        omp_only = omp_only + sum(detection_history{i}.source == 2);
        both_methods = both_methods + sum(detection_history{i}.source == 3);
    end
end

fprintf('平均每帧检测数: %.2f\n', total_detections/num_frames);
fprintf('CFAR专有检测: %.2f%%\n', 100 * cfar_only / max(1, total_detections));
fprintf('OMP专有检测: %.2f%%\n', 100 * omp_only / max(1, total_detections));
fprintf('两种方法共同检测: %.2f%%\n', 100 * both_methods / max(1, total_detections));

fprintf('\n=====================================================\n');

%% 绘制跟踪结果
% 提取轨迹数据
track_positions = zeros(num_frames, 3);
true_track = zeros(num_frames, 3);

for i = 1:num_frames
    if ~isempty(tracks{i})
        track_positions(i, :) = tracks{i}.position;
    end
    true_track(i, :) = track_to_truth{i};
end

% 创建图形
figure('Name', '3D目标跟踪结果', 'Position', [100, 100, 1000, 800]);

% 距离误差随时间变化
subplot(2, 2, 1);
plot(1:num_frames, position_errors(:, 1), 'LineWidth', 2);
hold on;
plot([1, num_frames], [rmse_range, rmse_range], 'r--', 'LineWidth', 1.5);
xlabel('帧索引');
ylabel('距离误差 (m)');
title(sprintf('距离估计误差 (RMSE: %.3f m)', rmse_range));
grid on;

% 方位角误差随时间变化
subplot(2, 2, 2);
plot(1:num_frames, position_errors(:, 2), 'LineWidth', 2);
hold on;
plot([1, num_frames], [rmse_azimuth, rmse_azimuth], 'r--', 'LineWidth', 1.5);
xlabel('帧索引');
ylabel('方位角误差 (度)');
title(sprintf('方位角估计误差 (RMSE: %.3f°)', rmse_azimuth));
grid on;

% 俯仰角误差随时间变化
subplot(2, 2, 3);
plot(1:num_frames, position_errors(:, 3), 'LineWidth', 2);
hold on;
plot([1, num_frames], [rmse_elevation, rmse_elevation], 'r--', 'LineWidth', 1.5);
xlabel('帧索引');
ylabel('俯仰角误差 (度)');
title(sprintf('俯仰角估计误差 (RMSE: %.3f°)', rmse_elevation));
grid on;

% 3D位置误差随时间变化
subplot(2, 2, 4);
plot(1:num_frames, position_errors_3d, 'LineWidth', 2);
hold on;
plot([1, num_frames], [rmse_3d, rmse_3d], 'r--', 'LineWidth', 1.5);
xlabel('帧索引');
ylabel('3D位置误差 (m)');
title(sprintf('3D位置估计误差 (RMSE: %.3f m)', rmse_3d));
grid on;

% 创建3D轨迹图
figure('Name', '目标轨迹', 'Position', [100, 100, 800, 600]);

% 转换到笛卡尔坐标
true_x = true_track(:, 1) .* cosd(true_track(:, 3)) .* cosd(true_track(:, 2));
true_y = true_track(:, 1) .* cosd(true_track(:, 3)) .* sind(true_track(:, 2));
true_z = true_track(:, 1) .* sind(true_track(:, 3));

est_x = track_positions(:, 1) .* cosd(track_positions(:, 3)) .* cosd(track_positions(:, 2));
est_y = track_positions(:, 1) .* cosd(track_positions(:, 3)) .* sind(track_positions(:, 2));
est_z = track_positions(:, 1) .* sind(track_positions(:, 3));

% 绘制3D轨迹
plot3(true_x, true_y, true_z, 'b-', 'LineWidth', 2);
hold on;
plot3(est_x, est_y, est_z, 'r--', 'LineWidth', 2);

% 添加起始点和结束点标记
plot3(true_x(1), true_y(1), true_z(1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(true_x(end), true_y(end), true_z(end), 'b*', 'MarkerSize', 10);
plot3(est_x(1), est_y(1), est_z(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot3(est_x(end), est_y(end), est_z(end), 'r*', 'MarkerSize', 10);

% 添加轴标签和标题
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('目标轨迹 (蓝色=真实, 红色=估计)');
grid on;
legend('真实轨迹', '估计轨迹', '真实起点', '真实终点', '估计起点', '估计终点', 'Location', 'best');
view(30, 30);

% 保存结果
save('tracking_results.mat', 'tracks', 'track_to_truth', 'position_errors', ...
     'position_errors_3d', 'rmse_range', 'rmse_azimuth', 'rmse_elevation', 'rmse_3d');

disp('跟踪演示完成，结果已保存。'); 