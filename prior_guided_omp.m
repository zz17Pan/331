function [estimated_range, estimated_azimuth, estimated_elevation] = prior_guided_omp(rx_signal, tx_array, rx_array, prior_info, prior_cov, params)
% PRIOR_GUIDED_OMP 高精度稀疏重建算法
% 目标：距离误差<1m，角度误差<1度

% 信号预处理
y = rx_signal(:);
if length(y) > 4000  % 降低到4000以提高计算速度
    down_factor = ceil(length(y) / 4000);
    y = y(1:down_factor:end);
end
params.current_signal_length = length(y);

% 提取CFAR检测结果
has_detection = isfield(params, 'measurements') && isfield(params.measurements, 'has_cfar') && params.measurements.has_cfar;
if has_detection
    cfar_range = params.measurements.cfar_range;
    cfar_velocity = params.measurements.cfar_velocity;
    % 检查CFAR结果是否具有可靠性信息
    if isfield(params.measurements, 'reliability')
        cfar_reliability = params.measurements.reliability;
    else
        cfar_reliability = 0.7; % 默认可靠性
    end
    fprintf('OMP算法使用CFAR检测: 距离=%.2f m, 速度=%.2f m/s, 可靠性=%.2f\n', ...
            cfar_range, cfar_velocity, cfar_reliability);
end

% 提取MUSIC角度估计结果
has_music = isfield(params, 'music_results') && isfield(params.music_results, 'reliability');
if has_music
    music_azimuth = params.music_results.azimuth;
    music_elevation = params.music_results.elevation;
    music_reliability = params.music_results.reliability;
    fprintf('OMP算法使用MUSIC角度估计: 方位角=%.2f°, 俯仰角=%.2f°, 可靠性=%.2f\n', ...
            music_azimuth, music_elevation, music_reliability);
end

% 获取先验标准差并计算基础搜索范围
prior_std = sqrt(diag(prior_cov));
base_range_dev = max(0.3, prior_std(1));    % 降低基础距离搜索范围
base_angle_dev = max(0.8, prior_std(2));    % 降低基础角度搜索范围

% 动态调整搜索范围
velocity_norm = norm(params.target_state.velocity);
velocity_factor = min(1.2, velocity_norm / 5.0);  % 考虑速度

% 计算自适应搜索范围
range_dev = base_range_dev * (1 + 0.2 * velocity_factor);  % 减小系数以收窄范围
angle_dev = base_angle_dev * (1 + 0.2 * velocity_factor);

% 设置网格密度
N_range = 30;  % 减少网格数以提高速度
N_az = 30;
N_el = 30;

% 设置搜索范围中心和宽度
% 根据CFAR结果和先验信息确定距离搜索中心
if has_detection
    % 根据CFAR可靠性加权计算搜索中心
    range_center = prior_info.range * (1 - cfar_reliability) + cfar_range * cfar_reliability;
    
    % 根据CFAR和先验的差异调整搜索范围
    range_dev = max(range_dev, abs(prior_info.range - cfar_range) * 0.6);
else
    range_center = prior_info.range;
end

% 如果有MUSIC结果，使用它来调整角度搜索中心
if has_music && music_reliability > 0.5
    % 按可靠性加权MUSIC结果
    music_weight = music_reliability * 0.8;  % 适当权重
    prior_weight = 1 - music_weight;
    
    az_center = prior_info.azimuth * prior_weight + music_azimuth * music_weight;
    el_center = prior_info.elevation * prior_weight + music_elevation * music_weight;
    
    % 根据MUSIC可靠性调整搜索范围
    angle_dev = angle_dev * (2.0 - music_reliability);  % 根据可靠性减小搜索范围
else
    az_center = prior_info.azimuth;
    el_center = prior_info.elevation;
end

% 定义搜索范围
range_min = max(0.1, range_center - range_dev);
range_max = range_center + range_dev;
az_min = az_center - angle_dev;
az_max = az_center + angle_dev;
el_min = el_center - angle_dev;
el_max = el_center + angle_dev;

% 确保角度在有效范围内
az_min = wrapTo180(az_min);
az_max = wrapTo180(az_max);
el_min = max(-90, min(90, el_min));
el_max = max(-90, min(90, el_max));

% 生成自适应网格
range_grid = generate_adaptive_grid(range_min, range_max, N_range, 'range');
az_grid = generate_adaptive_grid(az_min, az_max, N_az, 'angle');
el_grid = generate_adaptive_grid(el_min, el_max, N_el, 'angle');

% 初始化OMP算法变量
dict_size = N_range * N_az * N_el;
max_iter = min(params.omp.max_iter, 6);  % 减少最大迭代次数提高速度
residual = y;
selected_atoms = zeros(length(y), max_iter);
selected_indices = zeros(1, max_iter);
actual_selected_count = 0;
power_values = zeros(N_range, N_az, N_el);
iteration_powers = [];

% 优化块大小 - 分块处理减少内存消耗
block_size = min(40, dict_size);  % 减小块大小以提高内存效率
num_blocks = ceil(dict_size / block_size);

% OMP主迭代
for iter = 1:max_iter
    max_corr = 0;
    max_idx = 0;
    max_atom = [];
    
    for block = 1:num_blocks
        start_idx = (block-1) * block_size + 1;
        end_idx = min(block * block_size, dict_size);
        current_block_size = end_idx - start_idx + 1;
        
        block_atoms = zeros(length(y), current_block_size);
        block_indices = zeros(1, current_block_size);
        index = 1;
        
        for i_range = 1:N_range
            r = range_grid(i_range);
            for i_az = 1:N_az
                az = az_grid(i_az);
                for i_el = 1:N_el
                    el = el_grid(i_el);
                    
                    global_idx = i_range + (i_az-1)*N_range + (i_el-1)*N_range*N_az;
                    
                    if global_idx >= start_idx && global_idx <= end_idx
                        % 生成导向矢量
                        [a_tx, a_rx] = compute_steering_vector(tx_array, rx_array, r, az, el, params);
                        atom = generate_atom(a_tx, a_rx, r, params);
                        
                        % 确保原子长度匹配
                        if length(atom) ~= length(y)
                            if length(atom) > length(y)
                                atom = atom(1:length(y));
                            else
                                atom = [atom; zeros(length(y)-length(atom), 1)];
                            end
                        end
                        
                        block_atoms(:, index) = atom;
                        block_indices(index) = global_idx;
                        index = index + 1;
                    end
                end
            end
        end
        
        actual_block_size = index - 1;
        if actual_block_size > 0
            block_atoms = block_atoms(:, 1:actual_block_size);
            block_indices = block_indices(1:actual_block_size);
            
            % 计算相关性
            batch_correlation = abs(block_atoms' * residual);
            [batch_max_corr, batch_max_idx] = max(batch_correlation);
            
            if batch_max_corr > max_corr
                max_corr = batch_max_corr;
                max_idx = block_indices(batch_max_idx);
                max_atom = block_atoms(:, batch_max_idx);
            end
        end
    end
    
    % 改进的迭代终止条件
    if iter > 1
        power_gain_ratio = max_corr / iteration_powers(end);
        if (power_gain_ratio < 1.03 && iter > 2) || ...  % 降低阈值加速收敛
           (norm(residual) < params.omp.residual_tol * norm(y))
            fprintf('OMP迭代终止: 第%d次迭代，功率增益比=%.4f\n', iter, power_gain_ratio);
            break;
        end
    end
    
    iteration_powers(end+1) = max_corr;
    actual_selected_count = actual_selected_count + 1;
    
    % 更新选定的原子和系数
    selected_atoms(:, actual_selected_count) = max_atom;
    selected_indices(actual_selected_count) = max_idx;
    coeffs = selected_atoms(:, 1:actual_selected_count) \ y;
    residual = y - selected_atoms(:, 1:actual_selected_count) * coeffs;
    
    % 更新能量值
    [i_range, i_az, i_el] = ind2sub([N_range, N_az, N_el], max_idx);
    power_values(i_range, i_az, i_el) = abs(coeffs(end))^2;
    
    % 输出迭代信息
    fprintf('OMP迭代 %d: 选择原子 [距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°], 功率=%.2e\n', ...
        iter, range_grid(i_range), az_grid(i_az), el_grid(i_el), power_values(i_range, i_az, i_el));
end

% 高精度参数估计 - 通过亚栅格插值
[~, max_idx] = max(power_values(:));
[i_range, i_az, i_el] = ind2sub(size(power_values), max_idx);
fprintf('OMP原始最大功率点: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    range_grid(i_range), az_grid(i_az), el_grid(i_el));

% 二次插值优化参数
[estimated_range, range_power] = subgrid_refinement(range_grid, power_values, i_range, i_az, i_el, 1);
[estimated_azimuth, az_power] = subgrid_refinement(az_grid, power_values, i_range, i_az, i_el, 2);
[estimated_elevation, el_power] = subgrid_refinement(el_grid, power_values, i_range, i_az, i_el, 3);

fprintf('OMP插值优化后: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
    estimated_range, estimated_azimuth, estimated_elevation);

% 计算测量质量
quality = compute_measurement_quality(residual, y, power_values, range_power, az_power, el_power);
fprintf('OMP重建质量: %.2f\n', quality);

% 结果融合 - 科学地融合各个估计器的结果
% 考虑每个测量结果的可靠性
% 为每个测量设置权重，考虑质量和距离
if has_detection && has_music && music_reliability > 0.5
    % 构建科学的融合模型
    
    % 考虑OMP质量和距离的权重
    omp_quality = quality;
    range_factor = 1.0 / (1.0 + 0.05 * estimated_range); % 距离越远，OMP权重越小
    omp_weight_range = min(0.6, omp_quality * range_factor);
    omp_weight_angle = min(0.5, omp_quality);
    
    % CFAR权重依赖于测量可靠性
    cfar_weight = min(0.4, cfar_reliability);
    
    % MUSIC角度权重
    music_weight = min(0.6, music_reliability);
    
    % 先验权重 (剩余权重)
    prior_weight_range = max(0.1, 1 - omp_weight_range - cfar_weight);
    prior_weight_angle = max(0.1, 1 - omp_weight_angle - music_weight);
    
    % 融合距离估计
    estimated_range = omp_weight_range * estimated_range + ...
                     cfar_weight * cfar_range + ...
                     prior_weight_range * prior_info.range;
                     
    % 融合角度估计
    estimated_azimuth = omp_weight_angle * estimated_azimuth + ...
                       music_weight * music_azimuth + ...
                       prior_weight_angle * prior_info.azimuth;
                       
    estimated_elevation = omp_weight_angle * estimated_elevation + ...
                         music_weight * music_elevation + ...
                         prior_weight_angle * prior_info.elevation;
    
    fprintf('结果融合权重-距离: OMP=%.2f, CFAR=%.2f, Prior=%.2f\n', ...
            omp_weight_range, cfar_weight, prior_weight_range);
    fprintf('结果融合权重-角度: OMP=%.2f, MUSIC=%.2f, Prior=%.2f\n', ...
            omp_weight_angle, music_weight, prior_weight_angle);
            
elseif has_detection
    % 只有CFAR结果时
    omp_quality = quality;
    range_factor = 1.0 / (1.0 + 0.05 * estimated_range);
    omp_weight_range = min(0.5, omp_quality * range_factor);
    omp_weight_angle = min(0.6, omp_quality);
    
    cfar_weight = min(0.5, cfar_reliability);
    
    prior_weight_range = max(0.1, 1 - omp_weight_range - cfar_weight);
    prior_weight_angle = max(0.2, 1 - omp_weight_angle);
    
    % 融合估计
    estimated_range = omp_weight_range * estimated_range + ...
                     cfar_weight * cfar_range + ...
                     prior_weight_range * prior_info.range;
                     
    estimated_azimuth = omp_weight_angle * estimated_azimuth + ...
                       prior_weight_angle * prior_info.azimuth;
                       
    estimated_elevation = omp_weight_angle * estimated_elevation + ...
                         prior_weight_angle * prior_info.elevation;
                         
    fprintf('结果融合权重-距离: OMP=%.2f, CFAR=%.2f, Prior=%.2f\n', ...
            omp_weight_range, cfar_weight, prior_weight_range);
    fprintf('结果融合权重-角度: OMP=%.2f, Prior=%.2f\n', ...
            omp_weight_angle, prior_weight_angle);
            
elseif has_music && music_reliability > 0.5
    % 只有MUSIC结果时
    omp_quality = quality;
    
    omp_weight_range = min(0.6, omp_quality);
    omp_weight_angle = min(0.4, omp_quality);
    
    music_weight = min(0.6, music_reliability);
    
    prior_weight_range = max(0.2, 1 - omp_weight_range);
    prior_weight_angle = max(0.1, 1 - omp_weight_angle - music_weight);
    
    % 融合估计
    estimated_range = omp_weight_range * estimated_range + ...
                     prior_weight_range * prior_info.range;
                     
    estimated_azimuth = omp_weight_angle * estimated_azimuth + ...
                       music_weight * music_azimuth + ...
                       prior_weight_angle * prior_info.azimuth;
                       
    estimated_elevation = omp_weight_angle * estimated_elevation + ...
                         music_weight * music_elevation + ...
                         prior_weight_angle * prior_info.elevation;
                         
    fprintf('结果融合权重-距离: OMP=%.2f, Prior=%.2f\n', ...
            omp_weight_range, prior_weight_range);
    fprintf('结果融合权重-角度: OMP=%.2f, MUSIC=%.2f, Prior=%.2f\n', ...
            omp_weight_angle, music_weight, prior_weight_angle);
else
    % 无CFAR和MUSIC结果时
    omp_quality = quality;
    
    omp_weight = min(0.6, omp_quality);
    prior_weight = max(0.4, 1 - omp_weight);
    
    % 融合估计
    estimated_range = omp_weight * estimated_range + prior_weight * prior_info.range;
    estimated_azimuth = omp_weight * estimated_azimuth + prior_weight * prior_info.azimuth;
    estimated_elevation = omp_weight * estimated_elevation + prior_weight * prior_info.elevation;
    
    fprintf('结果融合权重: OMP=%.2f, Prior=%.2f\n', omp_weight, prior_weight);
end

% 应用动态变化率限制
persistent prev_az prev_el prev_range prev_frame
if isempty(prev_az)
    prev_az = prior_info.azimuth;
    prev_el = prior_info.elevation;
    prev_range = prior_info.range;
    prev_frame = 0;
end

% 计算时间增量
current_frame = params.frame_idx;
dt = 0.1;  % 默认帧间隔
if isfield(params, 'frame_interval')
    dt = params.frame_interval;
end

% 计算动态最大变化率
max_range_rate = velocity_norm * 1.5;  % 允许速度裕度
max_angle_rate = rad2deg(atan2(velocity_norm, estimated_range)) * 1.8;  % 考虑角速度

% 应用变化率限制
range_rate = (estimated_range - prev_range) / dt;
az_rate = wrapTo180(estimated_azimuth - prev_az) / dt;
el_rate = (estimated_elevation - prev_el) / dt;

if abs(range_rate) > max_range_rate
    estimated_range = prev_range + sign(range_rate) * max_range_rate * dt;
    fprintf('限制距离变化率: %.2f -> %.2f m/s\n', range_rate, sign(range_rate) * max_range_rate);
end
if abs(az_rate) > max_angle_rate
    estimated_azimuth = prev_az + sign(az_rate) * max_angle_rate * dt;
    fprintf('限制方位角变化率: %.2f -> %.2f °/s\n', az_rate, sign(az_rate) * max_angle_rate);
end
if abs(el_rate) > max_angle_rate
    estimated_elevation = prev_el + sign(el_rate) * max_angle_rate * dt;
    fprintf('限制俯仰角变化率: %.2f -> %.2f °/s\n', el_rate, sign(el_rate) * max_angle_rate);
end

% 更新历史值
prev_az = estimated_azimuth;
prev_el = estimated_elevation;
prev_range = estimated_range;
prev_frame = current_frame;

% 确保结果在物理有效范围内
estimated_range = max(0.1, estimated_range);
estimated_azimuth = wrapTo180(estimated_azimuth);
estimated_elevation = max(-90, min(90, estimated_elevation));

% 最终输出
fprintf('高精度估计结果: 距离=%.2f m, 方位角=%.2f°, 俯仰角=%.2f°\n', ...
        estimated_range, estimated_azimuth, estimated_elevation);
end

function grid = generate_adaptive_grid(min_val, max_val, N, type)
    if min_val == max_val
        grid = min_val * ones(1, N);
        return;
    end
    
    if strcmp(type, 'range')
        % 距离网格：使用非线性分布，中心更密集
        ratio = max_val / min_val;
        beta = 1.0;  % 控制网格密度分布
        t = linspace(0, 1, N).^beta;
        grid = min_val * (1 + t * (ratio - 1));
    else
        % 角度网格：中心更密集
        center = (min_val + max_val) / 2;
        half_width = (max_val - min_val) / 2;
        t = linspace(-1, 1, N);
        grid = center + half_width * sinh(2.0*t)/sinh(2.0);  % 使用双曲正弦函数增加中心密度
    end
    
    % 确保网格是行向量
    grid = reshape(grid, 1, []);
end

function [refined_value, refined_power] = subgrid_refinement(grid, power_values, i_range, i_az, i_el, dim)
    % 二次插值优化参数估计
    if dim == 1 && i_range > 1 && i_range < size(power_values,1)
        values = grid(i_range-1:i_range+1);
        powers = squeeze(power_values(i_range-1:i_range+1, i_az, i_el));
    elseif dim == 2 && i_az > 1 && i_az < size(power_values,2)
        values = grid(i_az-1:i_az+1);
        powers = squeeze(power_values(i_range, i_az-1:i_az+1, i_el));
    elseif dim == 3 && i_el > 1 && i_el < size(power_values,3)
        values = grid(i_el-1:i_el+1);
        powers = squeeze(power_values(i_range, i_az, i_el-1:i_el+1));
    else
        if dim == 1
            refined_value = grid(i_range);
        elseif dim == 2
            refined_value = grid(i_az);
        else
            refined_value = grid(i_el);
        end
        refined_power = power_values(i_range, i_az, i_el);
        return;
    end
    
    % 执行二次插值
    p = polyfit(values, powers, 2);
    
    if p(1) >= 0  % 确保二次曲线为凹函数
        refined_value = values(2);  % 使用中心点
        refined_power = powers(2);
    else
        % 计算最大值点
        refined_value = -p(2)/(2*p(1));
        
        % 确保结果在合理范围内
        if refined_value < values(1) || refined_value > values(3)
            refined_value = values(2);  % 回退到原始栅格点
        end
        
        % 计算对应功率
        refined_power = polyval(p, refined_value);
    end
end

function quality = compute_measurement_quality(residual, original_signal, power_values, range_power, az_power, el_power)
    % 计算残差质量
    relative_residual = norm(residual) / norm(original_signal);
    residual_quality = max(0, 1 - relative_residual);
    
    % 计算能量分布质量
    total_power = sum(power_values(:));
    if total_power == 0
        distribution_quality = 0;
    else
        peak_power = max([range_power, az_power, el_power]);
        power_ratio = peak_power / total_power;
        distribution_quality = min(1, power_ratio * 2.5);  % 强调峰值特性
    end
    
    % 计算空间一致性
    peak_power = max([range_power, az_power, el_power]);
    if peak_power == 0
        consistency = 0;
    else
        range_power_norm = range_power / peak_power;
        az_power_norm = az_power / peak_power;
        el_power_norm = el_power / peak_power;
        consistency = mean([range_power_norm, az_power_norm, el_power_norm]);
    end
    
    % 综合质量评估，采用加权平均
    quality = 0.5 * residual_quality + ...
             0.3 * distribution_quality + ...
             0.2 * consistency;
             
    % 限制输出范围
    quality = min(0.95, max(0.1, quality));
end