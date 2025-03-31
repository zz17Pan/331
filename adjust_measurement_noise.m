function R = adjust_measurement_noise(measurement, z_pred, base_r, frame_count)
%ADJUST_MEASUREMENT_NOISE 调整测量噪声
%   measurement: 当前测量值
%   z_pred: 预测测量值
%   base_r: 基础测量噪声协方差
%   frame_count: 当前帧计数
%   R: 调整后的测量噪声协方差

% 提取测量值
measured_dist = measurement(1);
measured_az = measurement(2);
measured_el = measurement(3);

% 提取预测测量值
predicted_dist = z_pred(1);
predicted_az = z_pred(2);
predicted_el = z_pred(3);

% 计算测量残差
dist_residual = abs(measured_dist - predicted_dist);
az_residual = abs(wrapTo180(measured_az - predicted_az));
el_residual = abs(wrapTo180(measured_el - predicted_el));

% 基础噪声
r_range = base_r(1,1);
r_azimuth = base_r(2,2);
r_elevation = base_r(3,3);

% 根据残差调整噪声
% 距离噪声调整
if dist_residual > 5.0
    r_range = r_range * (1.0 + dist_residual/10.0);
elseif dist_residual < 1.0
    r_range = r_range * 0.9;  % 减小噪声，更信任测量
end

% 方位角噪声调整
if az_residual > 3.0
    r_azimuth = r_azimuth * (1.0 + az_residual/5.0);
elseif az_residual < 0.5
    r_azimuth = r_azimuth * 0.9;
end

% 俯仰角噪声调整
if el_residual > 3.0
    r_elevation = r_elevation * (1.0 + el_residual/5.0);
elseif el_residual < 0.5
    r_elevation = r_elevation * 0.9;
end

% 距离对噪声的影响
distance_factor = min(2.0, max(0.8, 1.0 + (measured_dist - 50)/100));
r_range = r_range * distance_factor;

% 在初始帧更信任测量
if frame_count <= 3
    confidence_boost = 0.8;  % 增加对测量的信任
    r_range = r_range * confidence_boost;
    r_azimuth = r_azimuth * confidence_boost;
    r_elevation = r_elevation * confidence_boost;
end

% 构建最终噪声协方差矩阵
R = diag([r_range, r_azimuth, r_elevation]);

return;
end 