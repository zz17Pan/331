function [r, az, el] = cart2sph_sensing(x, y, z)
%CART2SPH_SENSING 将笛卡尔坐标系转换为球坐标系，适用于感知系统
%   [r, az, el] = cart2sph_sensing(x, y, z) 将笛卡尔坐标系转换为球坐标系
%   x, y, z: 笛卡尔坐标
%   r: 距离(m)
%   az: 方位角(度)，水平面内从x轴正方向逆时针旋转的角度
%   el: 俯仰角(度)，与水平面的夹角，向上为正

% 计算距离
r = sqrt(x^2 + y^2 + z^2);

% 计算方位角(弧度)
az_rad = atan2(y, x);

% 计算俯仰角(弧度)
if r < 1e-10  % 防止除以接近零的值
    el_rad = 0;
else
    el_rad = asin(z / r);
end

% 转换为度
az = az_rad * 180 / pi;
el = el_rad * 180 / pi;

% 确保方位角在[-180, 180]度范围内
az = wrapTo180(az);
end 