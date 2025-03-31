function [x, y, z] = sph2cart_sensing(r, az, el)
%SPH2CART_SENSING 将球坐标系转换为笛卡尔坐标系，适用于感知系统
%   [x, y, z] = sph2cart_sensing(r, az, el) 将球坐标系的坐标转换为笛卡尔坐标系
%   r: 距离(m)
%   az: 方位角(度)，水平面内从x轴正方向逆时针旋转的角度
%   el: 俯仰角(度)，与水平面的夹角，向上为正
%   x, y, z: 笛卡尔坐标

% 转换为弧度
az_rad = az * pi / 180;
el_rad = el * pi / 180;

% 球坐标系到笛卡尔坐标系的转换
% 注意：感知系统中，方位角是在水平面(xy平面)内，从x轴开始逆时针旋转
% 俯仰角是与xy平面的夹角，向上为正
x = r * cos(el_rad) * cos(az_rad);
y = r * cos(el_rad) * sin(az_rad);
z = r * sin(el_rad);
end 