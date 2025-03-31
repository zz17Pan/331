function x_next = nonlinear_state_transition(x, F, frame_count)
%NONLINEAR_STATE_TRANSITION 非线性状态转移函数
%   x: 当前状态向量
%   F: 状态转移矩阵
%   frame_count: 当前帧计数
%   x_next: 下一个状态向量

% 提取状态变量用于分析
r = x(1);         % 距离
vr = x(2);        % 径向速度
ar = x(3);        % 径向加速度
az = x(4);        % 方位角
vaz = x(5);       % 方位角速度
aaz = x(6);       % 方位角加速度
el = x(7);        % 俯仰角
vel = x(8);       % 俯仰角速度
ael = x(9);       % 俯仰角加速度

% 分析当前状态的运动特性 - 检测不合理的状态组合
is_inconsistent = false;

% 检查速度与加速度是否一致 (应该同向或加速度很小)
if abs(ar) > 0.5 && sign(vr) * sign(ar) < 0 && abs(vr) > 2.0
    % 速度大且加速度与速度方向相反，可能是不合理的状态
    is_inconsistent = true;
    fprintf('警告: 不一致的径向速度和加速度: vr=%.2f, ar=%.2f\n', vr, ar);
end

% 应用基本的线性状态转移
x_next = F * x;

% 检查并处理状态不一致问题
if is_inconsistent
    % 仅应用部分加速度变化，使运动更连贯
    x_next(3) = ar * 0.5;  % 减小加速度影响
    
    % 重新计算速度变化，确保一致性
    x_next(2) = vr + x_next(3) * F(2,3);
    
    fprintf('调整后: vr=%.2f, ar=%.2f\n', x_next(2), x_next(3));
end

% 确保状态在物理合理范围内
% 1. 处理角度的周期性
x_next(4) = wrapTo180(x_next(4));  % 方位角
x_next(7) = wrapTo180(x_next(7));  % 俯仰角

% 2. 位置约束
if x_next(1) < 0.1  % 距离不能为负或过小
    fprintf('警告: 预测距离过小(%.2f)，修正为最小值0.1m\n', x_next(1));
    x_next(1) = 0.1;
end

% 3. 速度方向修正 - 重点关注初始几帧
if frame_count <= 5
    % 强制确保前几帧的初始速度方向正确 (针对距离小于20米的情况)
    if x_next(1) < 20 && x_next(2) < 0
        x_next(2) = abs(x_next(2));  % 确保径向速度为正
        fprintf('初始阶段强制确保径向速度为正: %.2f\n', x_next(2));
    end
end

% 4. 速度约束
max_speed = 40.0;  % 最大允许速度 m/s
if abs(x_next(2)) > max_speed
    x_next(2) = sign(x_next(2)) * max_speed;
    fprintf('应用速度上限约束: |vr| = %.2f -> %.2f m/s\n', abs(x_next(2)), max_speed);
end

% 5. 加速度约束
max_accel = 20.0;  % 最大允许加速度 m/s²
if abs(x_next(3)) > max_accel
    x_next(3) = sign(x_next(3)) * max_accel;
    fprintf('应用加速度约束: |ar| = %.2f -> %.2f m/s²\n', abs(x_next(3)), max_accel);
end

% 6. 角速度约束
max_angular_vel = 30.0;  % 度/秒
if abs(x_next(5)) > max_angular_vel
    x_next(5) = sign(x_next(5)) * max_angular_vel;
end
if abs(x_next(8)) > max_angular_vel
    x_next(8) = sign(x_next(8)) * max_angular_vel;
end

% 7. 角加速度约束
max_angular_accel = 10.0;  % 度/秒²
if abs(x_next(6)) > max_angular_accel
    x_next(6) = sign(x_next(6)) * max_angular_accel;
end
if abs(x_next(9)) > max_angular_accel
    x_next(9) = sign(x_next(9)) * max_angular_accel;
end

% 运动学一致性检查 - 复杂的物理规则
% 检查径向距离变化与速度方向是否一致
if frame_count > 1 && abs(vr) > 1.0
    % 距离变化与速度方向应该基本一致
    expected_dr = vr * F(1,2) + 0.5 * ar * F(1,3)^2;
    actual_dr = x_next(1) - r;
    
    % 计算相对误差
    rel_error = abs((actual_dr - expected_dr) / (abs(expected_dr) + 1e-6));
    
    % 如果方向不一致且变化明显，可能是状态不合理
    if (sign(expected_dr) * sign(actual_dr) < 0 && abs(actual_dr) > 0.2) || rel_error > 0.5
        % 判断异常严重程度
        is_severe = rel_error > 1.0 || (sign(expected_dr) * sign(actual_dr) < 0 && abs(actual_dr) > 0.5);
        
        % 根据异常严重程度调整校正策略
        if is_severe
            % 严重异常：主要依赖预期变化，保留少量实际变化
            adjustment_factor = 0.85;  % 85%依赖预期变化
            new_r = r + expected_dr * adjustment_factor + actual_dr * (1 - adjustment_factor);
            
            % 根据修正的距离变化，调整速度
            corrected_vr = (new_r - r) / F(1,2) - 0.5 * ar * F(1,3)^2 / F(1,2);
            
            fprintf('严重运动学一致性异常: 预期dr=%.2f, 实际dr=%.2f, 修正后r=%.2f, vr=%.2f\n', ...
                  expected_dr, actual_dr, new_r, corrected_vr);
            
            x_next(1) = new_r;  % 更新距离
            x_next(2) = corrected_vr;  % 更新速度
        else
            % 轻度异常：使用加权平均校正
            adjustment_factor = 0.7;  % 调整强度
            
            % 计算基于距离变化推导的速度分量
            derived_vr = actual_dr / F(1,2) - 0.5 * ar * F(1,3) / F(1,2);
            
            % 融合当前速度估计和基于距离变化的速度估计
            corrected_vr = vr * (1 - adjustment_factor) + derived_vr * adjustment_factor;
            
            fprintf('轻度运动学一致性校正: 预期dr=%.2f, 实际dr=%.2f, 原始vr=%.2f, 校正后vr=%.2f\n', ...
                  expected_dr, actual_dr, vr, corrected_vr);
            
            x_next(2) = corrected_vr;  % 更新速度
        end
    end
end

% 添加加速度连续性检查
if frame_count > 2 && abs(x_next(3) - ar) > 5.0
    % 加速度变化过大，可能不合理
    max_accel_change = 5.0;  % 最大允许的加速度变化
    if abs(x_next(3) - ar) > max_accel_change
        % 限制加速度变化
        ar_change_limit = sign(x_next(3) - ar) * max_accel_change;
        x_next(3) = ar + ar_change_limit;
        fprintf('限制加速度变化: %.2f -> %.2f\n', ar, x_next(3));
    end
end

% 检查角速度与角度变化的一致性
if frame_count > 1
    % 方位角变化检查
    expected_daz = vaz * F(4,5) + 0.5 * aaz * F(4,6)^2;
    actual_daz = wrapTo180(x_next(4) - az);  % 考虑角度周期性
    
    if abs(wrapTo180(actual_daz - expected_daz)) > 5.0 && abs(vaz) > 2.0
        % 方位角变化与角速度不一致
        fprintf('方位角变化与角速度不一致: 预期=%.2f°, 实际=%.2f°\n', expected_daz, actual_daz);
        
        % 轻微调整角速度，使其更符合观察到的角度变化
        derived_vaz = actual_daz / F(4,5) - 0.5 * aaz * F(4,6) / F(4,5);
        x_next(5) = 0.7 * x_next(5) + 0.3 * derived_vaz;
    end
    
    % 俯仰角变化检查
    expected_del = vel * F(7,8) + 0.5 * ael * F(7,9)^2;
    actual_del = wrapTo180(x_next(7) - el);  % 考虑角度周期性
    
    if abs(wrapTo180(actual_del - expected_del)) > 5.0 && abs(vel) > 2.0
        % 俯仰角变化与角速度不一致
        fprintf('俯仰角变化与角速度不一致: 预期=%.2f°, 实际=%.2f°\n', expected_del, actual_del);
        
        % 轻微调整角速度，使其更符合观察到的角度变化
        derived_vel = actual_del / F(7,8) - 0.5 * ael * F(7,9) / F(7,8);
        x_next(8) = 0.7 * x_next(8) + 0.3 * derived_vel;
    end
end

return;
end 