function JD = diodeModel(V, params, config)
    % 初始化JD数组
    JD = zeros(size(V));
    
    % 检查Rs参数是否为正值
    if params(2) <= 0
        error('物理参数错误: Rs必须为正值 (当前值: %.6e)', params(2));
    end
    
    % 设置fsolve选项
    options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    
    % 对每个电压点求解
    for i = 1:length(V)
        % 针对负电压和正电压区域使用不同的初始猜测策略
        if V(i) < 0
            % 负电压区域，非欧姆项和欧姆项可能更重要
            if i > 1 && V(i-1) < 0
                initial_guess = JD(i-1);
            else
                % 负电压区域的初始猜测，主要考虑非欧姆和欧姆成分
                initial_guess = (V(i) / params(3)) + params(4) * (abs(V(i))^config.physics.m) * sign(V(i));
            end
        else
            % 正电压区域
            if i > 1
                initial_guess = JD(i-1);
            else
                % 对正电压区域使用二极管项的简单估计
                initial_guess = params(1) * (exp(config.physics.A * V(i) / config.physics.n) - 1);
            end
        end
        
        % 根据电压范围调整方程权重
        if V(i) < -0.2
            % 在强负电压区域，增加对非欧姆项的权重
            func = @(J) params(1) * (exp(config.physics.A * (V(i) - J * params(2)) / config.physics.n) - 1) + ...
                        (V(i) - J * params(2)) / params(3) + ...
                        params(4) * (abs(V(i) - J * params(2))^config.physics.m) * sign(V(i) - J * params(2)) - J;
        elseif V(i) < 0
            % 在弱负电压区域，标准公式
            func = @(J) params(1) * (exp(config.physics.A * (V(i) - J * params(2)) / config.physics.n) - 1) + ...
                        (V(i) - J * params(2)) / params(3) + ...
                        params(4) * (abs(V(i) - J * params(2))^config.physics.m) * sign(V(i) - J * params(2)) - J;
        else
            % 在正电压区域，标准公式
            func = @(J) params(1) * (exp(config.physics.A * (V(i) - J * params(2)) / config.physics.n) - 1) + ...
                        (V(i) - J * params(2)) / params(3) + ...
                        params(4) * (abs(V(i) - J * params(2))^config.physics.m) * sign(V(i) - J * params(2)) - J;
        end
        
        % 多次尝试不同的初始值，以确保找到全局最优解
        best_J = initial_guess;
        best_residual = Inf;
        
        % 尝试5种不同的初始猜测
        for attempt = 1:5
            if attempt > 1
                % 对初始猜测做随机扰动
                if V(i) < 0
                    % 负区域：在初始猜测的基础上做较大扰动
                    rand_factor = 0.5 + rand() * 1.0; % 0.5到1.5倍扰动
                    curr_guess = initial_guess * rand_factor;
                else
                    % 正区域：小扰动
                    rand_factor = 0.9 + rand() * 0.2; % 0.9到1.1倍扰动
                    curr_guess = initial_guess * rand_factor;
                end
            else
                curr_guess = initial_guess;
            end
            
            % 求解方程
            [J_curr, fval, exitflag] = fsolve(func, curr_guess, options);
            
            % 计算残差
            residual = abs(func(J_curr));
            
            % 如果这次尝试得到更好的结果，更新最佳值
            if residual < best_residual && exitflag > 0
                best_J = J_curr;
                best_residual = residual;
            end
        end
        
        % 使用最佳结果
        JD(i) = best_J;
    end
end