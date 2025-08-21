function params = initializeParameters(data_V, data_JD, config)
    % 初始参数 [J0, Rs, Rsh, k]
    % J0: 饱和电流，主要影响正向特性
    % Rs: 串联电阻，影响高电流区域的线性度
    % Rsh: 并联电阻，主要影响低电压区域的泄漏电流
    % k: 非欧姆系数，对负电压区域的拟合尤为重要
    
    % 使用LambertW函数方法估计初始参数
    fprintf('使用Lambert W函数估计初始参数...\n');
    
    % 估计初始Rs值 - 使用正向高电压区域的斜率
    pos_idx = find(data_V > 0.2);
    if length(pos_idx) >= 2
        % 使用高电压区域估计Rs
        [~, max_idx] = max(data_V);
        if max_idx > 1
            Rs_est = abs((data_V(max_idx) - data_V(max_idx-1)) / (data_JD(max_idx) - data_JD(max_idx-1)));
        else
            Rs_est = 1e3; % 默认值
        end
    else
        Rs_est = 1e3; % 默认值
    end
    
    % 估计Rsh - 使用负电压区域的斜率
    neg_idx = find(data_V < -0.2);
    if length(neg_idx) >= 2
        % 使用线性拟合估计Rsh
        p = polyfit(data_V(neg_idx), data_JD(neg_idx), 1);
        Rsh_est = abs(1/p(1));
    else
        Rsh_est = 1e7; % 默认值
    end
    
    % 设置理想因子n
    n = config.physics.n;
    
    % 计算热电压
    V_th = config.physics.kb * config.physics.T / config.physics.q;
    
    % 使用LambertW方法求解J0
    % 对于测量数据中的每个点计算J0，然后取平均值
    J0_values = [];
    
    % 选择适合估算J0的正向区域数据点
    pos_idx = find(data_V > 0.1 & data_V < 0.25); % 选择中等正向偏置区域
    
    if ~isempty(pos_idx)
        for i = pos_idx'
            V = data_V(i);
            I = data_JD(i);
            
            % 使用Lambert W函数公式计算J0
            % I = (n*V_th/Rs) * LambertW( (Rs*Rsh*J0)/(n*V_th*(Rs+Rsh)) * exp( (Rsh*(Rs*J0+V))/(n*V_th*(Rs+Rsh)) ) ) - (Rsh*J0-V)/(Rs+Rsh)
            
            % 简化计算，忽略I_ph (光生电流)
            % 使用迭代法估算J0
            J0_est = 1e-9; % 初始猜测
            
            % 定义目标函数：已知I和V，求J0
            function_handle = @(J0) solveForJ0(J0, V, I, Rs_est, Rsh_est, n, V_th);
            
            % 使用优化算法求解J0
            options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
            J0_opt = lsqnonlin(function_handle, J0_est, 1e-12, 1e-6, options);
            
            J0_values = [J0_values; J0_opt];
        end
        
        % 计算J0的平均值
        J0_est = median(J0_values);
    else
        J0_est = 1e-9; % 默认值
    end
    
    % 估计非欧姆系数k - 根据负电压区域的数据
    k_est = 5e-7; % 初始值，后续会优化
    
    % 输出估计的参数
    fprintf('Lambert W估计的参数:\n');
    fprintf('J0 = %.6e A\n', J0_est);
    fprintf('Rs = %.6e Ohm\n', Rs_est);
    fprintf('Rsh = %.6e Ohm\n', Rsh_est);
    
    % 设置参数
    params.x0 = [J0_est, Rs_est, Rsh_est, k_est];
    
    % 参数范围设置
    params.ub = [1e-6, 1e4, 1e10, 1e-5];    % 上界
    params.lb = [1e-12, 1e1, 1e5, 1e-10];    % 下界
    
    % 确保初始值在范围内
    params.x0 = min(max(params.x0, params.lb), params.ub);
    
    % 缩放因子，使各参数量级相近
    params.scaleFactors = [1e-9, 1e3, 1e7, 1e-8];
end

% 辅助函数：求解J0的目标函数
function residual = solveForJ0(J0, V, I, Rs, Rsh, n, V_th)
    % 计算电流
    I_calc = calculateCurrentWithLambertW(V, J0, Rs, Rsh, n, V_th);
    
    % 计算残差
    residual = I_calc - I;
end

% 使用Lambert W函数计算电流
function I = calculateCurrentWithLambertW(V, J0, Rs, Rsh, n, V_th)
    % 使用Lambert W函数计算电流
    % 简化公式：I = (n*V_th/Rs) * LambertW( (Rs*Rsh*J0)/(n*V_th*(Rs+Rsh)) * exp( (Rsh*(Rs*J0+V))/(n*V_th*(Rs+Rsh)) ) ) - (Rsh*J0-V)/(Rs+Rsh)
    
    % 计算Lambert W函数的参数
    x = (Rs*Rsh*J0)/(n*V_th*(Rs+Rsh)) * exp((Rsh*(Rs*J0+V))/(n*V_th*(Rs+Rsh)));
    
    % 使用MATLAB的lambertw函数，如果没有Symbolic Math Toolbox，可以使用近似实现
    if exist('lambertw', 'file')
        w = lambertw(x);
    else
        w = approximateLambertW(x);
    end
    
    % 计算电流
    I = (n*V_th/Rs) * w - (Rsh*J0-V)/(Rs+Rsh);
end

% Lambert W函数的近似实现
function w = approximateLambertW(x)
    % 对于x >= 0的近似计算
    if x < 0
        w = 0; % 负值没有定义，返回0作为默认值
    elseif x == 0
        w = 0;
    elseif x < 1
        % 使用级数展开近似
        w = x * (1 - x + 1.5*x^2 - 2.667*x^3 + 5.208*x^4);
    else
        % 使用迭代法
        % 初始猜测
        if x < 3
            w = 0.5;
        else
            w = log(x) - log(log(x));
        end
        
        % 牛顿迭代
        for i = 1:10
            exp_w = exp(w);
            w_next = w - (w*exp_w - x)/(exp_w + w*exp_w);
            
            % 检查收敛性
            if abs(w_next - w) < 1e-10
                w = w_next;
                break;
            end
            
            w = w_next;
        end
    end
end