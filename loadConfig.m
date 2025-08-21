function config = loadConfig()
    config.physics = struct(...
        'q', 1.602e-19, ...    % 电子电荷 (C)
        'kb', 1.38e-23, ...    % 波尔兹曼常数 (J/K)
        'T', 300, ...          % 温度 (K)
        'n', 1.8, ...          % Ideality factor - 增大以更好适应负区域
        'm', 2.4 ...           % Exponent for recombination - 对非欧姆项的指数调整
    );
    
    config.physics.A = config.physics.q / (config.physics.kb * config.physics.T);
    
    % 设置不同电压区域的参数
    config.fitting = struct(...
        'neg_voltage_threshold', -0.2, ... % 负电压区域阈值
        'pos_voltage_threshold', 0.1 ...   % 正电压区域阈值
    );
end