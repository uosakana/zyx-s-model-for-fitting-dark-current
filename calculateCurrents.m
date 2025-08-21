function currents = calculateCurrents(V, x, config)
    % 计算拟合的总电流
    currents.total = diodeModel(V, x, config);
    
    % 计算电压降
    V_drop = V - currents.total .* x(2);
    
    % 计算各分量电流
    currents.diode = x(1) * (exp(config.physics.A * V_drop / config.physics.n) - 1);
    currents.ohmic = V_drop / x(3);
    currents.nonohmic = x(4) * (abs(V_drop).^config.physics.m) .* sign(V_drop);
    
    % 计算总电流（验证用）
    currents.sum = currents.diode + currents.ohmic + currents.nonohmic;
end