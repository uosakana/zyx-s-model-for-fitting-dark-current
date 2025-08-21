function err = calculateError(x, data_V, data_JD, params, config)
    % 反缩放参数
    x_actual = x .* params.scaleFactors;
    
    % 计算模型预测值
    predicted = diodeModel(data_V, x_actual, config);
    
    % 对每个点计算相对误差
    relative_errors = zeros(size(data_JD));
    for i = 1:length(data_JD)
        % 使用相对误差，并给每个点相同的权重
        if abs(data_JD(i)) < 1e-12
            % 对零点附近的数据特殊处理
            relative_errors(i) = (predicted(i) - data_JD(i)) / 1e-12;
        else
            % 标准相对误差
            relative_errors(i) = (predicted(i) - data_JD(i)) / abs(data_JD(i));
        end
    end
    err = relative_errors;
end