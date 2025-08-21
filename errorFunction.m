function err = errorFunction(x, data_V, data_JD, params, config)
    % 反缩放参数
    x_actual = x .* params.scaleFactors;
    
    % 计算模型预测值
    predicted = diodeModel(data_V, x_actual, config);
    
    % 计算误差，使用对数空间误差
    % 这对于跨多个数量级的数据特别有效
    err = zeros(size(data_JD));
    
    for i = 1:length(data_JD)
        % 确保数据为正值以便取对数（对于二极管电流，大多数时候是正的）
        actual_abs = abs(data_JD(i));
        pred_abs = abs(predicted(i));
        
        % 阈值，避免对近零值取对数
        threshold = 1e-12;
        
        % 处理小于阈值的情况
        if actual_abs < threshold || pred_abs < threshold
            % 对小值直接使用标准化误差
            err(i) = (predicted(i) - data_JD(i)) / max(1e-12, max(max(abs(data_JD))));
        else
            % 计算对数空间的误差
            log_actual = log10(actual_abs);
            log_pred = log10(pred_abs);
            err(i) = log_pred - log_actual;
            
            % 保持符号一致性
            if sign(predicted(i)) ~= sign(data_JD(i)) && abs(predicted(i)) > threshold && abs(data_JD(i)) > threshold
                err(i) = err(i) * 3; % 对符号不同的情况增加误差权重
            end
        end
        
        % 对负电压区域的错误增加权重
        if data_V(i) < -0.2
            % 对强负电压区域增加权重
            err(i) = err(i) * 5.0;
        elseif data_V(i) < 0
            % 对弱负电压区域增加权重
            err(i) = err(i) * 3.0;
        end
    end
    
    % 对误差施加权重，使小电流区域和大电流区域的拟合同等重要
    % 将数据电压区间分为几个部分，对每部分施加不同权重
    neg_voltages = unique(data_V(data_V < 0));
    n_neg_segments = 4; % 将负电压范围分为4段
    
    if ~isempty(neg_voltages)
        neg_segment_size = ceil(length(neg_voltages) / n_neg_segments);
        
        for i = 1:n_neg_segments
            start_idx = (i-1) * neg_segment_size + 1;
            end_idx = min(i * neg_segment_size, length(neg_voltages));
            
            if start_idx <= end_idx
                segment_voltages = neg_voltages(start_idx:end_idx);
                for v = segment_voltages'
                    idx = find(data_V == v);
                    % 对第一段（最负的电压区域）给予更高权重
                    if i == 1
                        err(idx) = err(idx) * 2;
                    end
                end
            end
        end
    end
    
    % 特别处理电压为零附近的点
    zero_idx = find(abs(data_V) < 0.05);
    if ~isempty(zero_idx)
        err(zero_idx) = err(zero_idx) * 2;
    end
end