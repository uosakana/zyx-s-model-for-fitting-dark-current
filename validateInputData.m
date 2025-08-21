function validateInputData(V, JD)
    % 检查输入数据是否为空
    if isempty(V) || isempty(JD)
        error('输入数据不能为空');
    end
    
    % 检查数据长度是否匹配
    if length(V) ~= length(JD)
        error('电压和电流数据长度必须相同');
    end
    
    % 检查是否包含非数值数据
    if any(isnan(V)) || any(isnan(JD))
        error('输入数据包含NaN值');
    end
    
    % 检查是否包含无穷大
    if any(isinf(V)) || any(isinf(JD))
        error('输入数据包含无穷大值');
    end
    
    % 检查数据类型
    if ~isnumeric(V) || ~isnumeric(JD)
        error('输入数据必须为数值类型');
    end
end