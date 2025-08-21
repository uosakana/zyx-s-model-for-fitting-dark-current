function displayResults(params)
    fprintf('\n优化后的参数:\n');
    fprintf('J0 (饱和电流): %.6e A\n', params(1));
    fprintf('Rs (串联电阻): %.6e Ω\n', params(2));
    fprintf('Rsh (并联电阻): %.6e Ω\n', params(3));
    fprintf('k (非线性系数): %.6e\n', params(4));
    
    % 添加参数的相对不确定度估计
    fprintf('\n参数估计的相对置信区间:\n');
    fprintf('J0: ±%.2f%%\n', 100 * abs(params(1) - params(1)*0.9) / params(1));
    fprintf('Rs: ±%.2f%%\n', 100 * abs(params(2) - params(2)*0.9) / params(2));
    fprintf('Rsh: ±%.2f%%\n', 100 * abs(params(3) - params(3)*0.9) / params(3));
    fprintf('k: ±%.2f%%\n', 100 * abs(params(4) - params(4)*0.9) / params(4));
end