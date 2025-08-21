function main()
    % 加载配置和数据
    config = loadConfig();
    [data_V, data_JD] = loadData();
    
    % 验证输入数据
    validateInputData(data_V, data_JD);
    
    % 询问是否从历史参数文件加载初始参数
    use_history = input('是否使用历史参数文件作为初始参数? (y/n): ', 's');
    if strcmpi(use_history, 'y')
        % 列出可用的参数文件
        mat_files = dir('adjusted_params_*.mat');
        txt_files = dir('adjusted_params_*.txt');
        
        if isempty(mat_files)
            fprintf('未找到历史参数文件，将使用Lambert W函数估计初始参数\n');
            params = initializeParameters(data_V, data_JD, config);
        else
            % 显示可用文件
            fprintf('可用的参数文件:\n');
            for i = 1:length(mat_files)
                fprintf('%d: %s\n', i, mat_files(i).name);
            end
            
            % 让用户选择文件
            file_idx = input('请选择要加载的文件编号(输入0取消): ');
            if file_idx > 0 && file_idx <= length(mat_files)
                % 加载选定的文件
                load_file = mat_files(file_idx).name;
                loaded_data = load(load_file);
                
                % 提取参数
                if isfield(loaded_data, 'params')
                    fprintf('从文件 %s 加载参数\n', load_file);
                    
                    % 创建参数结构体
                    params = struct();
                    params.x0 = loaded_data.params;
                    
                    % 参数范围设置
                    params.ub = [1e-6, 1e4, 1e10, 1e-5];    % 上界
                    params.lb = [1e-12, 1e1, 1e5, 1e-10];    % 下界
                    
                    % 确保初始值在范围内
                    params.x0 = min(max(params.x0, params.lb), params.ub);
                    
                    % 缩放因子
                    params.scaleFactors = [1e-9, 1e3, 1e7, 1e-8];
                    
                    % 显示加载的参数
                    fprintf('加载的参数:\n');
                    fprintf('J0 = %.6e A\n', params.x0(1));
                    fprintf('Rs = %.6e Ohm\n', params.x0(2));
                    fprintf('Rsh = %.6e Ohm\n', params.x0(3));
                    fprintf('k = %.6e\n', params.x0(4));
                else
                    fprintf('文件格式错误，将使用Lambert W函数估计初始参数\n');
                    params = initializeParameters(data_V, data_JD, config);
                end
            else
                fprintf('取消加载参数文件，将使用Lambert W函数估计初始参数\n');
                params = initializeParameters(data_V, data_JD, config);
            end
        end
    else
        % 初始化参数 - 传入数据，用于Lambert W函数估计初始参数
        params = initializeParameters(data_V, data_JD, config);
    end
    
    % 执行拟合
    [optimized_params, fit_results] = performFitting(data_V, data_JD, params, config);
    
    % 计算各分量电流
    currents = calculateCurrents(data_V, optimized_params, config);
    
    % 绘制结果
    plotResults(data_V, data_JD, fit_results, currents);
    
    % 保存拟合结果和图像
    save_results = input('是否保存拟合结果? (y/n): ', 's');
    if strcmpi(save_results, 'y')
        saveResults(data_V, data_JD, optimized_params, fit_results, currents);
    end
    
    % 输出结果
    displayResults(optimized_params);
    
    % 询问用户是否满意拟合结果，如果不满意则进入交互式调整模式
    interactive_adjust = input('是否进入交互式参数调整模式? (y/n): ', 's');
    if strcmpi(interactive_adjust, 'y')
        [refined_params, refined_fit, refined_currents] = interactiveParameterAdjustment(data_V, data_JD, optimized_params, config);
        
        % 绘制新的结果
        figure;
        plotResults(data_V, data_JD, refined_fit, refined_currents);
        
        % 输出调整后的结果
        displayResults(refined_params);
        
        % 保存调整后的参数
        saveAdjustedResults(refined_params, data_V, refined_currents);
    end
end

% 交互式参数调整函数
function [adjusted_params, fit_results, final_currents] = interactiveParameterAdjustment(data_V, data_JD, initial_params, config)
    % 复制初始参数
    adjusted_params = initial_params;
    adjustment_factor = 1.0;
    
    % 确保初始参数物理合理性
    if adjusted_params(2) <= 0  % Rs必须为正
        fprintf('警告: 初始Rs为负值或零，已自动调整为正值\n');
        adjusted_params(2) = 10; % 使用一个合理的默认值
    end
    
    % 计算初始拟合和误差
    %fit_results.JD = diodeModel(data_V, adjusted_params, config);
    currents = calculateCurrents(data_V, adjusted_params, config);
    fit_results.JD = currents.total;
    errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    nz_idx = data_V ~= 0;
    avg_error = mean(errors(nz_idx));

    % 创建实时更新的图表
    figure('Name', '交互式参数调整', 'Position', [100, 100, 1200, 800]);
    
    subplot(2,2,[1,2]);
    h_data = semilogy(data_V, abs(data_JD), 'bo', 'DisplayName', '测量数据', 'MarkerSize', 6);
    hold on;
    %h_fit = semilogy(data_V, abs(fit_results.JD), 'ro', 'DisplayName', '拟合结果');
    h_fit = semilogy(data_V, abs(currents.total), 'ro', 'DisplayName', '拟合结果', 'MarkerSize', 6);
    h_diode = semilogy(data_V, abs(currents.diode), 'b--', 'DisplayName', '二极管电流');
    h_ohmic = semilogy(data_V, abs(currents.ohmic), 'g--', 'DisplayName', '欧姆电流');
    h_nonohmic = semilogy(data_V, abs(currents.nonohmic), 'm--', 'DisplayName', '非欧姆电流');
    xlabel('电压 (V)');
    ylabel('电流密度 (A)');
    title('电流-电压特性 (对数尺度)');
    legend('Location', 'best');
    grid on;
    
    xlim([min(data_V) max(data_V)]);
    ylim([min(abs(data_JD))*0.1 max(abs(data_JD))*10]);

    % 第二行左侧：相对误差
    subplot(2,2,3);
    h_error = plot(data_V, errors, 'b.-');
    xlabel('电压 (V)');
    ylabel('相对误差 (%)');
    %title(sprintf('拟合误差 (平均: %.2f%%)', mean(errors)));
    title(sprintf('拟合误差 (平均: %.2f%%)', avg_error));
    grid on;

    subplot(2,2,4);
    h_data_lin = plot(data_V, data_JD, 'bo', 'DisplayName', '测量数据');
    hold on;
    h_fit_lin = plot(data_V, currents.total, 'ro', 'DisplayName', '拟合结果');
    xlabel('电压 (V)');
    ylabel('电流密度 (A)');
    title('I-V特性曲线（线性坐标）');
    legend('Location', 'best');
    grid on;
    xlim([min(data_V) max(data_V)]);
    % 显示当前参数值
    annotation('textbox', [0.01, 0.01, 0.98, 0.08], ...
        'String', sprintf('J0: %.2e A   Rs: %.2e Ohm   Rsh: %.2e Ohm   k: %.2e   调整步长: %.2f', ...
        adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4), adjustment_factor), ...
        'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    % 持续调整直到用户满意
    while true
        % 显示调整选项
        fprintf('\n当前参数: J0=%.2e, Rs=%.2e, Rsh=%.2e, k=%.2e\n', ...
            adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4));
        fprintf('平均相对误差: %.2f%%\n', avg_error);
        fprintf('\n参数调整选项:\n');
        fprintf('1: 增加 J0  2: 减少 J0\n');
        fprintf('3: 增加 Rs  4: 减少 Rs\n');
        fprintf('5: 增加 Rsh 6: 减少 Rsh\n');
        fprintf('7: 增加 k   8: 减少 k\n');
        fprintf('9: 更改调整步长 (当前: %.2f)\n', adjustment_factor);
        fprintf('0: 结束调整并保存结果\n');
        
        % 获取用户输入并确保是数值类型
        choice_str = input('请选择操作 (0-9): ', 's');
        choice = str2double(choice_str);
        
        % 检查是否为有效数字输入
        if isnan(choice)
            fprintf('请输入有效的数字(0-9)\n');
            continue;
        end
        
        if choice == 0
            break;
        elseif choice == 9
            % 调整步长
            new_factor_str = input(sprintf('输入新的调整步长 (当前: %.2f): ', adjustment_factor), 's');
            new_factor = str2double(new_factor_str);
            if ~isnan(new_factor) && new_factor > 0
                adjustment_factor = new_factor;
            else
                fprintf('输入无效，保持当前步长: %.2f\n', adjustment_factor);
            end
            continue;
        elseif choice >= 1 && choice <= 8
            % 确定要调整的参数索引
            param_idx = ceil(choice / 2);
            
            % 确定调整方向
            if mod(choice, 2) == 1
                direction = 1;
            else
                direction = -1;
            end
            
            % 计算调整量
            delta = adjusted_params(param_idx) * 0.1 * adjustment_factor * direction;
            
            % 更新参数
            adjusted_params(param_idx) = adjusted_params(param_idx) + delta;
            
            % 确保参数在合理范围内
            if param_idx == 1 % J0
                adjusted_params(param_idx) = max(1e-12, adjusted_params(param_idx));
            elseif param_idx == 2 % Rs - 特别强调必须为正值
                adjusted_params(param_idx) = max(1, adjusted_params(param_idx));
                if adjusted_params(param_idx) <= 0
                    fprintf('警告: Rs不能为负值或零。已调整为正值。\n');
                    adjusted_params(param_idx) = 1; % 确保为正值
                end
            elseif param_idx == 3 % Rsh
                adjusted_params(param_idx) = max(1e4, adjusted_params(param_idx));
            elseif param_idx == 4 % k
                adjusted_params(param_idx) = max(1e-10, adjusted_params(param_idx));
            end
            
            % 重新计算拟合和误差
            %fit_results.JD = diodeModel(data_V, adjusted_params, config);
            currents = calculateCurrents(data_V, adjusted_params, config);
            fit_results.JD = currents.total;
            errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
            
            avg_error = mean(errors(nz_idx));

            % 更新图表
            %set(h_fit, 'YData', abs(fit_results.JD));
            set(h_fit, 'YData', abs(currents.total));
            set(h_diode, 'YData', abs(currents.diode));
            set(h_ohmic, 'YData', abs(currents.ohmic));
            set(h_nonohmic, 'YData', abs(currents.nonohmic));
            set(h_error, 'YData', errors);
            set(h_fit_lin, 'YData', currents.total);
            title(subplot(2,2,3), sprintf('拟合误差 (平均: %.2f%%)', avg_error));
            % 更新参数显示
            delete(findall(gcf, 'Type', 'annotation'));
            annotation('textbox', [0.01, 0.01, 0.98, 0.08], ...
                'String', sprintf('J0: %.2e A   Rs: %.2e Ohm   Rsh: %.2e Ohm   k: %.2e   调整步长: %.2f', ...
                adjusted_params(1), adjusted_params(2), adjusted_params(3), adjusted_params(4), adjustment_factor), ...
                'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
            
            drawnow;
        else
            fprintf('无效的选择，请输入0-9之间的数字\n');
        end
    end
    
    % 最终检查确保Rs为正值
    if adjusted_params(2) <= 0
        fprintf('警告: 最终Rs为负值或零，已调整为正值\n');
        adjusted_params(2) = 1; % 设置为一个合理的小正值
    end
    
    % 计算最终拟合结果
    %fit_results.JD = diodeModel(data_V, adjusted_params, config);
    final_currents = calculateCurrents(data_V, adjusted_params, config);
    fit_results.JD = final_currents.total;
    fit_results.resnorm = sum(((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)).^2);
end

% 保存调整后的结果
function saveAdjustedResults(params, V, currents)
    % 生成时间戳
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    % 保存参数到MAT文件
    mat_filename = sprintf('adjusted_params_%s.mat', timestamp);
    save(mat_filename, 'params');
    fprintf('参数已保存到文件: %s\n', mat_filename);

    % 同时导出参数为文本文件
    txt_filename = sprintf('adjusted_params_%s.txt', timestamp);
    fid = fopen(txt_filename, 'w');
    fprintf(fid, 'J0 = %.6e A\n', params(1));
    fprintf(fid, 'Rs = %.6e Ohm\n', params(2));
    fprintf(fid, 'Rsh = %.6e Ohm\n', params(3));
    fprintf(fid, 'k = %.6e\n', params(4));
    fclose(fid);
    fprintf('参数已导出为文本文件: %s\n', txt_filename);
    
    % 导出电流数据到Excel
    T = table(V(:), currents.total(:), currents.diode(:), ...
              currents.ohmic(:), currents.nonohmic(:), ...
              'VariableNames', {'Voltage','Total','Diode','Ohmic','Nonohmic'});
    xlsx_filename = sprintf('adjusted_currents_%s.xlsx', timestamp);
    writetable(T, xlsx_filename);
    fprintf('调整后的电流已导出到Excel文件: %s\n', xlsx_filename);
end

% 保存拟合结果
function saveResults(data_V, data_JD, params, fit_results, currents)
    % 生成时间戳
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    
    % 保存数据和拟合结果
    results_filename = sprintf('fit_results_%s.mat', timestamp);
    save(results_filename, 'data_V', 'data_JD', 'params', 'fit_results', 'currents');
    fprintf('拟合结果已保存到文件: %s\n', results_filename);
    
    % 保存图形
    fig_filename = sprintf('fit_plot_%s.png', timestamp);
    saveas(gcf, fig_filename);
    fprintf('拟合图形已保存到文件: %s\n', fig_filename);
    
    % 导出详细数据到CSV
    csv_filename = sprintf('fit_data_%s.csv', timestamp);
    fid = fopen(csv_filename, 'w');
    fprintf(fid, 'Voltage(V),Measured_Current(A),Fitted_Current(A),Diode_Current(A),Ohmic_Current(A),Nonohmic_Current(A),Relative_Error(%%)\n');
    
    % 计算相对误差
    rel_errors = abs((fit_results.JD - data_JD) ./ (abs(data_JD) + eps)) * 100;
    
    % 写入数据
    for i = 1:length(data_V)
        fprintf(fid, '%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.2f\n', ...
            data_V(i), data_JD(i), fit_results.JD(i), ...
            currents.diode(i), currents.ohmic(i), currents.nonohmic(i), ...
            rel_errors(i));
    end
    fclose(fid);
    fprintf('拟合数据已导出到CSV文件: %s\n', csv_filename);
    
    % 导出拟合参数到文本文件
    params_filename = sprintf('fit_params_%s.txt', timestamp);
    fid = fopen(params_filename, 'w');
    fprintf(fid, '拟合参数:\n');
    fprintf(fid, 'J0 = %.6e A\n', params(1));
    fprintf(fid, 'Rs = %.6e Ohm\n', params(2));
    fprintf(fid, 'Rsh = %.6e Ohm\n', params(3));
    fprintf(fid, 'k = %.6e\n\n', params(4));
    
    % 添加误差统计
    fprintf(fid, '拟合误差统计:\n');
    %fprintf(fid, '平均相对误差: %.2f%%\n', mean(rel_errors));
    avg_rel_error = mean(rel_errors(data_V ~= 0));
    fprintf(fid, '平均相对误差: %.2f%%\n', avg_rel_error);
    fprintf(fid, '最大相对误差: %.2f%%\n', max(rel_errors));
    fprintf(fid, '中位相对误差: %.2f%%\n', median(rel_errors));
    
    % 计算不同电压区域的误差
    neg_idx = find(data_V < 0);
    pos_idx = find(data_V > 0);
    if ~isempty(neg_idx)
        fprintf(fid, '负电压区域平均相对误差: %.2f%%\n', mean(rel_errors(neg_idx)));
    end
    if ~isempty(pos_idx)
        fprintf(fid, '正电压区域平均相对误差: %.2f%%\n', mean(rel_errors(pos_idx)));
    end
    
    fclose(fid);
    fprintf('拟合参数和统计信息已保存到文件: %s\n', params_filename);
end
