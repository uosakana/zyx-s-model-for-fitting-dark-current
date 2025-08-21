function plotResults(V, JD_measured, fit_results, currents)
    % 创建图形窗口 - 详细分析
    figure('Position', [100, 100, 1200, 800]);
    
    % 第一个子图：I-V特性（对数坐标）及各电流分量
    subplot(2,2,[1,2]);
    
    % 绘制测量数据和拟合结果
    semilogy(V, abs(JD_measured), 'bo', 'DisplayName', '测量数据', 'MarkerSize', 6);
    hold on;
    semilogy(V, abs(currents.total), 'ro', 'DisplayName', '总拟合电流', 'MarkerSize', 6);
    semilogy(V, abs(currents.diode), 'b--', 'DisplayName', '二极管电流', 'LineWidth', 1.5);
    semilogy(V, abs(currents.ohmic), 'g--', 'DisplayName', '欧姆电流', 'LineWidth', 1.5);
    semilogy(V, abs(currents.nonohmic), 'm--', 'DisplayName', '非欧姆电流', 'LineWidth', 1.5);
    
    % 设置坐标轴
    xlabel('电压 (V)', 'FontSize', 12);
    ylabel('电流密度 (A)', 'FontSize', 12);
    title('二极管I-V特性曲线及电流分量（对数坐标）', 'FontSize', 14);
    grid on;
    legend('Location', 'best');
    
    % 设置合适的坐标轴范围
    xlim([min(V) max(V)]);
    ylim([min(abs(JD_measured))*0.1 max(abs(JD_measured))*10]);
    
    % 第二个子图：相对误差（线性坐标）
    subplot(2,2,3);
    relative_error = abs((currents.total - JD_measured) ./ (abs(JD_measured) + eps)) * 100;
    plot(V, relative_error, 'b.-', 'LineWidth', 1);
    avg_rel_err = mean(relative_error(V ~= 0));
    xlabel('电压 (V)');
    ylabel('相对误差 (%)');
    title('拟合相对误差');
    grid on;
    
    % 第三个子图：I-V特性（线性坐标）
    subplot(2,2,4);
    plot(V, JD_measured, 'bo', 'DisplayName', '测量数据');
    hold on;
    plot(V, currents.total, 'ro', 'DisplayName', '拟合结果');
    xlabel('电压 (V)');
    ylabel('电流密度 (A)');
    title('I-V特性曲线（线性坐标）');
    grid on;
    legend('Location', 'best');
    
    % 添加总体信息
    sgtitle('二极管特性拟合结果分析', 'FontSize', 14);
    
    % 在图上添加拟合参数信息
    annotation('textbox', [0.02, 0.02, 0.4, 0.05], ...
        'String', sprintf('最大相对误差: %.2f%%  平均相对误差: %.2f%%', ...
        max(relative_error), avg_rel_err), ...
        'EdgeColor', 'none');
    fprintf('二极管电流占比: %.2f%%\n', mean(abs(currents.diode ./ currents.total)) * 100);
    fprintf('欧姆电流占比: %.2f%%\n', mean(abs(currents.ohmic ./ currents.total)) * 100);
    fprintf('非欧姆电流占比: %.2f%%\n', mean(abs(currents.nonohmic ./ currents.total)) * 100);
end