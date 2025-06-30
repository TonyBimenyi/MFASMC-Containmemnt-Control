clear all;
close all;

% Define parameters
file_path = 'dataplot.txt'; % Replace with your file's path

% Load the content of the data
file_content = fileread(file_path);

% Extract the CH3, CH1, CH2, CH4 data using regular expressions
ch3_data_match = regexp(file_content, 'CH3_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch1_data_match = regexp(file_content, 'CH1_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch2_data_match = regexp(file_content, 'CH2_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch4_data_match = regexp(file_content, 'CH4_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');
ch5_data_match = regexp(file_content, 'CH4_Data_OutPut\[\d+\]=\{(.*?)\};', 'tokens', 'dotall');

if ~isempty(ch3_data_match) && ~isempty(ch1_data_match) && ~isempty(ch2_data_match) && ~isempty(ch4_data_match)
    % Convert the extracted data to numerical arrays
    ch3_data = str2num(ch3_data_match{1}{1});
    ch1_data = str2num(ch1_data_match{1}{1});
    ch2_data = str2num(ch2_data_match{1}{1});
    ch4_data = str2num(ch4_data_match{1}{1});
    ch5_data = str2num(ch5_data_match{1}{1});
    
    % Calculate distributed errors (xi)
    xi1 = ch2_data - ch1_data; % Error for Motor 1
    xi2 = ch2_data - ch3_data; % Error for Motor 2
    xi3 = ch2_data - ch4_data; % Error for Motor 3
    
    % Calculate MSE for each error
    mse1 = mean(xi1.^2); % MSE for Motor 1
    mse2 = mean(xi2.^2); % MSE for Motor 2
    mse3 = mean(xi3.^2); % MSE for Motor 3
    
    % Display MSE values
    fprintf('MSE for Motor 1 (CH2 - CH1): %.10e\n', mse1);
    fprintf('MSE for Motor 2 (CH2 - CH3): %.10e\n', mse2);
    fprintf('MSE for Motor 3 (CH2 - CH4): %.10e\n', mse3);
    
    % Plot the data
    figure('Position', [100, 100, 800, 500]);
    plot(ch2_data, 'b--', 'LineWidth', 2, 'MarkerSize', 4.5, 'DisplayName', 'Leader 1');
    hold on;

    plot(ch3_data, '--', 'LineWidth', 2, 'Color', 'black', 'MarkerSize', 2.5, 'DisplayName', 'Leader 2');

    plot(ch1_data, '-.m', 'LineWidth', 2, 'MarkerSize', 5.5, 'DisplayName', 'Motor1');
  
    plot(ch4_data, 'r-', 'LineWidth', 2, 'MarkerSize', 2.5, 'DisplayName', 'Motor 2');

    plot(ch5_data, '--', 'LineWidth', 2, 'Color', 'green', 'MarkerSize', 5.5, 'DisplayName', 'Motor 3');
    hold off;
      % Add legend
      legend('show','location','north','orientation','horizontal');
    
    % Set plot limits and labels
    xlim([600, 1400]);
    ylim([20, 85]);
    % xlabel('Time step', 'FontSize', 14);
    % ylabel('Output', 'FontSize', 14);
    
    % ylim([-1.5 3]); % Y-axis limits for Agent 4
    set(gca, 'FontSize', 14);
    
  
    xlabel('Time step(k)', 'FontSize', 14);
    ylabel('Output', 'FontSize', 14);

    zoom_x_start = 940; % Start of zoomed x-range
    zoom_x_end = 980; % End of zoomed x-range
    axes('Position', [0.23,0.577,0.30,0.25]);
    box on; hold on;
    plot(ch2_data, '--b', 'LineWidth', 2.5);
    plot(ch3_data, '--k', 'LineWidth', 2.5);
    plot(ch1_data, '-.m', 'LineWidth', 2.5);

    plot(ch4_data, '-r', 'LineWidth', 2.5);
    plot(ch4_data, '--g', 'LineWidth', 2.5);
    xlim([zoom_x_start zoom_x_end]);

    % yticks([-0.4,0,0.2]);
    set(gca, 'FontSize', 14);
    
    
    % Save the plot
    print('-dpng', 'motor_plot.png');
else
    disp('CH3_Data_OutPut[1707] not found in the file.');
end