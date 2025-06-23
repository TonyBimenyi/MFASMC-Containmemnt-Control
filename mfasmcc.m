clc; clear;

% Parameters
rho = 10.5;
eta = 5;
lamda = 200;
mu = 80.5;
epsilon = 1e-5;
alpha = 1;
T = 0.1;
gamma1 = 0.05;
gamma2 = 0.05;
gamma3 = 0.05;
gamma4 = 0.05;
beta = 10;
sigma = 95;
tau = 1e-5;
nena = 1e-5;
rT = 1024;
L = 200;
m = 500;
n = 600;

% Lambda values for subplot
lambda_values = [75, 100, 120, 150];

% Time vector for plotting
k = 1:m;
t = (k-1) * T;

% First figure: 2x2 subplot for different lambda values with sinusoidal leaders
figure('Name', 'Outputs for Different \lambda (Piecewise Leaders)', 'Position', [100, 100, 1350, 750]);
w5_eq = '0.5 + 0.25 * sin(0.1 * t)';
w6_eq = '0.5 + 0.25 * sin(0.1 * t + pi/2)';
font_size = 14;

leader_type = 'sinusoidal';

for i = 1:4
    lambda = lambda_values(i);
    [y1, y2, y3, y4, w5, w6] = run_simulation(lambda, w5_eq, w6_eq, leader_type, m, T, rho, eta, mu, epsilon, alpha, gamma1, gamma2, gamma3, gamma4, beta, sigma, tau, nena, rT);
    
    subplot(2, 2, i);
    plot(k, y1(1:m), 'b-', 'LineWidth', 2, 'DisplayName', 'y1(k)'); hold on;
    plot(k, y2(1:m), 'r-', 'LineWidth', 2, 'DisplayName', 'y2(k)');
    plot(k, y3(1:m), 'g-', 'LineWidth', 2, 'DisplayName', 'y3(k)');
    plot(k, y4(1:m), 'm-', 'LineWidth', 2, 'DisplayName', 'y4(k)');
    plot(k, w5, 'c--', 'LineWidth', 2, 'DisplayName', 'w5(k)');
    plot(k, w6, 'k--', 'LineWidth', 2, 'DisplayName', 'w6(k)');
    ylim([0 0.9]);
    hold off;
    title(['\lambda = ', num2str(lambda)]);
    xlabel('Iteration (k)');
    ylabel('Value');
    legend('show', 'Location', 'north', 'Orientation', 'horizontal');
    set(gca, 'FontSize', font_size);
    grid off;
end


% === Setup ===
lambda_values = [75, 100, 115, 150];
t = 1:m+1;         % Time for y outputs
k = 1:m;           % Time for leader signals
font_size = 14;

zoom_x_start = 38;  % Start of zoomed-in x-range
zoom_x_end = 40;    % End of zoomed-in x-range

% === Create figure ===
figure('Name', 'Outputs for Different \lambda (Piecewise Leaders)', 'Position', [100, 100, 1350, 750]);

for i = 1:4
    lambda = lambda_values(i);
    w5_eq = ''; w6_eq = '';
    leader_type = 'piecewise';

    [y1, y2, y3, y4, w5, w6] = run_simulation(lambda, w5_eq, w6_eq, leader_type, ...
        m, T, rho, eta, mu, epsilon, alpha, gamma1, gamma2, gamma3, gamma4, ...
        beta, sigma, tau, nena, rT);

    subplot(2, 2, i);
    
    plot(t, y1(1:m+1), 'b-', 'LineWidth', 2, 'DisplayName', 'y1'); hold on;
    plot(t, y2(1:m+1), 'r-', 'LineWidth', 2, 'DisplayName', 'y2');
    plot(t, y3(1:m+1), 'g-', 'LineWidth', 2, 'DisplayName', 'y3');
    plot(t, y4(1:m+1), 'm-', 'LineWidth', 2, 'DisplayName', 'y4');
    plot(k, w5, 'c--', 'LineWidth', 2, 'DisplayName', 'w5'); 
    plot(k, w6, 'k--', 'LineWidth', 2, 'DisplayName', 'w6');
    hold off;

    title(['\lambda = ', num2str(lambda)], 'FontSize', font_size);
    xlabel('Time Step (k)', 'FontSize', font_size);
    ylabel('Output', 'FontSize', font_size);
    legend('show', 'Location', 'north', 'Orientation', 'horizontal');
    set(gca, 'FontSize', font_size);
    xlim([0 m]);
    ylim([0 1.9]);
    grid off;

    

    % % === Zoomed Axes Positions ===
    % if i == 1
    %     zoom_pos = [0.21, 0.62, 0.13, 0.09];
    % elseif i == 2
    %     zoom_pos = [0.62, 0.76, 0.13, 0.10];
    % elseif i == 3
    %     zoom_pos = [0.18, 0.287, 0.13, 0.10];
    % else
    %     zoom_pos = [0.62, 0.287, 0.13, 0.10];
    % end

    % % === Zoomed Axes ===
    % ax_zoom = axes('Position', zoom_pos);
    % box on; hold on;
    % plot(k, w5, 'c--', 'LineWidth', 2);
    % plot(k, w6, 'k--', 'LineWidth', 2);
    % plot(t, y1(1:m+1), 'b-', 'LineWidth', 2);
    % plot(t, y2(1:m+1), 'r-', 'LineWidth', 2);
    % plot(t, y3(1:m+1), 'g-', 'LineWidth', 2);
    % plot(t, y4(1:m+1), 'm-', 'LineWidth', 2);
    % xlim([zoom_x_start zoom_x_end]);
    % ylim([0 4.5])
    % set(gca, 'FontSize', font_size);
end


% Third figure: Sinusoidal leader signals
% figure('Name', 'System Outputs with Sinusoidal Leader Signals');
% lambda = 100; % Fixed lambda
% w5_eq = '0.5 + 0.25 * sin(0.1 * t)';
% w6_eq = '0.5 + 0.25 * sin(0.1 * t + pi/2)';
% leader_type = 'sinusoidal';
% [y1, y2, y3, y4, w5, w6] = run_simulation(lambda, w5_eq, w6_eq, leader_type, m, T, rho, eta, mu, epsilon, alpha, gamma1, gamma2, gamma3, gamma4, beta, sigma, tau, nena, rT);

% plot(k, y1(1:m), 'b-', 'LineWidth', 2, 'DisplayName', 'y1'); hold on;
% plot(k, y2(1:m), 'r-', 'LineWidth', 2, 'DisplayName', 'y2');
% plot(k, y3(1:m), 'g-', 'LineWidth', 2, 'DisplayName', 'y3');
% plot(k, y4(1:m), 'm-', 'LineWidth', 2, 'DisplayName', 'y4');
% plot(k, w5, 'c--', 'LineWidth', 2, 'DisplayName', 'w5 (sin, f=0.1)');
% plot(k, w6, 'k--', 'LineWidth', 2, 'DisplayName', 'w6 (sin, f=0.1, \phi=\pi/2)');
% hold off;
% title('System Outputs with Sinusoidal Leader Signals');
% xlabel('Iteration (k)');
% ylabel('Value');
% legend('show');
% grid on;

% Function to run simulation with given parameters
function [y1, y2, y3, y4, w5, w6] = run_simulation(lambda, w5_eq, w6_eq, leader_type, m, T, rho, eta, mu, epsilon, alpha, gamma1, gamma2, gamma3, gamma4, beta, sigma, tau, nena, rT)
    % Initialize arrays
    phi1 = zeros(m+1, 1); phi2 = zeros(m+1, 1); phi3 = zeros(m+1, 1); phi4 = zeros(m+1, 1);
    mfa1 = zeros(m+1, 1); mfa2 = zeros(m+1, 1); mfa3 = zeros(m+1, 1); mfa4 = zeros(m+1, 1);
    sm1 = zeros(m, 1); sm2 = zeros(m, 1); sm3 = zeros(m, 1); sm4 = zeros(m, 1);
    y1 = zeros(m+1, 1); y2 = zeros(m+1, 1); y3 = zeros(m+1, 1); y4 = zeros(m+1, 1);
    u1 = zeros(m, 1); u2 = zeros(m, 1); u3 = zeros(m, 1); u4 = zeros(m, 1);
    xi1 = zeros(m, 1); xi2 = zeros(m, 1); xi3 = zeros(m, 1); xi4 = zeros(m, 1);
    s1 = zeros(m, 1); s2 = zeros(m, 1); s3 = zeros(m, 1); s4 = zeros(m, 1);
    omega1 = zeros(m, 1); omega2 = zeros(m, 1); omega3 = zeros(m, 1); omega4 = zeros(m, 1);

    % Time vector for sinusoidal signals
    k = 1:m;
    t = (k-1) * T;

    % Compute w5 and w6 based on leader type
    if strcmp(leader_type, 'piecewise')
        w5 = zeros(1, m);
        w6 = zeros(1, m);
        w5(1:165) = 1.4; w5(166:330) = 1.6; w5(331:end) = 1.1;
        w6(1:165) = 0.6; w6(166:330) = 1.2; w6(331:end) = 0.8;
    else % sinusoidal
        w5 = eval(w5_eq);
        w6 = eval(w6_eq);
    end

    % Simulation loop
    for k = 1:m
        % Phi updates
        if k == 1
            phi1(k) = 1; phi2(k) = 1; phi3(k) = 1; phi4(k) = 1;
        elseif k == 2
            phi1(k) = phi1(k-1) + (eta * u1(k-1) / (mu + u1(k-1)^2)) * (y1(k) - phi1(k-1)*u1(k-1));
            phi2(k) = phi2(k-1) + (eta * u2(k-1) / (mu + u2(k-1)^2)) * (y2(k) - phi2(k-1)*u2(k-1));
            phi3(k) = phi3(k-1) + (eta * u3(k-1) / (mu + u3(k-1)^2)) * (y3(k) - phi3(k-1)*u3(k-1));
            phi4(k) = phi4(k-1) + (eta * u4(k-1) / (mu + u4(k-1)^2)) * (y4(k) - phi4(k-1)*u4(k-1));
        else
            phi1(k) = phi1(k-1) + (eta * (u1(k-1) - u1(k-2)) / (mu + (u1(k-1) - u1(k-2))^2)) * (y1(k) - y1(k-1) - phi1(k-1) * (u1(k-1) - u1(k-2)));
            phi2(k) = phi2(k-1) + (eta * (u2(k-1) - u2(k-2)) / (mu + (u2(k-1) - u2(k-2))^2)) * (y2(k) - y2(k-1) - phi2(k-1) * (u2(k-1) - u2(k-2)));
            phi3(k) = phi3(k-1) + (eta * (u3(k-1) - u3(k-2)) / (mu + (u3(k-1) - u3(k-2))^2)) * (y3(k) - y3(k-1) - phi3(k-1) * (u3(k-1) - u3(k-2)));
            phi4(k) = phi4(k-1) + (eta * (u4(k-1) - u4(k-2)) / (mu + (u4(k-1) - u4(k-2))^2)) * (y4(k) - y4(k-1) - phi4(k-1) * (u4(k-1) - u4(k-2)));
        end

        % Stability protection
        if k > 2 && (abs(phi1(k)) <= epsilon || abs(u1(k-1) - u1(k-2)) <= epsilon || sign(phi1(k)) ~= sign(phi1(1)))
            phi1(k) = phi1(1);
        end
        if k > 2 && (abs(phi2(k)) <= epsilon || abs(u2(k-1) - u2(k-2)) <= epsilon || sign(phi2(k)) ~= sign(phi2(1)))
            phi2(k) = phi2(1);
        end
        if k > 2 && (abs(phi3(k)) <= epsilon || abs(u3(k-1) - u3(k-2)) <= epsilon || sign(phi3(k)) ~= sign(phi3(1)))
            phi3(k) = phi3(1);
        end
        if k > 2 && (abs(phi4(k)) <= epsilon || abs(u4(k-1) - u4(k-2)) <= epsilon || sign(phi4(k)) ~= sign(phi4(1)))
            phi4(k) = phi4(1);
        end

        % Xi updates
        xi1(k) = y2(k) - 2 * y1(k) + w5(k);
        xi2(k) = y3(k) - y2(k);
        xi3(k) = y4(k) - 2 * y3(k) + y1(k);
        xi4(k) = y2(k) - 2 * y4(k) + w6(k);

        % Sliding surfaces
        if k == 1
            s1(k) = 0; s2(k) = 0; s3(k) = 0; s4(k) = 0;
        else
            s1(k) = alpha * xi1(k) - xi1(k-1);
            s2(k) = alpha * xi2(k) - xi2(k-1);
            s3(k) = alpha * xi3(k) - xi3(k-1);
            s4(k) = alpha * xi4(k) - xi4(k-1);
        end

        % Omega updates
        if k == 1
            omega1(k) = 0; omega2(k) = 0; omega3(k) = 0; omega4(k) = 0;
        else
            omega1(k) = s1(k) + tau * s1(k-1);
            omega2(k) = s2(k) + tau * s2(k-1);
            omega3(k) = s3(k) + tau * s3(k-1);
            omega4(k) = s4(k) + tau * s4(k-1);
        end

        % MFA updates
        if k == 1
            mfa1(k) = 0; mfa2(k) = 0; mfa3(k) = 0; mfa4(k) = 0;
        else
            mfa1(k) = mfa1(k-1) + (rho * phi1(k)) / (lambda + abs(phi1(k)^2)) * xi1(k);
            mfa2(k) = mfa2(k-1) + (rho * phi2(k)) / (lambda + abs(phi2(k)^2)) * xi2(k);
            mfa3(k) = mfa3(k-1) + (rho * phi3(k)) / (lambda + abs(phi3(k)^2)) * xi3(k);
            mfa4(k) = mfa4(k-1) + (rho * phi4(k)) / (lambda + abs(phi4(k)^2)) * xi4(k);
        end

        % SMC updates
        if k == 1
            sm1(k) = 0; sm2(k) = 0; sm3(k) = 0; sm4(k) = 0;
        else
            sm1(k) = sm1(k-1) + (beta * phi1(k)) / (sigma + (phi1(k)^2)) * ...
                ( (xi1(k) + (y4(k) - y4(k-1)) + (w5(k) - w5(k-1)) + (w6(k) - w6(k-1))) / 2 ...
                + (1-tau * alpha) * xi1(k) - tau * xi1(k-1) / (alpha * 2) + nena * sign(omega1(k)) );
            sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k)^2)) * ...
                ( (xi2(k) + (y3(k) - y3(k-1)) + (w5(k) - w5(k-1)) + (w6(k) - w6(k-1))) / 1 ...
                + (1-tau * alpha) * xi2(k) - tau * xi2(k-1) / (alpha * 1) + nena * sign(omega2(k)) );
            sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k)^2)) * ...
                ( (xi3(k) + (y1(k) - y1(k-1)) + (y4(k) - y4(k-1)) + (w5(k) - w5(k-1)) + (w6(k) - w6(k-1))) / 2 ...
                + (1-tau * alpha) * xi3(k) - tau * xi3(k-1) / (alpha * 2) + nena * sign(omega3(k)) );
            sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k)^2)) * ...
                ( (xi4(k) + (y2(k) - y3(k-1)) + (w5(k) - w5(k-1)) + (w6(k) - w6(k-1))) / 2 ...
                + (1-tau * alpha) * xi4(k) - tau * xi4(k-1) / (alpha * 2) + nena * sign(omega4(k)) );
        end

        % Control signals
        if k == 1
            u1(k) = 0; u2(k) = 0; u3(k) = 0; u4(k) = 0;
        else
            u1(k) = mfa1(k) + gamma1 * sm1(k);
            u2(k) = mfa2(k) + gamma2 * sm2(k);
            u3(k) = mfa3(k) + gamma3 * sm3(k);
            u4(k) = mfa4(k) + gamma4 * sm4(k);
        end

        % System dynamics
        if k == 1
            y1(k) = 0; y2(k) = 0; y3(k) = 0; y4(k) = 0;
        end
        y1(k+1) = m / (rT * 0.1) * u1(k);
        y2(k+1) = m / (rT * 0.1) * u2(k);
        y3(k+1) = m / (rT * 0.1) * u3(k);
        y4(k+1) = m / (rT * 0.3) * u4(k);
    end
% Compute and print MSE for xi_i(k)
mse_xi1 = mean((xi1).^2);
mse_xi2 = mean((xi2).^2);
mse_xi3 = mean((xi3).^2);
mse_xi4 = mean((xi4).^2);

fprintf('Lambda = %d:\n', lambda);
fprintf('  MSE of xi1 = %.6f\n', mse_xi1);
fprintf('  MSE of xi2 = %.6f\n', mse_xi2);
fprintf('  MSE of xi3 = %.6f\n', mse_xi3);
fprintf('  MSE of xi4 = %.6f\n\n', mse_xi4);


end
