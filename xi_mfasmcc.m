clc; clear;

%Parameters

rho = 10.5;
eta = 5;
lamda = 100;
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
n =600;


% initialization 

phi1 = zeros(m+1, 1);
phi2 = zeros(m+1, 1);
phi3 = zeros(m+1, 1);
phi4 = zeros(m+1, 1);

mfa1 = zeros(m+1, 1);
mfa2 = zeros(m+1, 1);
mfa3 = zeros(m+1, 1);
mfa4 = zeros(m+1, 1);

sm1 = zeros(m,1);
sm2 = zeros(m,1);
sm3 = zeros(m,1);
sm4 = zeros(m,1);

y1 = zeros(m+1, 1);
y2 = zeros(m+1, 1);
y3 = zeros(m+1, 1);
y4 = zeros(m+1, 1);

u1 = zeros(m,1);
u2 = zeros(m,1);
u3 = zeros(m,1);
u4 = zeros(m,1);

xi1 = zeros(m, 1);
xi2 = zeros(m, 1);
xi3 = zeros(m, 1);
xi4 = zeros(m, 1);

s1 = zeros(m, 1);
s2 = zeros(m, 1);
s3 = zeros(m, 1);
s4 = zeros(m, 1);

omega1 = zeros(m, 1);
omega2 = zeros(m, 1);
omega3 = zeros(m, 1);
omega4 = zeros(m, 1);

%Time Invarying Leaders signals (w5, w6)
% % Preallocate arrays (adjust size if needed)
w5 = zeros(1, m);  % example length, adjust as needed
w6 = zeros(1, m);  % example length, adjust as needed

% Set w5 values according to conditions
w5(1:165)     = 1.4;
w5(166:330)   = 1.6;
w5(331:end)   = 1.3;

% Set w6 values according to conditions
w6(1:165)     = 0.7;
w6(166:330)   = 1.2;
w6(331:end)   = 1.1;


% % Time vector for sinusoidal signals
% k = 1:m;
% t = (k-1) * T; % Time vector: t = (k-1)*T

% % Define sinusoidal leader signals (choose one of the variations below)

% % Variation 1: Low-frequency, same amplitude, phase-shifted signals
% w5 = 0.5 + 0.25 * sin(0.1 * t); % Amplitude = 0.25, frequency = 0.1 rad/s, offset = 1.15
% w6 = 0.5 + 0.25 * sin(0.1 * t + pi/2); % Same amplitude and frequency, phase shift = pi/2


for k = 1:m

    if k == 1
        phi1(k) = 1;
        phi2(k) = 1;
        phi3(k) = 1;
        phi4(k) = 1;
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
    if k > 2 && (abs(phi1(k)) <= epsilon || abs(u1(k - 1) - u1(k - 2)) <= epsilon || sign(phi1(k)) ~= sign(phi1(1)))
        phi1(k) = phi1(1);
    end
    
    if k > 2 && (abs(phi2(k)) <= epsilon || abs(u2(k - 1) - u2(k - 2)) <= epsilon || sign(phi2(k)) ~= sign(phi2(1)))
        phi2(k) = phi2(1);
    end
    
    if k > 2 && (abs(phi3(k)) <= epsilon || abs(u3(k - 1) - u3(k - 2)) <= epsilon || sign(phi3(k)) ~= sign(phi3(1)))
        phi3(k) = phi3(1);
    end
    
    if k > 2 && (abs(phi4(k)) <= epsilon || abs(u4(k - 1) - u4(k - 2)) <= epsilon || sign(phi4(k)) ~= sign(phi4(1)))
        phi4(k) = phi4(1);
    end

    % Example for one time step:
    xi1(k) = y2(k) - 2 * y1(k) + w5(k);
    xi2(k) = y3(k) - y2(k);
    xi3(k) = y4(k) - 2 * y3(k) + y1(k);
    xi4(k) = y2(k) - 2 * y4(k) + w6(k);

    % Fix: Handle k=1 case for sliding surfaces
    if k == 1
        s1(k) = 0;
        s2(k) = 0;
        s3(k) = 0;
        s4(k) = 0;
    else
        s1(k) = alpha * xi1(k) - xi1(k-1);
        s2(k) = alpha * xi2(k) - xi2(k-1);
        s3(k) = alpha * xi3(k) - xi3(k-1);
        s4(k) = alpha * xi4(k) - xi4(k-1);
    end

    % Fix: Handle k=1 case for sliding surfaces
    if k == 1
        omega1(k) = 0;
        omega2(k) = 0;
        omega3(k) = 0;
        omega4(k) = 0;
    else
        omega1(k) = s1(k) + tau * (s1(k-1));
        omega2(k) = s2(k) + tau * (s2(k-1));
        omega3(k) = s3(k) + tau * (s3(k-1));
        omega4(k) = s4(k) + tau * (s4(k-1));
    end

    if k == 1
        mfa1(k) = 0;
        mfa2(k) = 0;
        mfa3(k) = 0;
        mfa4(k) = 0;
    else
        mfa1(k) = mfa1(k-1) + (rho * phi1(k)) / (lamda + abs(phi1(k)^2)) * xi1(k);
        mfa2(k) = mfa2(k-1) + (rho * phi2(k)) / (lamda + abs(phi2(k)^2)) * xi2(k);
        mfa3(k) = mfa3(k-1) + (rho * phi3(k)) / (lamda + abs(phi3(k)^2)) * xi3(k);
        mfa4(k) = mfa4(k-1) + (rho * phi4(k)) / (lamda + abs(phi4(k)^2)) * xi4(k);
    end


    % SMC updates
    if k == 1 
        sm1(k) = 0;
        sm2(k) = 0;
        sm3(k) = 0;
        sm4(k) = 0;
    else
        sm1(k) = sm1(k-1) + (beta * phi1(k)) / (sigma + (phi1(k))^2) * ...
            ( (xi1(k) + (y4(k) - y4(k-1)) + (w5(k) - w5(k-1))+(w6(k) - w6(k-1))) / (1 + 1) ...
            + (1-tau * alpha) * xi1(k) - tau * xi1(k-1) / (alpha * (2)) + nena * sign(omega1(k)) );

        sm2(k) = sm2(k-1) + (beta * phi2(k)) / (sigma + (phi2(k))^2) * ...
            ( (xi2(k) + (y3(k) - y3(k-1)) + (w5(k) - w5(k-1))+(w6(k) - w6(k-1))) / (1) ...
            + (1-tau * alpha) * xi2(k) - tau * xi2(k-1) / (alpha * (1)) + nena * sign(omega2(k)) );


        sm3(k) = sm3(k-1) + (beta * phi3(k)) / (sigma + (phi3(k))^2) * ...
            ( (xi3(k) + (y1(k) - y1(k-1)) + (y4(k) - y4(k-1)) + (w5(k) - w5(k-1))+(w6(k) - w6(k-1))) / (2) ...
            + (1-tau * alpha) * xi3(k) - tau * xi3(k-1) / (alpha * (2)) + nena * sign(omega3(k)) );


        sm4(k) = sm4(k-1) + (beta * phi4(k)) / (sigma + (phi4(k))^2) * ...
            ( (xi4(k) + (y2(k) - y3(k-1)) + (w5(k) - w5(k-1))+(w6(k) - w6(k-1))) / (2) ...
            + (1-tau * alpha) * xi4(k) - tau * xi4(k-1) / (alpha * (2)) + nena * sign(omega4(k)) );
    end

    % Control signal
    if k == 1
        u1(k) = 0;
        u2(k) = 0;
        u3(k) = 0;
        u4(k) = 0;
    else
        u1(k) = mfa1(k) + gamma1 * sm1(k);
        u2(k) = mfa2(k) + gamma2 * sm2(k);
        u3(k) = mfa3(k) + gamma3 * sm3(k);
        u4(k) = mfa4(k) + gamma4 * sm4(k);
    end
    if k == 1
        y1(k) = 0;
        y2(k) = 0;
        y3(k) = 0;
        y4(k) = 0;
    end

    % System dynamics (example, replace with actual system equations)

    y1(k + 1) = m / (rT * 0.1) * u1(k);
    y2(k + 1) = m / (rT * 0.1) * u2(k);
    y3(k + 1) = m / (rT * 0.1) * u3(k);
    y4(k + 1) = m / (rT * 0.3) * u4(k);

    
    % Ensure y values do not exceed bounds (example, adjust as needed)
 
    % y1(k) = max(0, min(1, y1(k)));
    % y2(k) = max(0, min(1, y2(k)));
    % y3(k) = max(0, min(1, y3(k)));
    % y4(k) = max(0, min(1, y4(k)));
    % % Ensure u values do not exceed bounds (example, adjust as needed)
    % u1(k) = max(0, min(1, u1(k)));
    % u2(k) = max(0, min(1, u2(k)));
    % u3(k) = max(0, min(1, u3(k)));
    % u4(k) = max(0, min(1, u4(k)));
    % % Ensure phi values do not exceed bounds (example, adjust as needed)
    % phi1(k) = max(0, min(1, phi1(k)));
    % phi2(k) = max(0, min(1, phi2(k)));
    % phi3(k) = max(0, min(1, phi3(k)));
    % phi4(k) = max(0, min(1, phi4(k)));

    % Plotting



end

% Time vector for plotting (assuming k represents discrete time steps)
k = 1:m;

% Plot system outputs (y1, y2, y3, y4) and leaders (w5, w6)
figure;
plot(k, xi1(1:m), 'b-', 'LineWidth', 2.5, 'DisplayName', 'y1');
hold on;
plot(k, xi2(1:m), 'r-', 'LineWidth', 2.5, 'DisplayName', 'y2');
plot(k, xi3(1:m), 'g-', 'LineWidth', 2.5, 'DisplayName', 'y3');
plot(k, xi4(1:m), 'm-', 'LineWidth', 2.5, 'DisplayName', 'y4');
% plot(k, w5, 'c--', 'LineWidth', 2.5, 'DisplayName', 'w5');
% plot(k, w6, 'k--', 'LineWidth', 2.5, 'DisplayName', 'w6');
hold off;
title('System Outputs and Leader Signals');
xlabel('Iteration (k)');
ylabel('Value');
legend('show');
grid on;



    


