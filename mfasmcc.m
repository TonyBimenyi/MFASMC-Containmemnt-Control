clc; clear;

%Parameters

rho = 7.5;
eta = 1;
lamda = 350;
mu = 0.005;
epsilon = 1e-5;
alpha = 15;
T = 0.1;
gamma1 = 0.45;
gamma2 = 0.15;
gamma3 = 0.45;
gamma4 = 0.15;

beta = 10;
sigma = 95;
tau = 1e-5;

rT = 1024;
L = 200;
m = 200;
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

% Preallocate arrays (adjust size if needed)
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


end


% Example for one time step:
xi1(k) = y2(k) - 2 * y1(k) + w5(k);
xi2(k) = y3(k) - y2(k);
xi3(k) = y4(k) - 2 * y3(k) + y1(k);
xi4(k) = y2(k) - 2 * y4(k) + w6(k);

if k == 1
    sm1(k) = 0.1;
    sm2(k) = 0.1;
    sm3(k) = 0.1;
    sm4(k) = 0.1;
else
    
end


