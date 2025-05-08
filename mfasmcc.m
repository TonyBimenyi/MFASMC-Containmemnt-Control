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





