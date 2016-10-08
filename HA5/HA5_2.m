clc;clear
%% Kalman filter
% Parameters
Q = 1.5;
R = 2.5;
T = 1000; % time step
A = 1;
H = 1;

% Initialization
P = zeros(T,1);
x = zeros(T,1);
x_hat = zeros(T,1);
y = zeros(T,1);

% Initial state
x(1) = 2;
x_hat(1) = 2;
P(1) = 6;

% Generate samples and measurements
for i = 2:T
    x(i) = x(i-1) + sqrt(Q)*randn;
end

for i = 1:T
    y(i) = x(i) + sqrt(R)*randn;
end


% Kalman filtering
for i = 2:T
    x_hat(i) = A*x_hat(i-1);
    P(i) = A*P(i-1)*A + Q;
    
    S = H*P(i)*H + R;
    K = P(i)*H/S;
    v = y(i) - H*x_hat(i);
    
    x_hat(i) = x_hat(i) + K*v;
    P(i) = P(i) - K*S*K';
end

mse_kf = mean((x-x_hat).^2);

figure
subplot(2,1,1)
hold on
plot(x_hat)
plot(x)
xlabel('time step k')
ylabel('mean value')
legend('estimation','true value')
title('Expected state over time')
subplot(2,1,2)
plot(P)
xlabel('time step k')
ylabel('variance')
title('Variance of state over time')

%% Partical filter without resampling
% Draw samples from proposal (motion model)
N = 500;
x_hat = zeros(T,1);
P_hat = zeros(T,1);
x_hat(1) = 2;
P_hat(1) = 6;
q = sqrt(Q);
r = sqrt(R);
w = 1/N*ones(N,1);
for i = 2:T
    samples = normrnd(x_hat(i-1),q,N,1);
    pyx = normpdf(samples,y(i),r);
    w_tilde = w.*pyx;
    w = w_tilde/sum(w_tilde);
    x_hat(i) = sum(w.*samples);
    P_hat(i) = mean((x_hat(i)-samples).^2);
end

mse_pf = mean((x-x_hat).^2);

figure
subplot(2,1,1)
hold on
plot(x_hat)
plot(x)
xlabel('time step k')
ylabel('mean value')
legend('estimation','true value')
title('Expected state over time')
subplot(2,1,2)
plot(P_hat)
xlabel('time step k')
ylabel('variance')
title('Variance of state over time')

%% Partical filter with resampling
for i = 2:T
    samples = normrnd(x_hat(i-1),q,N,1);
    pyx = normpdf(samples,y(i),r);
    w_tilde = w.*pyx;
    w = w_tilde/sum(w_tilde);
    [samples,w] = resample(samples,w);
    x_hat(i) = sum(w.*samples);
    P_hat(i) = mean((x_hat(i)-samples).^2);
end

mse_pfre = mean((x-x_hat).^2);

figure
subplot(2,1,1)
hold on
plot(x_hat)
plot(x)
xlabel('time step k')
ylabel('mean value')
legend('estimation','true value')
title('Expected state over time')
subplot(2,1,2)
plot(P_hat)
xlabel('time step k')
ylabel('variance')
title('Variance of state over time')

