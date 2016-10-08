clc;clear
%% Kalman filter
% Parameters
Q = 1.5;
R = 2.5;
% The first time step is x_0
T = 20+1; % time step
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
    P(i) = A*P(i-1)*A' + Q;
    
    S = H*P(i)*H' + R;
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
N = 100;
x_sis = zeros(T,1);
P_hat = zeros(T,1);
x_sis(1) = 2;
P_hat(1) = 6;
q = sqrt(Q);
r = sqrt(R);
w = zeros(N,T);
w(:,1) = 1/N*ones(N,1);
samples = zeros(N,T);
for i = 2:T
    samples(:,i-1) = normrnd(samples(:,i-1),q);
    pyx = normpdf(samples(:,i-1),y(i),r);
    w_tilde = w(:,i-1).*pyx;
    w(:,i) = w_tilde/sum(w_tilde);
    x_sis(i) = sum(w(:,i).*samples(:,i-1));
    P_hat(i) = mean((x_sis(i)-samples(:,i-1)).^2);
end

mse_pf = mean((x-x_sis).^2);

figure
subplot(2,1,1)
hold on
plot(x_sis)
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

figure
hold on
plot(samples(:,1:end-1)','k-*')
plot(x_sis(2:end),'Linewidth',2)

%PostPlot(samples(:,1:end-1),w(:,2:end),x_hat(2:end)',P(2:end)',N,T-1,0)

%% Partical filter with resampling
x_sir = zeros(T,1);
x_sir(1) = 2;
for i = 2:T
    samples(:,i-1) = normrnd(samples(:,i-1),q);
    pyx = normpdf(samples(:,i-1),y(i),r);
    w_tilde = w(:,i-1).*pyx;
    w(:,i) = w_tilde/sum(w_tilde);
    [samples(:,i-1),w(:,i)] = resample(samples(:,i-1),w(:,i));
    x_sir(i) = sum(w(:,i).*samples(:,i-1));
    P_hat(i) = mean((x_sir(i)-samples(:,i-1)).^2);
end

mse_pfre = mean((x-x_sir).^2);

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

%PostPlot(samples(:,1:end-1),w(:,2:end),x_hat(2:end)',P(2:end)',N,T-1,1)
