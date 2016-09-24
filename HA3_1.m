%% Subtask b
clc;clear;clf
mu_x0 = [1;2];
P_x0 = [1 0;0 1]*0.002;
T = 2;
A = [1 T;0 1];
Qc = 0.1;
Q_euler = [0 0;0 1]*Qc*T;
Q_exact = [T^3/2 T^2/2;T^2/2 T]*Qc;

mu_x1 = A*mu_x0;
% Prediction using modified Euler
P_euler = A*P_x0*A' + Q_euler;
% Prediction using exact solution
P_exact = A*P_x0*A' + Q_exact;

% 3\Sigma ellipse
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi]
x_euler = repmat(mu_x1,1,n)+3*sqrtm(P_euler)*[cos(phi);sin(phi)];
x_exact = repmat(mu_x1,1,n)+3*sqrtm(P_exact)*[cos(phi);sin(phi)];
figure(1)
hold on
plot(x_euler(1,:),x_euler(2,:),'-g','LineWidth',2)
plot(x_exact(1,:),x_exact(2,:),'-r','LineWidth',2)
xlabel('position')
ylabel('velocity')
title('3-\sigma Level Curves for position & velocity')
legend('Euler','Exact solution')

%% Subtask c
K = 50;
K = 100;
T = 2;
A = [1 T;0 1];
Q_euler = [0 0;0 1]*Qc*T;
Q_exact = [T^3/2 T^2/2;T^2/2 T]*Qc;
P_euler = zeros(2,2,K);
P_euler(:,:,1) = P_x0;
P_exact = zeros(2,2,K);
P_exact(:,:,1) = P_x0;
x = zeros(2,K);
x(:,1) = mu_x0;
for i = 1:K-1
    x(:,i+1) = A*x(:,i);
    P_euler(:,:,i+1) = A*P_euler(:,:,i)*A' + Q_euler;
    P_exact(:,:,i+1) = A*P_exact(:,:,i)*A' + Q_exact;
end
x_euler = repmat(x(:,K),1,n)+3*sqrtm(P_euler(:,:,K))*[cos(phi);sin(phi)];
x_exact = repmat(x(:,K),1,n)+3*sqrtm(P_exact(:,:,K))*[cos(phi);sin(phi)];
figure(2)
hold on
plot(x_euler(1,:),x_euler(2,:),'-g','LineWidth',2)
plot(x_exact(1,:),x_exact(2,:),'-r','LineWidth',2)
xlabel('position')
ylabel('velocity')
title('3-\sigma Level Curves for position & velocity')
legend('Euler','Exact solution')

