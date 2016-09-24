clc;clear;clf
% Load data
run('StateSimulation.m')
% Generate measurements
s1 = [-1.5 0.5];
s2 = [1 1];
measurements = zeros(2,K);
for i = 1:K
measurements(:,i) = H_measurement( X(:,i),s1,s2 ) + 0.1^2*[randn;randn];
end

% Constant velocity
T = 0.01;
A_vel = [eye(2) T*eye(2);zeros(2,2) eye(2)];


% Constant acceleration
A_acc = [eye(2) T*eye(2) T^2/2*eye(2);...
    zeros(2,2) eye(2) T*eye(2);...
    zeros(2,2) zeros(2,2) eye(2)];


% Prediction
Q = [zeros(2,4);zeros(2,2) eye(2)]*T*10^2;
Q = [zeros(2,6);zeros(2,6);zeros(2,4) eye(2)]*T*100^2;
P0 = 0.01*eye(6);
x0 = [0 0 1 0 0 0]';
P = zeros(6,6,K);
x_hat = zeros(6,K);
P(:,:,1) = P0;
x_hat(:,1) = x0;
S = zeros(2,2,K);
KK = zeros(6,2,K);
R = 0.1^2*eye(2);
h_prime = zeros(2,6,K);

% Used for plotting
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 
xy_predict = zeros(2,n,K);
xy_posterior = zeros(2,n,K);

for i = 2:K
    % Predict
    x_hat(:,i) = A_acc*x_hat(:,i-1);
    P(:,:,i) = A_acc*P(:,:,i-1)*A_acc' + Q;
    xy_predict(:,:,i) = repmat(x_hat(1:2,i),1,n)+3*sqrtm(P(1:2,1:2,i))*[cos(phi);sin(phi)];
    % Update
    h_prime(:,:,i) = [H_derivative(x_hat(:,i),s1,1) H_derivative(x_hat(:,i),s1,0) 0 0 0 0;...
        H_derivative(x_hat(:,i),s2,1) H_derivative(x_hat(:,i),s2,0) 0 0 0 0];
    S(:,:,i) = h_prime(:,:,i)*P(:,:,i)*h_prime(:,:,i)' + R;
    KK(:,:,i) = P(:,:,i)*h_prime(:,:,i)'/S(:,:,i);
    x_hat(:,i) = x_hat(:,i) + KK(:,:,i)*(measurements(:,i) - H_measurement( x_hat(:,i),s1,s2 ));
    P(:,:,i) = P(:,:,i) - KK(:,:,i)*S(:,:,i)*KK(:,:,i)';
    xy_posterior(:,:,i) = repmat(x_hat(1:2,i),1,n)+3*sqrtm(P(1:2,1:2,i))*[cos(phi);sin(phi)];
end




hold on
plot(x_hat(1,:),x_hat(2,:),'-.','LineWidth',2)
plot(s1(1),s1(2),'*r')
plot(s2(1),s2(2),'*r')
title('True and estimated trajectories')
xlabel('X-Position')
ylabel('Y-Position')
legend('Real','Estimation')
figure
hold on
plot(s1(1),s1(2),'*')
plot(s2(1),s2(2),'*')
plot(xy_predict(1,:,50),xy_predict(2,:,50),'LineWidth',2)
plot(xy_posterior(1,:,50),xy_posterior(2,:,50),'LineWidth',2)
plot(xy_predict(1,:,100),xy_predict(2,:,100),'LineWidth',2)
plot(xy_posterior(1,:,100),xy_posterior(2,:,100),'LineWidth',2)
plot(xy_predict(1,:,250),xy_predict(2,:,250),'LineWidth',2)
plot(xy_posterior(1,:,250),xy_posterior(2,:,250),'LineWidth',2)
plot(X(1,50),X(2,50),'x')
plot(X(1,100),X(2,100),'x')
plot(X(1,250),X(2,250),'x')





