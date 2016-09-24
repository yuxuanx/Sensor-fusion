clc;clear;clf
% Load data
run('StateSimulation.m')
% Generate measurements
s1 = [-1.5 0.5];
s2 = [1 1];
measurements = zeros(2,K);
for i = 1:K
measurements(:,i) = H_measurement( X(:,i),s1,s2 ) + 0.1*[randn;randn];
end

% Constant velocity
T = 0.01;
A_vel = [eye(2) T*eye(2);zeros(2,2) eye(2)];


% Constant acceleration
A_acc = [eye(2) T*eye(2) T^2/2*eye(2);...
    zeros(2,2) eye(2) T*eye(2);...
    zeros(2,2) zeros(2,2) eye(2)];
% 
% err = zeros(length(0.05:0.001:0.5),1);
% idx = 1;
% for Q_c = 0.05:0.001:0.5
% Prediction
Q_c = 0.0145;
Q_c = 0.4580;
Q = [zeros(2,4);zeros(2,2) eye(2)]*T*Q_c;
%Q = [zeros(2,6);zeros(2,6);zeros(2,4) eye(2)]*T*100^2;
P0 = 0.01*eye(4);
x0 = [0 0 1 0]';
P = zeros(4,4,K);
x_hat = zeros(4,K);
P(:,:,1) = P0;
x_hat(:,1) = x0;
S = zeros(2,2,K);
KK = zeros(4,2,K);
R = 0.1^2*eye(2);
h_prime = zeros(2,4,K);
v = zeros(2,K);

% Used for plotting
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 
xy_predict = zeros(2,n,K);
xy_posterior = zeros(2,n,K);
x_predict = zeros(4,K);

for i = 2:K
    % Predict
    x_hat(:,i) = A_vel*x_hat(:,i-1);
    P(:,:,i) = A_vel*P(:,:,i-1)*A_vel' + Q;
    x_predict(:,i) = x_hat(:,i);
    xy_predict(:,:,i) = repmat(x_hat(1:2,i),1,n)+3*sqrtm(P(1:2,1:2,i))*[cos(phi);sin(phi)];
    % Update
    h_prime(:,:,i) = [H_derivative(x_hat(:,i),s1,1) H_derivative(x_hat(:,i),s1,0) 0 0;...
        H_derivative(x_hat(:,i),s2,1) H_derivative(x_hat(:,i),s2,0) 0 0];
    S(:,:,i) = h_prime(:,:,i)*P(:,:,i)*h_prime(:,:,i)' + R;
    KK(:,:,i) = P(:,:,i)*h_prime(:,:,i)'/S(:,:,i);
    v(:,i) = measurements(:,i) - H_measurement( x_hat(:,i),s1,s2 );
    x_hat(:,i) = x_hat(:,i) + KK(:,:,i)*v(:,i);
    P(:,:,i) = P(:,:,i) - KK(:,:,i)*S(:,:,i)*KK(:,:,i)';
    xy_posterior(:,:,i) = repmat(x_hat(1:2,i),1,n)+3*sqrtm(P(1:2,1:2,i))*[cos(phi);sin(phi)];
end

% for i = 1:K
%     err(idx) = err(idx) + (X(:,i)-x_hat(:,i))'*(X(:,i)-x_hat(:,i));
% end
% idx = idx + 1;
% % err = (X - x_hat)*(X - x_hat)^T;
% % 
% % lms(idx) = sum(err(:));
% % idx = idx + 1;
% end

x1 = linspace(-1.5,1,100);
y1_80 = s1(2)-(x1-s1(1))*tan(-measurements(1,80));
y1_240 = s1(2)-(x1-s1(1))*tan(-measurements(1,240));
y1_400 = s1(2)-(x1-s1(1))*tan(-measurements(1,400));
y2_80 = s2(2) + (s2(1)-x1)*-tan(measurements(2,80));
y2_240 = s2(2) + (s2(1)-x1)*-tan(measurements(2,240));
y2_400 = s2(2) + (s2(1)-x1)*-tan(measurements(2,400));



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
plot(xy_predict(1,:,80),xy_predict(2,:,80),'LineWidth',2)
plot(xy_posterior(1,:,80),xy_posterior(2,:,80),'LineWidth',2)
plot(xy_predict(1,:,240),xy_predict(2,:,240),'LineWidth',2)
plot(xy_posterior(1,:,240),xy_posterior(2,:,240),'LineWidth',2)
plot(xy_predict(1,:,400),xy_predict(2,:,400),'LineWidth',2)
plot(xy_posterior(1,:,400),xy_posterior(2,:,400),'LineWidth',2)
plot(s1(1),s1(2),'^')
plot(s2(1),s2(2),'^')
plot(X(1,80),X(2,80),'*')
plot(X(1,240),X(2,240),'*')
plot(X(1,400),X(2,400),'*')
plot(x1,y1_80,'--','LineWidth',1.5);
plot(x1,y1_240,'--','LineWidth',1.5);
plot(x1,y1_400,'--','LineWidth',1.5);
plot(x1,y2_80,'--','LineWidth',1.5);
plot(x1,y2_240,'--','LineWidth',1.5);
plot(x1,y2_400,'--','LineWidth',1.5);
xlabel('X-position')
ylabel('Y-position')
legend('predicted density k = 80','posterior density k = 80',...
    'predicted density k = 240','posterior density k = 240',...
    'predicted density k = 400','posterior density k = 400',...
    'sensor 1','sensor 2','true position k = 80','true position k = 240',...
    'true position k = 400','s1 measurement k = 80','s1 measurement k = 240',...
    's1 measurement k = 400','s2 measurement k = 80','s2 measurement k = 240',...
    's2 measurement k = 400')


% cov(X(1,:),x_hat(1,:))
% cov(X(2,:),x_hat(2,:))
% cov(X(3,:),x_hat(3,:))
% cov(X(4,:),x_hat(4,:))
err = (X - x_hat).^2;
figure
hold on
title('Estimation error square over 500 time steps')
xlabel('time step k')
ylabel('logarithm')
plot(err(1,:),'LineWidth',2)
plot(err(2,:),'LineWidth',2)
plot(err(3,:),'LineWidth',2)
plot(err(4,:),'LineWidth',2)
legend('estimation error square p_x','estimation error square p_y',...
    'estimation error square v_x','estimation error square v_y')
figure
hold on
xlabel('time step k')
ylabel('logarithm')
title('Posterior variance of position and velocity over 500 time steps')
plot(squeeze(P(1,1,:)),'LineWidth',2)
plot(squeeze(P(2,2,:)),'LineWidth',2)
plot(squeeze(P(3,3,:)),'LineWidth',2)
plot(squeeze(P(4,4,:)),'LineWidth',2)
legend('posterior variance p_x','posterior variance p_y'...
    ,'posterior variance v_x','posterior variance v_y')

figure
hold on
xlabel('X-position')
ylabel('Y-position')
plot(s1(1),s1(2),'^')
plot(s2(1),s2(2),'^')
plot(x1,y1_80,'--','LineWidth',1.5);
plot(x1,y1_240,'--','LineWidth',1.5);
plot(x1,y1_400,'--','LineWidth',1.5);
plot(x1,y2_80,'--','LineWidth',1.5);
plot(x1,y2_240,'--','LineWidth',1.5);
plot(x1,y2_400,'--','LineWidth',1.5);
quiver(x_predict(1,80),x_predict(2,80),KK(1,1,80).*v(1,80),KK(2,1,80).*v(1,80),50,'LineWidth',2)
quiver(x_predict(1,80),x_predict(2,80),KK(1,2,80).*v(2,80),KK(2,2,80).*v(2,80),50,'LineWidth',2)
quiver(x_predict(1,240),x_predict(2,240),KK(1,1,240).*v(1,240),KK(2,1,240).*v(1,240),50,'LineWidth',2)
quiver(x_predict(1,240),x_predict(2,240),KK(1,2,240).*v(2,240),KK(2,2,240).*v(2,240),50,'LineWidth',2)
quiver(x_predict(1,400),x_predict(2,400),KK(1,1,400).*v(1,400),KK(2,1,400).*v(1,400),50,'LineWidth',2)
quiver(x_predict(1,400),x_predict(2,400),KK(1,2,400).*v(2,400),KK(2,2,400).*v(2,400),50,'LineWidth',2)
quiver(x_predict(1,80),x_predict(2,80),v(1,80),v(1,80),50,'LineWidth',2)
% quiver(x_predict(1,400),x_predict(2,400),KK(1,1,400).*v(1,400),KK(2,1,400).*v(1,400))
legend('sensor 1','sensor 2','s1 measurement k = 80','s1 measurement k = 240',...
    's1 measurement k = 400','s2 measurement k = 80',...
    's2 measurement k = 240','s2 measurement k = 400',...
    'Kalman gain * innnovation s1 k = 80','Kalman gain * innnovation s2 k = 80',...
    'Kalman gain * innnovation s1 k = 240','Kalman gain * innnovation s2 k = 240',...
    'Kalman gain * innnovation s1 k = 400','Kalman gain * innnovation s2 k = 400',...
    'mm')