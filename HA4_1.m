clc;clear
close all
%% Monte Carlo integration
% Generate samples of x
mu = [2;2];
sigma = [2 -1.8;-1.8 2];
num = 10000;
x_hat = mvnrnd(mu,sigma,num);

% Generate samples of d,theta
d = sqrt(x_hat(:,1).^2+x_hat(:,2).^2)';
theta = atan(x_hat(:,2)./x_hat(:,1))';
for i = 1:num
    if theta(i) < -pi/4
        theta(i) = theta(i) + pi;
    end
end

% Calculate mean and covariance
mu_x = mean(x_hat,1);
cov_x = cov(x_hat,1);
mu_dtheta = mean([d;theta],2);
cov_dtheta = cov([d;theta]');

%% Monte Carlo illustration
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 

figure
hold on
grid on
% Samples
plot(x_hat(:,1),x_hat(:,2),'*')
% Mean
plot(mu_x(1),mu_x(2),'r+','MarkerSize',10)
% 3\Sigma ellipse
x = repmat(mu,1,n)+3*sqrtm(sigma)*[cos(phi);sin(phi)];
dtheta = repmat(mu_dtheta,1,n)+3*sqrtm(cov_dtheta)*[cos(phi);sin(phi)];
plot(x(1,:),x(2,:),'LineWidth',2)
plot(dtheta(1,:),dtheta(2,:),'LineWidth',2)
% Samples d,theta
plot(d,theta,'*')
% Mean d,theta
plot(mu_dtheta(1),mu_dtheta(2),'g+','MarkerSize',10);
xlabel('X1')
ylabel('X2')
title('Samples, expected means and 3-\sigma level curves for x and [d;\theta] using moment matching')
legend('Samples x','Mean value x','3-\sigma level curve x',...
    '3-\sigma level curve [d;\theta]','Samples [d;\theta]','Mean value [d;\theta]')

%% Sigma point integration and illustration
% Gauss-Hermite quadrature
n = 2;
p = 3;
sqrt_P = chol(sigma)';
epsilon = zeros(2,6);
epsilon(:,1) = [0;0];
epsilon(:,2) = [sqrt(3);0];
epsilon(:,3) = [-sqrt(3);0];
epsilon(:,4) = [0;sqrt(3)];
epsilon(:,5) = [sqrt(3);sqrt(3)];
epsilon(:,6) = [-sqrt(3);sqrt(3)];
epsilon(:,7) = [0;-sqrt(3)];
epsilon(:,8) = [sqrt(3);-sqrt(3)];
epsilon(:,9) = [-sqrt(3);-sqrt(3)];
sigmaPoints = mu + sqrt_P*epsilon;
d_hat = sqrt(sigmaPoints(1,:).^2+sigmaPoints(2,:).^2)';
theta_hat = atan(sigmaPoints(2,:)./sigmaPoints(1,:))';
for i = 1:p^n
    if theta_hat(i) < -pi/4
        theta_hat(i) = theta_hat(i) + pi;
    end
end

w1 = [2/3;1/6;1/6];
w1 = [2/3*w1;1/6*w1;1/6*w1];
dt1 = [d_hat theta_hat];
mu_Gauss = sum(w1.*dt1);
covGauss = zeros(2,2,p^n);
for i = 1:p^n
    covGauss(:,:,i) = w1(i).*(dt1(i,:) - mu_Gauss)'*(dt1(i,:) - mu_Gauss);
end
cov_Gauss = sum(covGauss,3);
Gauss = repmat(mu_Gauss',1,100)+3*sqrtm(cov_Gauss)*[cos(phi);sin(phi)];

figure
hold on
grid on
plot(d,theta,'*')
plot(d_hat,theta_hat,'ro','MarkerSize',10)
plot(Gauss(1,:),Gauss(2,:),'LineWidth',2)
plot(dtheta(1,:),dtheta(2,:),'LineWidth',2)
xlabel('d')
ylabel('\theta')
title('Samples and \sigma points plot of [d;\theta] using Gauss-Hermite quadrature')
legend('Samples','\sigma points','3-\sigma level curve [d;\theta] Gauss-Hermite',...
    '3-\sigma level curve [d;\theta] Monte Carlo')

% Cubature rule
n = 2;
chi = zeros(2,2*n);
chi(:,1) = mu + sqrt(2)*sqrt_P(:,1);
chi(:,2) = mu - sqrt(2)*sqrt_P(:,1);
chi(:,3) = mu + sqrt(2)*sqrt_P(:,2);
chi(:,4) = mu - sqrt(2)*sqrt_P(:,2);
d_hat2 = sqrt(chi(1,:).^2+chi(2,:).^2)';
theta_hat2 = atan(chi(2,:)./chi(1,:))';
for i = 1:2*n
    if theta_hat2(i) < -pi/4
        theta_hat2(i) = theta_hat2(i) + pi;
    end
end

dt2 = [d_hat2 theta_hat2];
mu_Cubature = 1/4*sum(dt2);
covCubature = zeros(2,2,2*n);
for i = 1:2*n
    covCubature(:,:,i) = 1/4*(dt2(i,:) - mu_Cubature)'*(dt2(i,:) - mu_Cubature);
end
cov_Cubature = sum(covCubature,3);
Cubature = repmat(mu_Cubature',1,100)+3*sqrtm(cov_Cubature)*[cos(phi);sin(phi)];

figure
hold on
grid on
plot(d,theta,'*')
plot(d_hat2,theta_hat2,'ro','MarkerSize',10)
plot(Cubature(1,:),Cubature(2,:),'LineWidth',2)
plot(dtheta(1,:),dtheta(2,:),'LineWidth',2)
xlabel('d')
ylabel('\theta')
title('Samples and \sigma points plot of [d;\theta] using cubature rule')
legend('Samples','\sigma points','3-\sigma level curve [d;\theta] Cubature',...
    '3-\sigma level curve [d;\theta] Monte Carlo')

% EKF linearization
Jaco = [mu(1)/sqrt(mu(1)^2+mu(2)^2) mu(2)/sqrt(mu(1)^2+mu(2)^2);...
    mu(1)/(mu(1)^2+mu(2)^2) -mu(2)/(mu(1)^2+mu(2)^2)];
y_hat = [sqrt(mu(1)^2+mu(2)^2);atan(mu(2)/mu(1))] + Jaco*(x_hat'-mu);
mu_y = mean(y_hat,2);
cov_y = cov(y_hat');
ekf = repmat(mu_y,1,100)+3*sqrtm(cov_y)*[cos(phi);sin(phi)];
figure
hold on
grid on
plot(d,theta,'*')
plot(ekf(1,:),ekf(2,:),'LineWidth',2)
plot(dtheta(1,:),dtheta(2,:),'LineWidth',2)
plot(mu_dtheta(1),mu_dtheta(2),'r+','MarkerSize',10);
plot(mu_y(1),mu_y(2),'g+','MarkerSize',10);
xlabel('d')
ylabel('\theta')
title('Samples and \sigma points plot of [d;\theta] using EKF linearization')
legend('Samples','3-\sigma level curve [d;\theta] EKF',...
    '3-\sigma level curve [d;\theta] Monte Carlo','Mean value of Samples',...
    'Mean value of EKF linearization')