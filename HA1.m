%% Problem 2.1
% Sample code to plot 3?sigma ellipsoids
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 
mu = [1;1]; % Mean of a 2D normal random variable
P = [1 0;0 1]; % Covariance of a 2?D normal random variable

x_hat = mvnrnd(mu, P, n)';

figure;plot(x_hat(1,:),x_hat(2,:),'*');hold on

% 3\Sigma ellipse
x = repmat(mu,1,n)+3*sqrtm(P)*[cos(phi);sin(phi)]; % Equation(12) in HA1 document 
plot(x(1,:),x(2,:),'-g','LineWidth',2)
xlabel('X1')
ylabel('X2')
title('3-\sigma Level Curves for a 2D Normal random variable')
legend('Samples','3-\sigma level curve')

%% Problem 2.2
A1 = [1 0;0 2];
A2 = [0.1 1.8;0 2];
b1 = [0;0];
b2 = [5;5];
 
mx = [1;1];
Qx = [1 0;0 1];
 
my1 = A1*mx + b1
my2 = A2*mx + b2
 
Qy1 = A1 * Qx * A1'
Qy2 = A2 * Qx * A2'
 
y1_hat = A1*x_hat+repmat(b1,1,n);
y2_hat = A2*x_hat+repmat(b2,1,n);
 
y1 = repmat(my1,1,n)+3*sqrtm(Qy1)*[cos(phi);sin(phi)]; % Equation(12) in HA1 document 
y2 = repmat(my2,1,n)+3*sqrtm(Qy2)*[cos(phi);sin(phi)]; % Equation(12) in HA1 document 

figure;plot(y1_hat(1,:),y1_hat(2,:),'*');hold on
plot(y2_hat(1,:),y2_hat(2,:),'*');
plot(y1(1,:),y1(2,:),'-','LineWidth',2);
plot(y2(1,:),y2(2,:),'-','LineWidth',2);
xlabel('Y1')
ylabel('Y2')
title('3-\sigma Level Curves for a 2D Normal random variable')
legend('Samples-y1','Samples-y2','3-\sigma level curve-y1','3-\sigma level curve-y2')


%% Problem 2.3
mu1 = mean(y1_hat,2)
mu2 = mean(y2_hat,2)
cov1 = cov(y1_hat(1,:),y1_hat(2,:))
cov2 = cov(y2_hat(1,:),y2_hat(2,:))

%% Problem 4.1
n = 100; % Number of grid points
phi = linspace(0,2*pi,100); % Create grid in interval [0,2*pi] 
mu = [19;19];
P = [1;1]*25*[1 1] + [0;1]*4*[0 1];
% 3\Sigma ellipse
xy = repmat(mu,1,n)+3*sqrtm(P)*[cos(phi);sin(phi)]; % Equation(12) in HA1 document 
figure
plot(xy(1,:),xy(2,:),'-g','LineWidth',2)
xlabel('x')
ylabel('y')
title('3-\sigma Level Curves for a 2D Normal random variable')
legend('3-\sigma level curve')

muy1 = 17.3;
muy2 = 25.9;
sigma = sqrt(100/29);
y = linspace(0,100,1000);
x1 = normpdf(y,muy1,sigma);
x2 = normpdf(y,muy2,sigma);
figure
plot(y,x1);hold on
plot(y,x2);
legend('p(x|y_1)','p(x|y_2)');title('Posterior distribution')
