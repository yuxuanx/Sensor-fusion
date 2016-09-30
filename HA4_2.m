clc;clear
close all
%% Generate Samples
load('Xk.mat')
s = [300;0];
len = size(Xk,2);
d = zeros(1,len);
theta = zeros(1,len);
for i = 1:len
    d(i) = sqrt((Xk(1,i)-s(1))^2+(Xk(2,i)-s(2))^2);
    theta(i) = atan2(Xk(2,i)-s(2),Xk(1,i)-s(1));
    if theta(i) < 0
        theta(i) = theta(i) + 2*pi;
    end
end

y_hat = [d + 3*randn(1,200);theta + 0.05*pi*randn(1,200)];

figure
plot(Xk(1,:),Xk(2,:),'*')
title('Trajectory samples');
xlabel('X direction');
ylabel('Y direction')

figure
plot(y_hat(1,:),y_hat(2,:),'.')
title('Sensor measurement')
xlabel('Distance');
ylabel('Angle (radians)');

%% Cubature Kalman filter (CA model)
v_x = diff(Xk(1,:));
v_y = diff(Xk(2,:));

sigmaCV = 0.1; % motion noise (setting manually)
T = 0.3;
A_cv = kron([1 T;0 1],eye(2));
Q_cv = sigmaCV^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));

dim = 4;
P = zeros(dim,dim,len);
sqrt_P = zeros(dim,dim,len); % storing P^(1/2);
x = zeros(dim,len);
P(:,:,1) = eye(dim); % initial covariance
x(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1)]; % initial state


for i = 2:len
    chi = sigmaPoints(x(:,i-1),P(:,:,i-1),dim); % generate sigma points
    x(:,i) = 1/(2*dim)*sum(A_cv*chi,2); % predicted state
    P(:,:,i) = Q_cv + 1/(2*dim)*motionCov(chi,x(:,i),dim); % predicted covariance
    
end







