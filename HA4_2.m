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
% y_hat = [d;theta];

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

%% Cubature Kalman filter (CV model)
T = 0.3;
v_x = diff(Xk(1,:))/T;
v_y = diff(Xk(2,:))/T;

sigmaCV = 20; % motion noise (setting manually)
A_cv = kron([1 T;0 1],eye(2));
Q_cv = sigmaCV^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));
R = [3^2 0;0 (0.05*pi)^2];


dim = 4;
P = zeros(dim,dim,len);
x = zeros(dim,len);
P(:,:,1) = eye(dim); % initial covariance
x(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1)]; % initial state


for i = 2:len
    chi = sigmaPoints(x(:,i-1),P(:,:,i-1),dim); % generate sigma points
    x(:,i) = 1/(2*dim)*sum(A_cv*chi,2); % predicted state
    P(:,:,i) = Q_cv + 1/(2*dim)*motionCov(chi,x(:,i),dim); % predicted covariance
    % Update step
    chi = sigmaPoints(x(:,i),P(:,:,i),dim);
    y = 1/(2*dim)*sum(H_meas(chi,dim,s),2);
    Pxy = 1/(2*dim)*sum(P_xy( chi, x(:,i), y, s, dim ),3);
    S = R + 1/(2*dim)*sum(measCov( chi, dim, s, y ),3);
    x(:,i) = x(:,i) + Pxy/S*(y_hat(:,i)-y);
    P(:,:,i) = P(:,:,i) - Pxy/S*Pxy';
end

figure
plot(x(1,:),x(2,:),'*')
title('Estimation of CV model using cubature rule');
xlabel('X direction')
ylabel('Y direction')

figure
hold on
plot(x(3,:),x(4,:),'*')
plot(v_x,v_y)
title('Velocity estimation of CV model using cubature rule');
xlabel('X direction')
ylabel('Y direction')
legend('estimation','true value')

%% Cubature Kalman filter (CA model)
a_x = diff(v_x)/T;
a_y = diff(v_y)/T;

sigmaCA = 50; % motion noise (setting manually)
A_ca = kron([1 T T^2/2;0 1 T;0 0 1],eye(2));
Q_ca = sigmaCA^2*kron([T^5/20 T^4/8 T^3/6;T^4/8 T^3/3 T^2/2;...
    T^3/6 T^2/2 T],eye(2));

dim = 6;
P = zeros(dim,dim,len);
x = zeros(dim,len);
P(:,:,1) = eye(dim); % initial covariance
x(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1);a_x(1);a_y(1)]; % initial state

for i = 2:len
    chi = sigmaPoints(x(:,i-1),P(:,:,i-1),dim); % generate sigma points
    x(:,i) = 1/(2*dim)*sum(A_ca*chi,2); % predicted state
    P(:,:,i) = Q_ca + 1/(2*dim)*motionCov(chi,x(:,i),dim); % predicted covariance
    % Update step
    chi = sigmaPoints(x(:,i),P(:,:,i),dim);
    y = 1/(2*dim)*sum(H_meas(chi,dim,s),2);
    Pxy = 1/(2*dim)*sum(P_xy( chi, x(:,i), y, s, dim ),3);
    S = R + 1/(2*dim)*sum(measCov( chi, dim, s, y ),3);
    x(:,i) = x(:,i) + Pxy/S*(y_hat(:,i)-y);
    P(:,:,i) = P(:,:,i) - Pxy/S*Pxy';
end

figure
plot(x(1,:),x(2,:),'*')
title('Estimation of CA model using cubature rule');
xlabel('X direction')
ylabel('Y direction')

figure
hold on
plot(x(3,:),x(4,:),'*')
plot(v_x,v_y)
title('Velocity estimation of CV model using cubature rule');
xlabel('X direction')
ylabel('Y direction')
legend('estimation','true value')

%% 


