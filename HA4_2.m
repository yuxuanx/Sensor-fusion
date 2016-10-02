clc;clear
close all
%% Generate Samples
load('Xk.mat')
s = [50;-150];
len = size(Xk,2);
d = zeros(1,len);
theta = zeros(1,len);
for i = 1:len
    d(i) = sqrt((Xk(1,i)-s(1))^2+(Xk(2,i)-s(2))^2);
    theta(i) = atan2(Xk(2,i)-s(2),Xk(1,i)-s(1));
    %         if theta(i) < 0
    %             theta(i) = theta(i) + pi;
    %         end
end

y_hat = [d + 3*randn(size(d));theta + 0.015*pi*randn(size(theta))];
%y_hat = [d;theta];

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

sigmaCV = 2; % motion noise (setting manually)

% dist = zeros(96,1);
% j = 1;
% for sigmaCV = 0.5:0.1:10

A_cv = kron([1 T;0 1],eye(2));
Q_cv = sigmaCV^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));
R = [3^2 0;0 (0.015*pi)^2];
%R = zeros(2,2);

dim = 4;
P = zeros(dim,dim,len);
x = zeros(dim,len);
P(:,:,1) = eye(dim); % initial covariance
x(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1)]; % initial state


for i = 2:len
    chi = sigmaPoints(x(:,i-1),P(:,:,i-1),dim); % generate sigma points
    x(:,i) = 1/(2*dim)*sum(A_cv*chi,2); % predicted state
    P(:,:,i) = Q_cv + 1/(2*dim)*motionCov(A_cv,chi,x(:,i),dim); % predicted covariance
    % Update step
    chi = sigmaPoints(x(:,i),P(:,:,i),dim);
    measurements = H_meas(chi,dim,s);
    y = 1/(2*dim)*sum(measurements,2);
    Pxy = 1/(2*dim)*sum(P_xy(measurements, chi, x(:,i), y, dim ),3);
    S = R + 1/(2*dim)*sum(measCov( measurements, dim, y ),3);
    x(:,i) = x(:,i) + Pxy/S*(y_hat(:,i)-y);
    P(:,:,i) = P(:,:,i) - Pxy/S*Pxy';
end

% for idx = 1:len-1
%     dist(j) = dist(j) + (x(:,idx)-[Xk(:,idx);v_x(idx);v_y(idx)])'*(x(:,idx)-[Xk(:,idx);v_x(idx);v_y(idx)]);
% end
% j = j + 1;
% end

figure
hold on
plot(x(1,:),x(2,:),'*')
plot(Xk(1,:),Xk(2,:),'*')
title('Estimation of CV model using cubature rule');
legend('estimation','true trajectory')
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

sigmaCA = 0.6; % motion noise (setting manually)
% dist = zeros(96,1);
% j = 1;
% for sigmaCA = 0.5:0.1:10

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
    P(:,:,i) = Q_ca + 1/(2*dim)*motionCov(A_ca,chi,x(:,i),dim); % predicted covariance
    % Update step
    chi = sigmaPoints(x(:,i),P(:,:,i),dim);
    measurements = H_meas(chi,dim,s);
    y = 1/(2*dim)*sum(measurements,2);
    Pxy = 1/(2*dim)*sum(P_xy( measurements, chi, x(:,i), y, dim ),3);
    S = R + 1/(2*dim)*sum(measCov( measurements, dim, y ),3);
    x(:,i) = x(:,i) + Pxy/S*(y_hat(:,i)-y);
    P(:,:,i) = P(:,:,i) - Pxy/S*Pxy';
end

% for idx = 1:len-2
%     dist(j) = dist(j) + (x(1:4,idx)-[Xk(:,idx);v_x(idx);v_y(idx)])'*(x(1:4,idx)-[Xk(:,idx);v_x(idx);v_y(idx)]);
% end
% j = j + 1;
% end

figure
plot(x(1,:),x(2,:),'*')
title('Estimation of CA model using cubature rule');
xlabel('X direction')
ylabel('Y direction')

figure
hold on
plot(x(3,:),x(4,:),'*')
plot(v_x,v_y)
title('Velocity estimation of CA model using cubature rule');
xlabel('X direction')
ylabel('Y direction')
legend('estimation','true value')

%% Cubature RTS smoother (CV model)
sigmaCV = 2; % motion noise (setting manually)
A_cv = kron([1 T;0 1],eye(2));
Q_cv = sigmaCV^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));
R = [3^2 0;0 (0.015*pi)^2];
%R = zeros(2,2);

dim = 4;
P_predict = zeros(dim,dim,len-1);
P_update = zeros(dim,dim,len);
x_predict = zeros(dim,len-1);
x_update = zeros(dim,len);
P_update(:,:,1) = eye(dim); % initial covariance
x_update(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1)]; % initial state

% Forward filtering
for i = 2:len
    chi = sigmaPoints(x_update(:,i-1),P_update(:,:,i-1),dim); % generate sigma points
    x_predict(:,i-1) = 1/(2*dim)*sum(A_cv*chi,2); % predicted state
    P_predict(:,:,i-1) = Q_cv + 1/(2*dim)*motionCov(A_cv,chi,x_predict(:,i-1),dim); % predicted covariance
    % Update step
    chi = sigmaPoints(x_predict(:,i-1),P_predict(:,:,i-1),dim);
    measurements = H_meas(chi,dim,s);
    y = 1/(2*dim)*sum(measurements,2);
    Pxy = 1/(2*dim)*sum(P_xy(measurements,chi,x_predict(:,i-1),y,dim),3);
    S = R + 1/(2*dim)*sum(measCov(measurements,dim,y),3);
    x_update(:,i) = x_predict(:,i-1) + Pxy/S*(y_hat(:,i)-y);
    P_update(:,:,i) = P_predict(:,:,i-1) - Pxy/S*Pxy';
end


x_smooth = zeros(dim,len);
P_smooth = zeros(dim,dim,len);
x_smooth(:,end) = x_update(:,end);
P_smooth(:,:,end) = P_update(:,:,end);
G = zeros(dim,dim,len-1);
% Backward recursion
for i = 1:len-1
    G(:,:,i) = P_update(:,:,end-i)*A_cv'/P_predict(:,:,end-i+1);
    x_smooth(:,end-i) = x_update(:,end-i) + G(:,:,i)*...
        (x_smooth(:,end-i+1) - x_predict(:,end-i+1));
    P_smooth(:,:,end-i) = P_update(:,:,end-i) + G(:,:,i)*...
        (P_smooth(:,:,end-i+1) - P_predict(:,:,end-i+1))*G(:,:,i)';
end


figure
hold on
plot(x_smooth(1,:),x_smooth(2,:),'*')
title('Estimation of CV model using RTS smoother');
xlabel('X direction')
ylabel('Y direction')

%% No measurements between time step 150~175
sigmaCV = 2; % motion noise (setting manually)
A_cv = kron([1 T;0 1],eye(2));
Q_cv = sigmaCV^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));
R = [3^2 0;0 (0.015*pi)^2];
%R = zeros(2,2);

dim = 4;
P_predict = zeros(dim,dim,len-1);
P_update = zeros(dim,dim,len);
x_predict = zeros(dim,len-1);
x_update = zeros(dim,len);
P_update(:,:,1) = eye(dim); % initial covariance
x_update(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1)]; % initial state

% Forward filtering
for i = 2:len
    chi = sigmaPoints(x_update(:,i-1),P_update(:,:,i-1),dim); % generate sigma points
    x_predict(:,i-1) = 1/(2*dim)*sum(A_cv*chi,2); % predicted state
    P_predict(:,:,i-1) = Q_cv + 1/(2*dim)*motionCov(A_cv,chi,x_predict(:,i-1),dim); % predicted covariance
    % Update step
    chi = sigmaPoints(x_predict(:,i-1),P_predict(:,:,i-1),dim);
    measurements = H_meas(chi,dim,s);
    y = 1/(2*dim)*sum(measurements,2);
    Pxy = 1/(2*dim)*sum(P_xy(measurements,chi,x_predict(:,i-1),y,dim),3);
    S = R + 1/(2*dim)*sum(measCov(measurements,dim,y),3);
    if i >= 150 && i <= 175
        x_update(:,i) = x_predict(:,i-1);
        P_update(:,:,i) = P_predict(:,:,i-1);
    else
        x_update(:,i) = x_predict(:,i-1) + Pxy/S*(y_hat(:,i)-y);
        P_update(:,:,i) = P_predict(:,:,i-1) - Pxy/S*Pxy';
    end
end


x_smooth = zeros(dim,len);
P_smooth = zeros(dim,dim,len);
x_smooth(:,end) = x_update(:,end);
P_smooth(:,:,end) = P_update(:,:,end);
G = zeros(dim,dim,len-1);
% Backward recursion
for i = 1:len-1
    G(:,:,i) = P_update(:,:,end-i)*A_cv'/P_predict(:,:,end-i+1);
    x_smooth(:,end-i) = x_update(:,end-i) + G(:,:,i)*...
        (x_smooth(:,end-i+1) - x_predict(:,end-i+1));
    P_smooth(:,:,end-i) = P_update(:,:,end-i) + G(:,:,i)*...
        (P_smooth(:,:,end-i+1) - P_predict(:,:,end-i+1))*G(:,:,i)';
end


figure
hold on
plot(x_update(1,:),x_update(2,:),'*')
plot(x_smooth(1,:),x_smooth(2,:),'*')
title('Estimation of CV model using RTS smoother');
xlabel('X direction')
ylabel('Y direction')
legend('Kalman filtering','RTS smoother')

%% MSE for individual state elements
trial = 100;
se_px = zeros(trial,len);
se_py = zeros(trial,len);
se_vx = zeros(trial,len-1);
se_vy = zeros(trial,len-1);
covPosterior = zeros(dim,len,trial);
for t = 1:trial
    y_hat = [d + 3*randn(size(d));theta + 0.015*pi*randn(size(theta))];
    dim = 4;
    P = zeros(dim,dim,len);
    x = zeros(dim,len);
    P(:,:,1) = eye(dim); % initial covariance
    x(:,1) = [Xk(1,1);Xk(2,1);v_x(1);v_y(1)]; % initial state
    
    
    for i = 2:len
        chi = sigmaPoints(x(:,i-1),P(:,:,i-1),dim); % generate sigma points
        x(:,i) = 1/(2*dim)*sum(A_cv*chi,2); % predicted state
        P(:,:,i) = Q_cv + 1/(2*dim)*motionCov(A_cv,chi,x(:,i),dim); % predicted covariance
        % Update step
        chi = sigmaPoints(x(:,i),P(:,:,i),dim);
        measurements = H_meas(chi,dim,s);
        y = 1/(2*dim)*sum(measurements,2);
        Pxy = 1/(2*dim)*sum(P_xy(measurements, chi, x(:,i), y, dim ),3);
        S = R + 1/(2*dim)*sum(measCov( measurements, dim, y ),3);
        x(:,i) = x(:,i) + Pxy/S*(y_hat(:,i)-y);
        P(:,:,i) = P(:,:,i) - Pxy/S*Pxy';
    end
    
    for i = 1:len
        covPosterior(:,i,t) = diag(P(:,:,i));
    end
    se_px(t,:) = (x(1,:) - Xk(1,:)).^2;
    se_py(t,:) = (x(2,:) - Xk(2,:)).^2;
    se_vx(t,:) = (x(3,1:end-1) - v_x).^2;
    se_vy(t,:) = (x(4,1:end-1) - v_y).^2;
    
end
mse_px = mean(se_px);
mse_py = mean(se_py);
mse_vx = mean(se_vx);
mse_vy = mean(se_vy);

averageCov = mean(covPosterior,3);

figure
hold on
plot(averageCov');
xlabel('time step')
ylabel('posterior covariance')
title('Average Posterior covariance of individual state over 100 simulations')
legend('position x','position y','velocity x','velocity y')

figure
hold on
plot(mse_px)
plot(mse_py)
plot(mse_vx)
plot(mse_vy)
xlabel('time step')
ylabel('mean square error')
title('Mean square error of individual state over 100 simulations')
legend('position x','position y','velocity x','velocity y')





