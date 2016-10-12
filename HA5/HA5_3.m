clc;clear
close all
%% Generate velocity measurements
load('Xk.mat')
x_v = diff(Xk')';
sigma_r = 0.05;
R = sigma_r^2*eye(2);
y_v = x_v + sigma_r*randn(size(x_v));

% figure
% plot(Xk(1,:),Xk(2,:))
% figure
% plot(x_v')

%% Particle filter
T = 1;
A = kron([1 T;0 1],eye(2));
sigma_q = 0.1;
Q = sigma_q^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));
H = [0 0 1 0;0 0 0 1];

N = 1000;
K = size(x_v,2);
x = zeros(4,K);
P = zeros(4,4,K);
x(:,1) = [Xk(:,1);x_v(1,1);x_v(2,1)];
P(:,:,1) = 10*eye(4);
samples = zeros(N,4,K);
samples(:,:,1) = mvnrnd(x(:,1),P(:,:,1),N,1);
% xx = 1 + 10*rand(N,1);
% yy = 1 + 8*rand(N,1);
% samples(:,:,1) = [xx yy zeros(N,2)];

w = 1/N*ones(N,1);

for i = 2:K+1
    samples(:,:,i) = mvnrnd(samples(:,:,i-1)*A',Q);
    temp = samples(:,:,i);
    idx = isOnRoad(temp(:,1),temp(:,2));
    samplesOnRoad = temp(idx==1,:);
    w = w(idx==1);
    pyx = mvnpdf(samplesOnRoad*H',y_v(:,i-1)',R);
    w_tilde = w.*pyx;
    w = w_tilde/(sum(w_tilde));
    % Resample
    [samples(:,:,i),w] = resample2(samplesOnRoad,w,N);
    x(:,i) = sum(w.*samples(:,:,i));
end

% for i = 2:K+1
%     samples(:,:,i) = mvnrnd(samples(:,:,i-1)*A',Q);
%     temp = samples(:,:,i);
%     idx = isOnRoad(temp(:,1),temp(:,2));
%     w(idx==0) = 0;
%     pyx = mvnpdf(temp*H',y_v(:,i-1)',R);
%     w_tilde = w.*pyx;
%     w = w_tilde/(sum(w_tilde));
%     % Resample
%     [samples(:,:,i),w] = resample2(temp,w,N);
%     x(:,i) = sum(w.*samples(:,:,i));
% end
% figure
% plot(x(1,:),x(2,:),'*')

figure
clf
hold on
axis([0.8 11.2 0.8 9.2])
h1 = plot(Xk(1,:),Xk(2,:));
h2 = plot(x(1,2:end),x(2,2:end),'-*');
plot([1+1i 1+9*1i 5+9*1i])
plot([7+9*1i 11+9*1i 11+1i 7+1i]);plot([5+1i 1+1i])
plot([2+5.2*1i 2+8.3*1i 4+8.3*1i 4+5.2*1i 2+5.2*1i])%House 1
plot([2+3.7*1i 2+4.4*1i 4+4.4*1i 4+3.7*1i 2+3.7*1i])%House 2
plot([2+2*1i 2+3.2*1i 4+3.2*1i 4+2*1i 2+2*1i])%House 3
plot([5+1i 5+2.2*1i 7+2.2*1i 7+1i])%House 4
plot([5+2.8*1i 5+5.5*1i 7+5.5*1i 7+2.8*1i 5+2.8*1i])%House 5
plot([5+6.2*1i 5+9*1i]);plot([7+9*1i 7+6.2*1i 5+6.2*1i])%House 6
plot([8+4.6*1i 8+8.4*1i 10+8.4*1i 10+4.6*1i 8+4.6*1i])%House 7
plot([8+2.4*1i 8+4*1i 10+4*1i 10+2.4*1i 8+2.4*1i])%House 8
plot([8+1.7*1i 8+1.8*1i 10+1.8*1i 10+1.7*1i 8+1.7*1i])%House 9
title('Particle filter estimation along with true trajectory')
xlabel('X position')
ylabel('Y position')
legend([h1,h2],'true trajectory','estimated trajectory')
for i = 2:K+1
    h = plot(samples(:,1,i),samples(:,2,i),'r.');
    pause
    delete(h)
    drawnow;
end

figure
subplot(2,1,1)
hold on
plot(x_v(1,:));
plot(x(3,2:end));
xlabel('time step')
ylabel('velocity')
title('X velocity')
legend('true velocity','estimation')
subplot(2,1,2)
hold on
plot(x_v(2,:));
plot(x(4,2:end));
xlabel('time step')
ylabel('velocity')
title('Y velocity')
legend('true velocity','estimation')





