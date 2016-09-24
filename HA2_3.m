clc;clear
% parameters
Q = 0.01;
R = 2.5;
T = 100; % time step
A = [1 1;0 1];
H = [1 0];
q = [0 0;0 Q];

% initialization
P = zeros(2,2,T);
K = zeros(2,T);
S = zeros(T,1);
v = zeros(T,1);
x = zeros(2,T);
x_hat = zeros(2,T);
y = zeros(T,1);

% initial position
x(:,1) = [0;0];
x_hat(:,1) = [0;0];
P(:,:,1) = zeros(2,2); % initial covariance

for i = 2:T
    x(:,i) = A*x(:,i-1) + [0;1]*sqrt(0.5)*randn;
end

for i = 1:T
    y(i) = x(1,i) + sqrt(R)*randn;
end

N = 200;
xc = linspace(-50,50,N);
posterior = zeros(N,6);
predict = zeros(N,3);
k = 1;
j = 1;

% Kalman filtering
for i = 2:T
    x_hat(:,i) = A*x_hat(:,i-1);
    P(:,:,i) = A*P(:,:,i-1)*A' + q;
    if (i == 5)||(i == 10)||(i == 20)
        predict(:,k) = normpdf(xc,x_hat(1,i),sqrt(P(1,1,i)));
        k = k + 1;
    end
    
    S(i) = H*P(:,:,i)*H' + R;
    K(:,i) = P(:,:,i)*H'/S(i);
    v(i) = y(i) - H*x_hat(:,i);
    
    x_hat(:,i) = x_hat(:,i) + K(:,i)*v(i);
    P(:,:,i) = P(:,:,i) - K(:,i)*S(i)*K(:,i)';
    if (i==4)||(i==5)||(i==9)||(i==10)||(i==19)||(i==20)
        posterior(:,j) = normpdf(xc,x_hat(1,i),sqrt(P(1,1,i)));
        j = j + 1;
    end
end

% a = mean(x(1,10:end)-x_hat(1,10:end))
% b = mean(x(2,10:end)-x_hat(2,10:end))
% c = cov(x_hat(1,10:end)-x(1,10:end))
% d = cov(x_hat(2,10:end)-x(2,10:end))

% figure(1)
% subplot(2,1,1)
% plot(x(1,:),'-.','LineWidth',1.5);hold on
% plot(x_hat(1,:),'LineWidth',1.5);
% plot(y,'--','LineWidth',1.5)
% title('Data generated, observation and filter estimation')
% xlabel('time step k')
% ylabel('position')
% legend('real position','estimated position','observation')
% 
% subplot(2,1,2)
% plot(x(2,:),'-.','LineWidth',1.5);hold on
% plot(x_hat(2,:),'LineWidth',1.5);
% title('Data generated, observation and filter estimation')
% xlabel('time step k')
% ylabel('velocity')
% legend('real velocity','estimated velocity')
% 
% figure(2)
% clf
% hold on
% plot(xc,predict(:,1),'LineWidth',1.5)
% plot(xc,predict(:,2),'LineWidth',1.5)
% plot(xc,predict(:,3),'LineWidth',1.5)
% plot(xc,posterior(:,2),'--','LineWidth',1.5)
% plot(xc,posterior(:,4),'--','LineWidth',1.5)
% plot(xc,posterior(:,6),'--','LineWidth',1.5)
% plot(xc,posterior(:,1),'-.','LineWidth',1.5)
% plot(xc,posterior(:,3),'-.','LineWidth',1.5)
% plot(xc,posterior(:,5),'-.','LineWidth',1.5)
% plot(y(5),0,'*')
% plot(y(10),0,'*')
% plot(y(20),0,'*')
% title('Predicted and posterior distribution at different time steps')
% xlabel('x_k')
% ylabel('p(x_k)')
% legend('predicted (k=5)','predicted (k=10)','predicted (k=20)'...
%     ,'posterior (k=5)','posterior (k=10)','posterior (k=20)'...
%     ,'posterior (k=4)','posterior (k=9)','posterior (k=19)'...
%     ,'measurement y_{5}','measurement y_{10}','measurement y_{20}')
% 
% 
figure(3)
plot(xcov(v)/max(xcov(v)));
title('Autocorrelation of innovation')
xlabel('lag')
ylabel('sample autocorrelation')


figure
subplot(2,2,1)
plot(a/max(a))
xlabel('lag')
ylabel('sample autocorrelation')
title('Q = 0.05')
subplot(2,2,2)
plot(b/max(b))
xlabel('lag')
ylabel('sample autocorrelation')
title('Q = 0.5')
subplot(2,2,3)
plot(c/max(c))
xlabel('lag')
ylabel('sample autocorrelation')
title('Q = 5')
subplot(2,2,4)
plot(d/max(d))
xlabel('lag')
ylabel('sample autocorrelation')
title('Q = 0.8')