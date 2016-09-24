clc;clear
Q = 1.5;
R = 2.5;
T = 10000; % time step
A = 1;
H = 1;

% initialization
P = zeros(T,1);
K = zeros(T,1);
S = zeros(T,1);
v = zeros(T,1);
x = zeros(T,1);
x_hat = zeros(T,1);
y = zeros(T,1);

% initial position
x(1) = 0;
x_hat(1) = 0;
P(1) = 0; % initial covariance

% generate samples and measurements
for i = 2:T
    x(i) = x(i-1) + sqrt(Q)*randn;
end

for i = 1:T
    y(i) = x(i) + sqrt(R)*randn;
end

N = 200;
xc = linspace(-10,10,N);
posterior = zeros(N,6);
predict = zeros(N,3);
j = 1;
k = 1;

% Kalman filtering  
for i = 2:T
    x_hat(i) = A*x_hat(i-1);
    P(i) = A*P(i-1)*A + Q;
    if (i == 5)||(i == 20)||(i == 50)
        predict(:,k) = normpdf(xc,x_hat(i),sqrt(P(i)));
        k = k + 1;
    end
    
    S(i) = H*P(i)*H + R;
    K(i) = P(i)*H/S(i);
    v(i) = y(i) - H*x_hat(i);
    
    x_hat(i) = x_hat(i) + K(i)*v(i);
    P(i) = P(i) - K(i)*S(i)*K(i);
    if i==4||5||19||20||49||50
        posterior(:,j) = normpdf(xc,x_hat(i),sqrt(P(i)));
        j = j + 1;
    end
end

figure
plot(x,'--','LineWidth',1.5);hold on
plot(x_hat,'LineWidth',1.5);
plot(y,'-.','LineWidth',1.5)
title('Data generated, observation and filter estimates over 100 steps')
xlabel('time step k');ylabel('position')
legend('real data','estimation','observation')

figure
hold on
plot(xc,predict(:,1),'LineWidth',1.5)
plot(xc,predict(:,2),'LineWidth',1.5)
plot(xc,predict(:,3),'LineWidth',1.5)
plot(xc,posterior(:,2),'--','LineWidth',1.5)
plot(xc,posterior(:,4),'--','LineWidth',1.5)
plot(xc,posterior(:,6),'--','LineWidth',1.5)
plot(xc,posterior(:,1),'-.','LineWidth',1.5)
plot(xc,posterior(:,3),'-.','LineWidth',1.5)
plot(xc,posterior(:,5),'-.','LineWidth',1.5)
plot(y(5),0,'*')
plot(y(20),0,'*')
plot(y(50),0,'*')
title('Predicted and posterior distribution at different time steps')
xlabel('x_k')
ylabel('p(x_k)')
legend('predicted (k=5)','predicted (k=20)','predicted (k=50)'...
    ,'posterior (k=5)','posterior (k=20)','posterior (k=50)'...
    ,'posterior (k=4)','posterior (k=19)','posterior (k=49)'...
    ,'measurement y_5','measurement y_{20}','measurement y_{50}')

figure
p1 = normpdf(xc,x(1),sqrt(Q));
p2 = normpdf(xc,x_hat(2),sqrt(P(2)));
hold on
plot(xc,p1,'LineWidth',1.5)
plot(xc,p2,'LineWidth',1.5)
plot(y(2),0,'*r')
plot([x(1) x(1)],[0,max(p1)],'--r')
plot([x_hat(2) x_hat(2)],[0,max(p2)],'--r')
xlabel('position')
ylabel('p')
title('Probability density distribution')
legend('p(x_1)','p(x_1|y_1)','position y_1');

a = mean(x(20:end)-x_hat(20:end))
b = 3*sqrt(P(T)/(T-20))
c = cov(x_hat(20:end)-x(20:end))
d = P(20)