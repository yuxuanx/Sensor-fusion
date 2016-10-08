%% Generate velocity measurements
load('Xk.mat')
x_v = diff(Xk')';
sigma_r = 0.05;
R = sigma_r^2*eye(2);
y_v = x_v + sigma_r*randn(size(x_v));

figure
plot(Xk(1,:),Xk(2,:))
figure
plot(x_v')

%% Particle filter
T = 1;
A = kron([1 T;0 1],eye(2));
sigma_q = 0.1;
Q = sigma_q^2*kron([T^3/3 T^2/2;T^2/2 T],eye(2));
H = [0 0 1 0;0 0 0 1];

N = 100;


