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
K = size(x_v,2);
x = zeros(4,K);
P = zeros(4,4,K);
x(:,1) = [Xk(:,1);0;0];
P(:,:,1) = eye(4);
samples = zeros(N,4,K);
samples(:,:,1) = mvnrnd(x(:,1),P(:,:,1),N,1);
w = 1/N*ones(N,1);

for i = 2:K+1
    samples(:,:,i) = mvnrnd(samples(:,:,i-1)*A',Q);
    pyx = mvnpdf(samples(:,:,i)*H',y_v(:,i-1)',R);
    w_tilde = w.*pyx;
    w = w_tilde/(sum(w_tilde));
    % Resample
    [samples(:,:,i),w] = resample(samples(:,:,i),w);
    x(:,i) = sum(w.*samples(:,:,i));
end







