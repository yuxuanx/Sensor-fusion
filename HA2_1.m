clc;clear
Q = 1.5;
R = 2.5;
T = 50;

A = 1;
H = 1;
P = zeros(T,1);
K = zeros(T,1);
S = zeros(T,1);
v = zeros(T,1);
x = zeros(T,1);
x(1) = 0;
P(1) = 1; % initial value

figure(1);hold on
xc = linspace(-20,20,200);

for i = 2:T
    x(i) = A*x(i-1);
    P(i) = A*P(i-1)*A + Q;
    if mod(i,10) == 0
        px = normpdf(xc,x(i),sqrt(P(i)));
        plot(xc,px);
    end
end

xlabel('x_k');ylabel('p(x_k)');
title('Probability density distribution of x_k')
legend('k = 10','k = 20','k = 30','k = 40','k = 50')

