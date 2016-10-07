clc;clear
%% Task a
mu_x = 3;
sigma_x = 4;
x = linspace(mu_x - 3*sigma_x,mu_x + 3*sigma_x,100);
h_x = 0.01*x.^3;
y = 3.5*ones(length(x),1);
figure
xlim([mu_x - 3*sigma_x,mu_x + 3*sigma_x])
hold on
plot(x,h_x);
plot(x,y)
xlabel('x')
legend('h(x)','y')

p_x = normpdf(x,mu_x,sigma_x);
p_y = normpdf(0.01*x.^3,3.5,sqrt(2));
figure
hold on
plot(x,p_x)
plot(x,p_y)
title('Distribution of p(x) and p(y|x)')
xlabel('x')
ylabel('p')
legend('p(x)','p(y|x)')

%% Task b
% Generate samples from proposal
n = 200;
samples = normrnd(mu_x,sigma_x,n,1);
pyx = normpdf(0.01*samples.^3,3.5,sqrt(2));
y_samples = 0.01*samples.^3 + sqrt(2)*randn(n,1);
w_tilde = 1/n*pyx;
weight = w_tilde/(sum(w_tilde));
E_xy = sum(y_samples.*weight);
Var_xy = sum(y_samples.^2.*weight) - E_xy^2;











