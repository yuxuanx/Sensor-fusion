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
% y_samples = 0.01*samples.^3 + sqrt(2)*randn(n,1);
w_tilde = 1/n*pyx;
w = w_tilde/(sum(w_tilde));
E_xy = sum(samples.*w);
Var_xy = sum(samples.^2.*w) - E_xy^2;

%% Task c
u = sort(rand(n,1));
px = normpdf(samples,mu_x,sigma_x);
px_norm = px/sum(px);
px_cum = cumsum(px_norm);
[~,idx1] = sort([u;px_cum]);
idx2 = find(idx1<=n);
idx = idx2 - (0:n-1)';
resamples = samples(idx);
w_resample = ones(n,1)/n;
E_resample = sum(resamples.*w_resample);








