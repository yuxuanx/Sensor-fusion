clc;clear
%% Task a
mu_x = 3;
sigma_x = 4;
x = linspace(mu_x - 3*sigma_x,mu_x + 3*sigma_x,200);
h_x = 0.01*x.^3;
y = 3.5*ones(length(x),1);
figure
xlim([mu_x - 3*sigma_x,mu_x + 3*sigma_x])
hold on
grid on
plot(x,h_x);
plot(x,y)
xlabel('x')
legend('h(x)','y')
title('h(x) and y')

p_x = normpdf(x,mu_x,sigma_x);
p_y = normpdf(h_x,3.5,sqrt(2));
figure
hold on
grid on
plot(x,p_x)
plot(x,p_y)
title('Distribution of p(x) and p(y|x)')
xlabel('x')
ylabel('p')
legend('p(x)','p(y|x)')

%% Task b
% Generate samples from proposal
n = 50;
samples = normrnd(mu_x,sigma_x,n,1);
pyx = normpdf(0.01*samples.^3,3.5,sqrt(2));
% y_samples = 0.01*samples.^3 + sqrt(2)*randn(n,1);
w_tilde = 1/n*pyx;
w = w_tilde/(sum(w_tilde));
E_xy = sum(samples.*w);
Var_xy = sum(samples.^2.*w) - E_xy^2;


%% Task c
% Resampling
u = sort(rand(n,1));
px_cum = cumsum(w);
[~,idx1] = sort([u;px_cum]);
idx2 = find(idx1<=n);
idx = idx2 - (0:n-1)';
resamples = samples(idx);
w_resample = ones(n,1)/n;
E_resample = sum(resamples.*w_resample);

p = find([n;diff(idx);n]);
instances = diff(p);
[~,idx_unique] = unique(idx);
resamples = resamples(idx_unique);
n_unique = length(instances);
w_resample = ones(n_unique,1)/n.*instances;

figure
subplot(2,1,1)
stem(samples,w)
xlabel('samples')
ylabel('p')
title('Original samples after resampling along with their weights')
subplot(2,1,2)
stem(resamples,w_resample)
xlabel('samples')
ylabel('p')
title('New samples after resampling along with their weights')







