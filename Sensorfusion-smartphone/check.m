
alpha = 30;
v = [0,0,1];

q = [cos(1/2*alpha);
     sin(1/2*alpha)*v(1);
     sin(1/2*alpha)*v(2);
     sin(1/2*alpha)*v(3)];
   
Qq(q)

%% Collect Measurements

[xhat, meas] = filterTemplate();

%% Visualize


nSamples = length(meas.acc);
slc = 1:nSamples;
f = figure('Position', [0,0,600,900]);
subplot(3,1,1)
plot(meas.t(slc), meas.acc(:,slc)')
xlabel('time [s]')
ylabel('acceleration [m/s^2]')
title('Accelerometer data')
legend('X', 'Y', 'Z')

subplot(3,1,2)
plot(meas.t(slc), meas.gyr(:,slc)')
xlabel('time [s]')
ylabel('rotational velocity [rad/s]')
title('Gyrometer data')
legend('X', 'Y', 'Z')

subplot(3,1,3)
plot(meas.t(slc), meas.mag(:,slc)', '.')
xlabel('time [s]')
ylabel('Tesla [T]')
title('Magnetic field data')
legend('X', 'Y', 'Z')


%% Statistics
slc = 100:nSamples - 100;

f = figure('Position', [0,0,600,900]);
subplot(3,1,1)
histogram(meas.acc(1,slc)); hold on;
histogram(meas.acc(2,slc))
histogram(meas.acc(3,slc))
xlabel('acceleration [m/s^2]')
title('Accelerometer data')
legend('X', 'Y', 'Z')

subplot(3,1,2)
histogram(meas.gyr(1,slc)); hold on;
histogram(meas.gyr(2,slc))
histogram(meas.gyr(3,slc))
xlabel('rotational velocity [rad/s]')
title('Gyrometer data')
legend('X', 'Y', 'Z')

subplot(3,1,3)
histogram(meas.mag(1,slc)); hold on;
histogram(meas.mag(2,slc))
histogram(meas.mag(3,slc))
xlabel('Tesla [T]')
title('Magnetic field data')
legend('X', 'Y', 'Z')

% noise
acc_cov = arrayfun(@(x)cov(meas.acc(x,slc), 'omitrows'), 1:3);
mag_cov = arrayfun(@(x)cov(meas.mag(x,slc), 'omitrows'), 1:3);
gyr_cov = arrayfun(@(x)cov(meas.gyr(x,slc), 'omitrows'), 1:3);

% noise
acc_mean = arrayfun(@(x)mean(meas.acc(x,slc), 'omitnan'), 1:3);
mag_mean = arrayfun(@(x)mean(meas.mag(x,slc), 'omitnan'), 1:3);
gyr_mean = arrayfun(@(x)mean(meas.gyr(x,slc), 'omitnan'), 1:3);

%% Autocorrelation

close all
direction = {'X', 'Y', 'Z'};
f = figure('Position', [0,0,700,220]);
for i = 1:3
  subplot(1,3,i)
  quant = meas.gyr(i,slc);
  autocorr(quant(~isnan(quant)));
  title(direction{i})
end
suptitle('Sample Autocorrelation Function, Angular velocity')
direction = {'X', 'Y', 'Z'};
f = figure('Position', [0,0,700,220]);
for i = 1:3
  subplot(1,3,i)
  quant = meas.acc(i,slc);
  autocorr(quant(~isnan(quant)));
  title(direction{i})
end
suptitle('Sample Autocorrelation Function, Accelerations')
direction = {'X', 'Y', 'Z'};
f = figure('Position', [0,0,700,220]);
for i = 1:3
  subplot(1,3,i)
  quant = meas.mag(i,slc);
  autocorr(quant(~isnan(quant)));
  title(direction{i})
end
suptitle('Sample Autocorrelation Function, Magnetic field')

%% Calculate covariance
meas.gyr(:,any(isnan(meas.gyr))) = [];
cov_gyr = cov(meas.gyr');
meas.acc(:,any(isnan(meas.acc))) = [];
cov_acc = cov(meas.acc');
meas.mag(:,any(isnan(meas.mag))) = [];
cov_mag = cov(meas.mag');
