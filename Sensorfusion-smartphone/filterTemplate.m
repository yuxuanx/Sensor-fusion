function [xhat, meas] = filterTemplate(calAcc, calGyr, calMag)
% FILTERTEMPLATE  Filter template
%
% This is a template function for how to collect and filter data
% sent from a smartphone live.  Calibration data for the
% accelerometer, gyroscope and magnetometer assumed available as
% structs with fields m (mean) and R (variance).
%
% The function returns xhat as an array of structs comprising t
% (timestamp), x (state), and P (state covariance) for each
% timestamp, and meas an array of structs comprising t (timestamp),
% acc (accelerometer measurements), gyr (gyroscope measurements),
% mag (magnetometer measurements), and orint (orientation quaternions
% from the phone).  Measurements not availabe are marked with NaNs.
%
% As you implement your own orientation estimate, it will be
% visualized in a simple illustration.  If the orientation estimate
% is checked in the Sensor Fusion app, it will be displayed in a
% separate view.
%
% Note that it is not necessary to provide inputs (calAcc, calGyr, calMag).

%% Setup necessary infrastructure
import('com.liu.sensordata.*');  % Used to receive data.

%% Filter settings
t0 = [];  % Initial time (initialize on first data received)
nx = 4;   % Assuming that you use q as state variable.
% Add your filter settings here.
f = 100; % Sampling frequency
alpha = 0.05; % Used in AR filter to rstimate magnetic field
g0 = [-0.7103;-0.0161;9.7996]; % Nominal gravity vector
% Should be estimated from training data
m0 = [0;sqrt(10.8348^2+19.1457^2);-38.0160]; % Earth magnetic field
L = norm(m0);
% Process noise (Gyrometer)
Rw = 10e-7*[35.171 -7.6 -1.76;...
    -7.6 59.475 3.07;...
    -1.76 3.07 47.350];
% Measurement noise (Accelerator)
Ra = 10e-7*[158.74 2.1232 9.3710;...
    2.1232 140.39 3.1530;...
    9.3710 3.1530 655.80];
% Should be estimated from training data
% Measurement noise (Magnetometer)
Rm = [0.1567 -0.0123 0.0366;...
    -0.0123 0.1264 0.0035;...
    0.0366 0.0035 1.2410];

% Current filter state.
x = [1; 0; 0 ;0];
P = eye(nx, nx);

% Saved filter states.
xhat = struct('t', zeros(1, 0),...
    'x', zeros(nx, 0),...
    'P', zeros(nx, nx, 0));

meas = struct('t', zeros(1, 0),...
    'acc', zeros(3, 0),...
    'gyr', zeros(3, 0),...
    'mag', zeros(3, 0),...
    'orient', zeros(4, 0));
try
    %% Create data link
    server = StreamSensorDataReader(3400);
    % Makes sure to resources are returned.
    sentinel = onCleanup(@() server.stop());
    
    server.start();  % Start data reception.
    
    % Used for visualization.
    figure(1);
    subplot(1, 2, 1);
    ownView = OrientationView('Own filter', gca);  % Used for visualization.
    googleView = [];
    counter = 0;  % Used to throttle the displayed frame rate.
    
    %% Filter loop
    while server.status()  % Repeat while data is available
        % Get the next measurement set, assume all measurements
        % within the next 5 ms are concurrent (suitable for sampling
        % in 100Hz).
        data = server.getNext(5);
        
        if isnan(data(1))  % No new data received
            continue;        % Skips the rest of the look
        end
        t = data(1)/1000;  % Extract current time
        
        if isempty(t0)  % Initialize t0
            t0 = t;
        end
        
        acc = data(1, 2:4)';
        
        if ~any(isnan(acc))  % Acc measurements are available.
            if norm(acc) > 9.81
                accOut = 1;
            else
                accOut = 0;
                [x, P] = mu_g(x, P, acc, Ra, g0);
                [x, P] = mu_normalizeQ(x, P);
            end
        end
        gyr = data(1, 5:7)';
        if ~any(isnan(gyr))  % Gyro measurements are available.
            [x, P] = tu_qw(x, P, gyr, 1/f, Rw);
            [x, P] = mu_normalizeQ(x, P);
        end
        
        mag = data(1, 8:10)';
        if ~any(isnan(mag))  % Mag measurements are available.
            if abs(L-norm(mag)) > 10
                magOut = 1;
            else
                magOut = 0;
            [x, P] = mu_m(x, P, mag, m0, Rm);
            [x, P] = mu_normalizeQ(x, P);
            L = (1-alpha)*L + alpha*norm(mag);
            end
        end
        
        orientation = data(1, 18:21)';  % Google's orientation estimate.
        
        % Visualize result
        if rem(counter, 10) == 0
            ownView.setAccDist(accOut);
            setOrientation(ownView, x(1:4));
            title(ownView, 'OWN', 'FontSize', 16);
            if ~any(isnan(orientation))
                if isempty(googleView)
                    subplot(1, 2, 2);
                    % Used for visualization.
                    googleView = OrientationView('Google filter', gca);
                end
                googleView.setMagDist(magOut);
                setOrientation(googleView, orientation);
                title(googleView, 'GOOGLE', 'FontSize', 16);
            end
        end
        counter = counter + 1;
        
        % Save estimates
        xhat.x(:, end+1) = x;
        xhat.P(:, :, end+1) = P;
        xhat.t(end+1) = t - t0;
        
        meas.t(end+1) = t - t0;
        meas.acc(:, end+1) = acc;
        meas.gyr(:, end+1) = gyr;
        meas.mag(:, end+1) = mag;
        meas.orient(:, end+1) = orientation;
    end
catch e
    fprintf(['Unsuccessful connecting to client!\n' ...
        'Make sure to start streaming from the phone *after*'...
        'running this function!']);
end
end
