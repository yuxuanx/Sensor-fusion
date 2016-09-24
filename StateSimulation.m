clear all
T  = 0.01;          % Sampling period
x0 = [0;0;1;0];     % Initial state
K = 500;            % Length of time sequence
a = zeros(1,K);
a(1,50:100)  = pi/2/51/T + 0.01*randn(1,51);
a(1,200:250) = pi/2/51/T + 0.01*randn(1,51);
a(1,350:400) = pi/2/51/T + 0.01*randn(1,51);
x = x0;
X = [];
for i=1:K
  F = [0 0  1    0;...
       0 0  0    1;...
       0 0  0   a(i);...
       0 0 -a(i) 0];
  x = expm(F*T)*x;
  X = [X x];
end
close all; figure(1);clf;plot(X(1,:),X(2,:),'.','LineWidth',2);
title('Trajectory of the moving target')
xlabel('X position')
ylabel('Y position')
%axis([-1 1 -1.8 0.2])

