function [x, P] = tu_qw(x, P, omega, T, Rw)
%Time update function with angular rate measurement
A = eye(4) + 0.5*T*Somega(omega);
B = 0.5*T*Sq(x);
x = A*x;
P = A*P*A' + B*Rw*B';


end

