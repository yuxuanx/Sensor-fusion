function [x, P] = tu_q(x, P, T, Rw)
%Time update function without angular rate measurement
A = eye(4);
x = A*x;
B = 0.5*T*Sq(x);
P = A*P*A' + B*Rw*B';

end

