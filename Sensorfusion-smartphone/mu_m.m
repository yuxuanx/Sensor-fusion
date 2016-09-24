function [x, P] = mu_m(x, P, mag, m0, Rm)
%Magnetometer measurement update 
[dQ1,dQ2,dQ3,dQ4] = dQqdq(x);
J = [dQ1'*m0 dQ2'*m0 dQ3'*m0 dQ4'*m0];
S = J*P*J' + Rm;
Vs= chol(S); inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
K = P*J'*iS;
x = x + K*(mag - Qq(x)'*m0);
P = P - K*S*K';


end

