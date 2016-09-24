function [x, P] = mu_g(x, P, yacc, Ra, g0)
%Accelerometer measurement update
[dQ1,dQ2,dQ3,dQ4] = dQqdq(x);
J = [dQ1'*g0 dQ2'*g0 dQ3'*g0 dQ4'*g0];
S = J*P*J' + Ra;
Vs= chol(S); inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';
K = P*J'*iS;
x = x + K*(yacc - Qq(x)'*g0);
P = P - K*S*K';

end

