function [ chi ] = sigmaPoints( x, P, dim )
%Generating sigma points using cubature rule
sqrt_P = chol(P)';
chi = zeros(dim,2*dim);
for i = 1:dim
    chi(:,i) = x + sqrt(dim)*sqrt_P(:,i);
end
for i = dim+1:2*dim
    chi(:,i) = x - sqrt(dim)*sqrt_P(:,i-dim);
end

