function [ summCov ] = motionCov( A, chi, x, dim )
%Calculate the covariance summation in prediction step
summCov = zeros(dim,dim);
for i = 1:2*dim
    summCov = summCov + (A*chi(:,i) - x)*(A*chi(:,i) - x)';
end

